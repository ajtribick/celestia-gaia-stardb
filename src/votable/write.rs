use std::{
    cmp,
    convert::TryInto,
    io::{self, ErrorKind, Write},
    mem,
};

use bitvec::prelude::*;
use byteorder::{BigEndian, WriteBytesExt};
use quick_xml::{
    events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event},
    Writer,
};

use super::{FieldGetter, VotableRecord};
use crate::error::Error;

pub struct VotableWriter<W: Write> {
    writer: Option<Writer<W>>,
    data_writer: Option<Base64Writer<W>>,
    mask_data: BitVec<Msb0, u8>,
    field_data: Vec<u8>,
    field_count: usize,
    min_length: usize,
}

impl<W: Write> VotableWriter<W> {
    pub fn new(inner: W) -> Result<Self, Error> {
        let mut writer = Writer::new(inner);
        writer.write_event(Event::Decl(BytesDecl::new(b"1.0", Some(b"UTF-8"), None)))?;
        writer.write(b"\n")?;
        let doc_element = BytesStart::borrowed(b"VOTABLE", b"VOTABLE".len()).with_attributes([
            (b"version".as_ref(), b"1.3".as_ref()),
            (b"xmlns", b"http://www.ivoa.net/xml/VOTable/v1.3"),
            (b"xmlns:xsi", b"http://www.w3.org/2001/XMLSchema-instance"),
            (
                b"xsi:schemaLocation",
                b"http://www.ivoa.net/xml/VOTable/v1.3 http://www.ivoa.net/xml/VOTable/v1.3",
            ),
        ]);
        writer.write_event(Event::Start(doc_element))?;
        writer.write(b"\n")?;
        writer.write_event(Event::Start(BytesStart::borrowed(
            br#"RESOURCE type="results""#,
            b"RESOURCE".len(),
        )))?;
        writer.write(b"\n")?;
        writer.write_event(Event::Start(BytesStart::borrowed(b"TABLE", b"TABLE".len())))?;
        writer.write(b"\n")?;
        Ok(Self {
            writer: Some(writer),
            data_writer: None,
            mask_data: BitVec::new(),
            field_data: Vec::new(),
            field_count: 0,
            min_length: 0,
        })
    }

    pub fn write_fields<T: VotableRecord>(&mut self) -> Result<(), Error> {
        let writer = self
            .writer
            .as_mut()
            .expect("Cannot add fields after data is written");
        for field in T::fields() {
            let (data_type, array_size, data_size): (&[u8], Option<&str>, usize) =
                match field.getter {
                    FieldGetter::Short(_) => (b"short", None, mem::size_of::<i16>()),
                    FieldGetter::Int(_) => (b"int", None, mem::size_of::<i32>()),
                    FieldGetter::Long(_) => (b"long", None, mem::size_of::<i64>()),
                    FieldGetter::Float(_) => (b"float", None, mem::size_of::<f32>()),
                    FieldGetter::Double(_) => (b"double", None, mem::size_of::<f64>()),
                    FieldGetter::Char(_) => (b"char", Some("*"), mem::size_of::<u32>()),
                };

            self.field_count += 1;
            self.min_length += data_size;

            let mut element = BytesStart::borrowed(b"FIELD", b"FIELD".len()).with_attributes([
                (b"name".as_ref(), field.name.as_bytes()),
                (b"datatype", data_type),
            ]);

            if let Some(a) = array_size {
                element.push_attribute(("arraysize", a));
            }

            if let Some(u) = field.ucd {
                element.push_attribute(("ucd", u));
            }

            if let Some(u) = field.unit {
                element.push_attribute(("unit", u));
            }

            writer.write_event(Event::Start(element))?;
            writer.write(b"\n")?;

            writer.write_event(Event::Start(BytesStart::borrowed(
                b"DESCRIPTION",
                b"DESCRIPTION".len(),
            )))?;
            writer.write_event(Event::Text(BytesText::from_plain_str(field.description)))?;
            writer.write_event(Event::End(BytesEnd::borrowed(b"DESCRIPTION")))?;
            writer.write(b"\n")?;

            writer.write_event(Event::End(BytesEnd::borrowed(b"FIELD")))?;
            writer.write(b"\n")?;
        }

        Ok(())
    }

    pub fn add_data<T: VotableRecord>(&mut self, data: &T) -> Result<(), Error> {
        for field in T::fields() {
            match field.getter {
                FieldGetter::Short(f) => {
                    let (value, mask) = f(data).map_or((Default::default(), true), |x| (x, false));
                    self.field_data.write_i16::<BigEndian>(value)?;
                    self.mask_data.push(mask);
                }
                FieldGetter::Int(f) => {
                    let (value, mask) = f(data).map_or((Default::default(), true), |x| (x, false));
                    self.field_data.write_i32::<BigEndian>(value)?;
                    self.mask_data.push(mask);
                }
                FieldGetter::Long(f) => {
                    let (value, mask) = f(data).map_or((Default::default(), true), |x| (x, false));
                    self.field_data.write_i64::<BigEndian>(value)?;
                    self.mask_data.push(mask);
                }
                FieldGetter::Float(f) => {
                    let value = f(data);
                    self.field_data.write_f32::<BigEndian>(value)?;
                    self.mask_data.push(value.is_nan());
                }
                FieldGetter::Double(f) => {
                    let value = f(data);
                    self.field_data.write_f64::<BigEndian>(value)?;
                    self.mask_data.push(value.is_nan());
                }
                FieldGetter::Char(f) => match f(data) {
                    Some(s) => {
                        let length = s.len().try_into().unwrap_or(u32::MAX);
                        self.field_data.write_u32::<BigEndian>(length)?;
                        self.field_data.extend_from_slice(&s[..length as usize]);
                        self.mask_data.push(length == 0);
                    }
                    None => {
                        self.field_data.write_u32::<BigEndian>(0)?;
                        self.mask_data.push(true);
                    }
                },
            }
        }

        Ok(())
    }

    pub fn write_row(&mut self) -> Result<(), Error> {
        assert_eq!(
            self.mask_data.len(),
            self.field_count,
            "Number of fields written does not match field count"
        );
        let writer = match self.data_writer.as_mut() {
            Some(w) => w,
            None => {
                let mut writer = self.writer.take().expect("Invalid state");

                self.mask_data.reserve(self.field_count);
                self.field_data.reserve(self.min_length);

                writer.write_event(Event::Start(BytesStart::borrowed(b"DATA", b"DATA".len())))?;
                writer.write(b"\n")?;
                writer.write_event(Event::Start(BytesStart::borrowed(
                    b"BINARY2",
                    b"BINARY2".len(),
                )))?;
                writer.write(b"\n")?;
                writer.write_event(Event::Start(BytesStart::borrowed(
                    b"STREAM encoding='base64'",
                    b"STREAM".len(),
                )))?;
                writer.write(b"\n")?;
                self.data_writer.insert(Base64Writer::new(writer))
            }
        };

        writer.write_all(self.mask_data.as_raw_slice())?;
        writer.write_all(&self.field_data)?;

        self.field_data.clear();
        self.mask_data.clear();
        Ok(())
    }

    pub fn finish(self) -> Result<W, Error> {
        let mut writer = match self.data_writer {
            Some(data_writer) => {
                let mut writer = data_writer.finish()?;
                writer.write(b"\n")?;
                writer.write_event(Event::End(BytesEnd::borrowed(b"STREAM")))?;
                writer.write(b"\n")?;
                writer.write_event(Event::End(BytesEnd::borrowed(b"BINARY2")))?;
                writer.write(b"\n")?;
                writer.write_event(Event::End(BytesEnd::borrowed(b"DATA")))?;
                writer.write(b"\n")?;
                writer
            }
            None => self.writer.expect("Invalid state"),
        };

        writer.write_event(Event::End(BytesEnd::borrowed(b"TABLE")))?;
        writer.write(b"\n")?;
        writer.write_event(Event::End(BytesEnd::borrowed(b"RESOURCE")))?;
        writer.write(b"\n")?;
        writer.write_event(Event::End(BytesEnd::borrowed(b"VOTABLE")))?;
        writer.write(b"\n")?;

        Ok(writer.into_inner())
    }
}

struct Base64Writer<W: Write> {
    writer: Writer<W>,
    buffer: Vec<u8>,
    output_buffer: Vec<u8>,
    position: usize,
}

impl<W: Write> Base64Writer<W> {
    pub fn new(writer: Writer<W>) -> Self {
        let mut output_buffer = vec![0; 65];
        output_buffer[64] = b'\n';
        Self {
            writer,
            buffer: vec![0; 48],
            output_buffer,
            position: 0,
        }
    }

    pub fn finish(mut self) -> Result<Writer<W>, Error> {
        if self.position == 0 {
            return Ok(self.writer);
        }
        let output_length = base64::encode_config_slice(
            &self.buffer[..self.position],
            base64::STANDARD,
            &mut self.output_buffer[..64],
        );
        self.output_buffer[output_length] = b'\n';
        self.writer
            .write(&self.output_buffer[..output_length + 1])?;
        self.writer.write(b"\n")?;
        Ok(self.writer)
    }
}

impl<W: Write> Write for Base64Writer<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let length = cmp::min(buf.len(), self.buffer.len() - self.position);
        let new_position = self.position + length;
        self.buffer[self.position..new_position].copy_from_slice(&buf[..length]);
        self.position = new_position;
        if self.position == self.buffer.len() {
            base64::encode_config_slice(
                &self.buffer,
                base64::STANDARD,
                &mut self.output_buffer[..64],
            );
            self.writer
                .write(&self.output_buffer)
                .map_err(|e| match e {
                    quick_xml::Error::Io(e) => e,
                    e => io::Error::new(ErrorKind::InvalidData, e),
                })?;
            self.position = 0;
        }

        Ok(length)
    }

    fn flush(&mut self) -> io::Result<()> {
        let length = self.position - self.position % 3;
        if length == 0 {
            return Ok(());
        }

        let output_length = base64::encode_config_slice(
            &self.buffer[..length],
            base64::STANDARD,
            &mut self.output_buffer[..64],
        );
        self.writer
            .write(&self.output_buffer[..output_length])
            .map_err(|e| match e {
                quick_xml::Error::Io(e) => e,
                e => io::Error::new(ErrorKind::InvalidData, e),
            })?;

        self.buffer.copy_within(length..self.position, 0);
        self.position -= length;

        Ok(())
    }
}
