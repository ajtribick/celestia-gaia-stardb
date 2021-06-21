use std::{
    borrow::Cow,
    cmp,
    collections::HashMap,
    io::{self, BufRead, BufReader, ErrorKind, Read},
};

use bitvec::prelude::*;
use byteorder::{BigEndian, ReadBytesExt};
use flate2::read::GzDecoder;
use quick_xml::{
    events::{attributes::Attributes, Event},
    Reader as XmlReader,
};

use super::error::Error;

const VOTABLE_NS: &[u8] = b"http://www.ivoa.net/xml/VOTable/v1.3";

enum DataType {
    Short,
    Int,
    Long,
    Float,
    Double,
}

impl DataType {
    fn parse_bytes(bstr: &[u8]) -> Result<Self, Error> {
        match bstr {
            b"short" => Ok(Self::Short),
            b"int" => Ok(Self::Int),
            b"long" => Ok(Self::Long),
            b"float" => Ok(Self::Float),
            b"double" => Ok(Self::Double),
            _ => Err(Error::parse("Unsupported data type")),
        }
    }

    fn width(&self) -> usize {
        match self {
            Self::Short => 2,
            Self::Int => 4,
            Self::Long => 8,
            Self::Float => 4,
            Self::Double => 8,
        }
    }
}

fn parse_field(attributes: Attributes) -> Result<(Vec<u8>, DataType), Error> {
    let mut name = None;
    let mut datatype = None;
    for attribute_result in attributes {
        let attribute = attribute_result?;
        match attribute.key {
            b"name" => name = Some(attribute.value.into_owned()),
            b"datatype" => datatype = Some(DataType::parse_bytes(&attribute.value)?),
            b"arraysize" => {
                return Err(Error::parse("Array types not supported"))
            }
            _ => (),
        }
    }

    match (name, datatype) {
        (Some(n), Some(dt)) => Ok((n, dt)),
        _ => Err(Error::parse("Field missing name and datatype")),
    }
}

pub struct VotableReader<R: Read> {
    reader: Binary2Reader<BufReader<GzDecoder<R>>>,
    field_types: Vec<DataType>,
    field_offsets: Vec<usize>,
    field_names: HashMap<Vec<u8>, usize>,
    mask_width: usize,
    buffer: Vec<u8>,
}

impl<R: Read> VotableReader<R> {
    pub fn new(source: R) -> Result<Self, Error> {
        let decoder = GzDecoder::new(source);
        let buf_reader = BufReader::new(decoder);
        let mut xml_reader = XmlReader::from_reader(buf_reader);
        let mut buf = Vec::with_capacity(1024);
        let mut ns_buf = Vec::with_capacity(256);
        let mut field_names = HashMap::new();
        let mut field_types = Vec::new();
        let mut field_offsets = Vec::new();
        let mut offset = 0;
        let mut is_binary2 = false;
        loop {
            match xml_reader.read_namespaced_event(&mut buf, &mut ns_buf)? {
                (Some(VOTABLE_NS), Event::Start(ref e)) => match e.name() {
                    b"FIELD" => {
                        let (name, data_type) = parse_field(e.attributes())?;
                        field_names.insert(name, field_names.len());
                        field_offsets.push(offset);
                        offset += data_type.width();
                        field_types.push(data_type);
                    }
                    b"BINARY2" => is_binary2 = true,
                    b"STREAM" => break,
                    _ => (),
                },
                (_, Event::Eof) => return Err(ErrorKind::UnexpectedEof.into()),
                _ => (),
            }
        }

        if field_offsets.len() == 0 {
            return Err(Error::parse("No fields found in file"));
        }

        if !is_binary2 {
            return Err(Error::parse("Format not supported"));
        }

        let mask_width = (field_offsets.len() + 7) / 8;

        let buf_reader = xml_reader.into_underlying_reader();
        let reader = Binary2Reader::new(buf_reader);

        Ok(Self {
            reader,
            field_types,
            field_offsets,
            field_names,
            mask_width,
            buffer: vec![0; mask_width + offset],
        })
    }

    pub fn ordinal(&self, name: &[u8]) -> Option<usize> {
        self.field_names.get(name).copied()
    }

    pub fn read(&mut self) -> Result<Option<RecordAccessor>, Error> {
        let mut length = 0;
        while length < self.buffer.len() {
            length += match self.reader.read(&mut self.buffer[length..]) {
                Ok(0) if length == 0 => return Ok(None),
                Ok(0) => return Err(Error::Io(ErrorKind::UnexpectedEof.into())),
                Ok(n) => n,
                Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                Err(e) => return Err(e.into()),
            };
        }

        Ok(Some(RecordAccessor {
            mask: (&self.buffer[..self.mask_width]).view_bits(),
            field_types: &self.field_types,
            field_offsets: Cow::from(&self.field_offsets),
            data: &self.buffer[self.mask_width..],
        }))
    }
}

pub struct RecordAccessor<'a> {
    mask: &'a BitSlice<Msb0, u8>,
    field_types: &'a [DataType],
    field_offsets: Cow<'a, [usize]>,
    data: &'a [u8],
}

impl<'a> RecordAccessor<'a> {
    pub fn read_i16(&self, ordinal: usize) -> Result<Option<i16>, Error> {
        if !matches!(self.field_types[ordinal], DataType::Short) {
            return Err(Error::field_type("Field type mismatch"));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some((&self.data[offset..offset + std::mem::size_of::<i16>()])
            .read_i16::<BigEndian>()?))
    }

    pub fn read_i32(&self, ordinal: usize) -> Result<Option<i32>, Error> {
        if !matches!(self.field_types[ordinal], DataType::Int) {
            return Err(Error::field_type("Field type mismatch"));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some((&self.data[offset..offset + std::mem::size_of::<i32>()])
            .read_i32::<BigEndian>()?))
    }

    pub fn read_i64(&self, ordinal: usize) -> Result<Option<i64>, Error> {
        if !matches!(self.field_types[ordinal], DataType::Long) {
            return Err(Error::field_type("Field type mismatch"));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some((&self.data[offset..offset + std::mem::size_of::<i64>()])
            .read_i64::<BigEndian>()?))
    }

    pub fn read_f32(&self, ordinal: usize) -> Result<Option<f32>, Error> {
        if !matches!(self.field_types[ordinal], DataType::Float) {
            return Err(Error::field_type("Field type mismatch"));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some((&self.data[offset..offset + std::mem::size_of::<f32>()])
            .read_f32::<BigEndian>()?))
    }

    pub fn read_f64(&self, ordinal: usize) -> Result<Option<f64>, Error> {
        if !matches!(self.field_types[ordinal], DataType::Float) {
            return Err(Error::field_type("Field type mismatch"));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some((&self.data[offset..offset + std::mem::size_of::<f64>()])
            .read_f64::<BigEndian>()?))
    }
}

struct Binary2Reader<R: BufRead> {
    reader: R,
    buffer: Vec<u8>,
    line_buffer: Vec<u8>,
    position: usize,
}

impl<R: BufRead> Binary2Reader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buffer: Vec::with_capacity(64),
            line_buffer: Vec::with_capacity(64),
            position: 0,
        }
    }
}

impl<R: BufRead> Read for Binary2Reader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let data = self.fill_buf()?;
        let length = cmp::min(buf.len(), data.len());
        buf[..length].copy_from_slice(&data[..length]);
        self.consume(length);
        Ok(length)
    }
}

impl<R: BufRead> BufRead for Binary2Reader<R> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.position < self.buffer.len() {
            return Ok(&self.buffer[self.position..]);
        }

        self.buffer.clear();
        self.line_buffer.clear();
        self.position = 0;
        loop {
            let data = self.reader.fill_buf()?;
            let mut start_pos = 0;
            let mut end_pos = memchr::memchr(b'<', data).unwrap_or(data.len());
            if end_pos == 0 {
                return Ok(&self.buffer);
            }

            for p in memchr::memchr3_iter(b' ', b'\n', b'\t', &data[..end_pos]) {
                if p == start_pos {
                    start_pos += 1;
                } else {
                    end_pos = p;
                    break;
                }
            }

            if start_pos >= end_pos {
                self.reader.consume(end_pos);
                continue;
            }

            self.line_buffer
                .extend_from_slice(&data[start_pos..end_pos]);
            if end_pos == data.len() && (self.line_buffer.len() & 0b11) != 0 {
                self.reader.consume(end_pos);
                continue;
            }

            base64::decode_config_buf(&self.line_buffer, base64::STANDARD, &mut self.buffer)
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;
            self.reader.consume(end_pos);
            return Ok(&self.buffer);
        }
    }

    fn consume(&mut self, amt: usize) {
        self.position += amt;
    }
}
