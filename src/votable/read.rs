/*
* gaia-stardb: Processing Gaia DR2 for celestia.Sci/Celestia
* Copyright (C) 2019–2021  Andrew Tribick
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

use std::{
    cmp,
    collections::HashMap,
    convert::TryInto,
    io::{self, BufRead, BufReader, ErrorKind, Read},
    mem,
};

use arrayvec::ArrayVec;
use bitvec::prelude::*;
use byteorder::{BigEndian, ReadBytesExt};
use flate2::read::GzDecoder;
use quick_xml::{
    events::{attributes::Attributes, Event},
    Reader as XmlReader,
};

use super::{DataType, VOTABLE_NS};
use crate::error::AppError;

fn parse_field(attributes: Attributes) -> Result<(Vec<u8>, DataType), AppError> {
    let mut name = None;
    let mut datatype = None;
    let mut is_variable_length_array = false;
    for attribute_result in attributes {
        let attribute = attribute_result?;
        match attribute.key {
            b"name" => name = Some(attribute.value.into_owned()),
            b"datatype" => datatype = Some(DataType::parse_bytes(&attribute.value)?),
            b"arraysize" => {
                if attribute.value.as_ref() == b"*" {
                    is_variable_length_array = true;
                } else {
                    return Err(AppError::parse("Fixed size arrays not supported"));
                }
            }
            _ => (),
        }
    }

    match (name, datatype) {
        (Some(n), Some(DataType::Char)) if is_variable_length_array => Ok((n, DataType::Char)),
        (Some(_), Some(DataType::Char)) => Err(AppError::parse("Char fields not supported")),
        (Some(_), Some(_)) if is_variable_length_array => {
            Err(AppError::parse("Non-string arrays not supported"))
        }
        (Some(n), Some(dt)) => Ok((n, dt)),
        _ => Err(AppError::parse("Field must have name and datatype")),
    }
}

pub struct VotableReader<R: Read> {
    reader: Binary2Reader<BufReader<GzDecoder<R>>>,
    field_types: Vec<DataType>,
    field_offsets: Vec<usize>,
    field_names: HashMap<Vec<u8>, usize>,
    mask_width: usize,
    has_dynamic_lengths: bool,
    buffer: Vec<u8>,
}

impl<R: Read> VotableReader<R> {
    pub fn new(source: R) -> Result<Self, AppError> {
        let decoder = GzDecoder::new(source);
        let buf_reader = BufReader::new(decoder);
        let mut xml_reader = XmlReader::from_reader(buf_reader);
        let mut buf = Vec::with_capacity(1024);
        let mut ns_buf = Vec::with_capacity(256);
        let mut field_names = HashMap::new();
        let mut field_types = Vec::new();
        let mut has_dynamic_lengths = false;
        let mut field_offsets = Vec::new();
        let mut offset = 0;
        let mut is_binary2 = false;
        loop {
            match xml_reader.read_namespaced_event(&mut buf, &mut ns_buf)? {
                (Some(VOTABLE_NS), Event::Start(ref e))
                | (Some(VOTABLE_NS), Event::Empty(ref e)) => match e.name() {
                    b"FIELD" => {
                        let (name, data_type) = parse_field(e.attributes())?;
                        field_names.insert(name, field_names.len());
                        field_offsets.push(offset);
                        match data_type.width() {
                            Some(w) => offset += w.get(),
                            None => {
                                has_dynamic_lengths = true;
                                offset += 4;
                            }
                        }
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
            return Err(AppError::parse("No fields found in file"));
        }

        if !is_binary2 {
            return Err(AppError::parse("Format not supported"));
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
            has_dynamic_lengths,
            buffer: vec![0; mask_width + offset],
        })
    }

    pub fn ordinal(&self, name: &[u8]) -> Result<usize, AppError> {
        self.field_names
            .get(name)
            .copied()
            .ok_or_else(|| AppError::MissingField(name.to_owned().into()))
    }

    pub fn read(&mut self) -> Result<Option<RecordAccessor>, AppError> {
        let has_record = if self.has_dynamic_lengths {
            self.read_dynamic()?
        } else {
            self.read_fixed()?
        };

        let result = if has_record {
            Some(RecordAccessor {
                mask: (&self.buffer[..self.mask_width]).view_bits(),
                field_types: &self.field_types,
                field_offsets: &self.field_offsets,
                data: &self.buffer[self.mask_width..],
            })
        } else {
            None
        };

        Ok(result)
    }

    fn read_fixed(&mut self) -> Result<bool, AppError> {
        let mut length = 0;
        while length < self.buffer.len() {
            length += match self.reader.read(&mut self.buffer[length..]) {
                Ok(0) if length == 0 => return Ok(false),
                Ok(0) => return Err(AppError::Io(ErrorKind::UnexpectedEof.into())),
                Ok(n) => n,
                Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                Err(e) => return Err(e.into()),
            };
        }

        Ok(true)
    }

    fn read_dynamic(&mut self) -> Result<bool, AppError> {
        let mut length = 0;
        self.buffer.resize(self.mask_width, 0);
        while length < self.mask_width {
            length += match self.reader.read(&mut self.buffer[length..]) {
                Ok(0) if length == 0 => return Ok(false),
                Ok(0) => return Err(AppError::Io(ErrorKind::UnexpectedEof.into())),
                Ok(n) => n,
                Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                Err(e) => return Err(e.into()),
            };
        }

        let mut position = self.buffer.len();
        for (field_type, field_offset) in self.field_types.iter().zip(self.field_offsets.iter_mut())
        {
            *field_offset = position - self.mask_width;
            match field_type.width() {
                Some(w) => {
                    self.buffer.resize(position + w.get(), 0);
                    self.reader.read_exact(&mut self.buffer[position..])?;
                    position += w.get();
                }
                None => {
                    self.buffer.resize(position + mem::size_of::<u32>(), 0);
                    self.reader.read_exact(&mut self.buffer[position..])?;
                    let length = (&self.buffer[position..]).read_u32::<BigEndian>()? as usize;
                    position += mem::size_of::<u32>();
                    if length > 0 {
                        self.buffer.resize(position + length, 0);
                        self.reader.read_exact(&mut self.buffer[position..])?;
                        position += length;
                    }
                }
            }
        }

        Ok(true)
    }
}

pub struct RecordAccessor<'a> {
    mask: &'a BitSlice<Msb0, u8>,
    field_types: &'a [DataType],
    field_offsets: &'a [usize],
    data: &'a [u8],
}

impl<'a> RecordAccessor<'a> {
    pub fn read_i16(&self, ordinal: usize) -> Result<Option<i16>, AppError> {
        let field_type = self.field_types[ordinal];
        if field_type != DataType::Short {
            return Err(AppError::field_type(ordinal, DataType::Short, field_type));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some(
            (&self.data[offset..offset + mem::size_of::<i16>()]).read_i16::<BigEndian>()?,
        ))
    }

    pub fn read_i32(&self, ordinal: usize) -> Result<Option<i32>, AppError> {
        let field_type = self.field_types[ordinal];
        if field_type != DataType::Int {
            return Err(AppError::field_type(ordinal, DataType::Int, field_type));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some(
            (&self.data[offset..offset + mem::size_of::<i32>()]).read_i32::<BigEndian>()?,
        ))
    }

    pub fn read_i64(&self, ordinal: usize) -> Result<Option<i64>, AppError> {
        let field_type = self.field_types[ordinal];
        if field_type != DataType::Long {
            return Err(AppError::field_type(ordinal, DataType::Long, field_type));
        }

        if self.mask[ordinal] {
            return Ok(None);
        }

        let offset = self.field_offsets[ordinal];
        Ok(Some(
            (&self.data[offset..offset + mem::size_of::<i64>()]).read_i64::<BigEndian>()?,
        ))
    }

    pub fn read_f32(&self, ordinal: usize) -> Result<f32, AppError> {
        let field_type = self.field_types[ordinal];
        if field_type != DataType::Float {
            return Err(AppError::field_type(ordinal, DataType::Float, field_type));
        }

        if self.mask[ordinal] {
            return Ok(f32::NAN);
        }

        let offset = self.field_offsets[ordinal];
        Ok((&self.data[offset..offset + mem::size_of::<f32>()]).read_f32::<BigEndian>()?)
    }

    pub fn read_f64(&self, ordinal: usize) -> Result<f64, AppError> {
        let field_type = self.field_types[ordinal];
        if field_type != DataType::Double {
            return Err(AppError::field_type(ordinal, DataType::Double, field_type));
        }

        if self.mask[ordinal] {
            return Ok(f64::NAN);
        }

        let offset = self.field_offsets[ordinal];
        Ok((&self.data[offset..offset + mem::size_of::<f64>()]).read_f64::<BigEndian>()?)
    }

    pub fn read_char<const CAP: usize>(
        &self,
        ordinal: usize,
    ) -> Result<ArrayVec<u8, CAP>, AppError> {
        let field_type = self.field_types[ordinal];
        if field_type != DataType::Char {
            return Err(AppError::field_type(ordinal, DataType::Char, field_type));
        }

        if self.mask[ordinal] {
            return Ok(ArrayVec::new());
        }

        let offset = self.field_offsets[ordinal];
        let data_offset = offset + mem::size_of::<u32>();
        let length = (&self.data[offset..data_offset]).read_u32::<BigEndian>()? as usize;
        Ok(self.data[data_offset..data_offset + length].try_into()?)
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

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn binary2_read() {
        let source: &[u8] = b"AgMFBwsNERMXHR8l\nKSsvNTs9Q0c=\n</STREAM>";
        let expected = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        ];
        let buf_reader = BufReader::with_capacity(5, source);
        let mut reader = Binary2Reader::new(buf_reader);
        let mut actual = Vec::new();
        let length = reader.read_to_end(&mut actual).unwrap();
        assert_eq!(length, 20);
        assert_eq!(&expected, actual.as_slice());
    }
}
