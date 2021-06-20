use std::io::{self, BufRead, BufReader, ErrorKind, Read};

use base64::read::DecoderReader;
use bitvec::vec::BitVec;
use byteorder::{BigEndian, ReadBytesExt};
use flate2::read::GzDecoder;
use quick_xml::{events::Event, Reader as XmlReader};

use crate::error::Error;

const VOTABLE_NS: &[u8] = b"http://www.ivoa.net/xml/VOTable/v1.3";

pub fn parse_votable(source: impl Read, hip_ids: &mut BitVec) -> Result<(), Error> {
    let decoder = GzDecoder::new(source);
    let buf_reader = BufReader::new(decoder);
    let mut reader = XmlReader::from_reader(buf_reader);
    let mut xml_buf = Vec::with_capacity(1024);
    let mut ns_buf = Vec::with_capacity(256);
    let mut record_count = 0;
    let mut record_size = 0;
    let mut hip_offset = 0;
    loop {
        match reader.read_namespaced_event(&mut xml_buf, &mut ns_buf)? {
            (Some(VOTABLE_NS), Event::Start(ref e)) => match e.name() {
                b"FIELD" => {
                    let mut name = None;
                    let mut datatype = None;
                    for attribute_result in e.attributes() {
                        let attribute = attribute_result?;
                        match attribute.key {
                            b"name" => name = Some(attribute.value),
                            b"datatype" => datatype = Some(attribute.value),
                            b"arraysize" => unimplemented!(),
                            _ => (),
                        }
                    }

                    record_count += 1;

                    if matches!(name.as_deref(), Some(b"hip")) {
                        hip_offset = record_size;
                    }

                    match datatype.as_deref() {
                        Some(b"int") => record_size += 4,
                        Some(b"long") => record_size += 8,
                        Some(b"float") => record_size += 4,
                        Some(b"double") => record_size += 8,
                        Some(b"short") => record_size += 2,
                        _ => unimplemented!(),
                    }
                }

                b"STREAM" => {
                    let mask_size = (record_count + 7) / 8;
                    let mut buffer = vec![0; record_size + mask_size];
                    let bytes_reader = reader.into_underlying_reader();
                    let mut lines_reader = Base64DataReader::new(bytes_reader);
                    let mut decoder = DecoderReader::new(&mut lines_reader, base64::STANDARD);
                    'outer: loop {
                        let mut p = 0;
                        while p < buffer.len() {
                            p += match decoder.read(&mut buffer[p..]) {
                                Ok(0) if p == 0 => break 'outer,
                                Ok(0) => {
                                    return Err(Error::IoError(ErrorKind::UnexpectedEof.into()))
                                }
                                Ok(n) => n,
                                Err(ref e) if e.kind() == ErrorKind::Interrupted => continue,
                                Err(e) => return Err(e.into()),
                            }
                        }
                        let hip = (&buffer[hip_offset + mask_size..hip_offset + mask_size + 4])
                            .read_i32::<BigEndian>()? as usize;
                        hip_ids.set(hip - 1, true);
                    }

                    break;
                }

                _ => (),
            },
            (_, Event::Eof) => break,
            _ => (),
        }

        xml_buf.clear();
    }

    Ok(())
}

struct Base64DataReader<R: BufRead> {
    inner: R,
}

impl<R: BufRead> Base64DataReader<R> {
    pub fn new(inner: R) -> Self {
        Self { inner }
    }
}

impl<R: BufRead> Read for Base64DataReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        loop {
            let src = self.inner.fill_buf()?;
            if src.is_empty() {
                return Ok(0);
            }

            let close_pos = memchr::memchr(b'<', src).unwrap_or(src.len());
            if close_pos == 0 {
                return Ok(0);
            }

            let mut start_pos = 0;
            let mut end_pos = close_pos;
            for p in memchr::memchr3_iter(b' ', b'\n', b'\t', &src[..close_pos]) {
                if p == start_pos {
                    start_pos += 1;
                } else {
                    end_pos = p;
                    break;
                }
            }

            if start_pos >= end_pos {
                self.inner.consume(end_pos);
                continue;
            }

            let length = std::cmp::min(buf.len(), end_pos - start_pos);
            end_pos = start_pos + length;
            buf[..length].copy_from_slice(&src[start_pos..end_pos]);
            self.inner.consume(end_pos);

            return Ok(length);
        }
    }
}
