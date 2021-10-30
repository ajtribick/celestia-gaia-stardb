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
    collections::HashMap,
    io::{self, BufRead},
    ops::Range,
};

pub struct CsvReader<B: BufRead> {
    reader: B,
    fields: HashMap<String, usize>,
    line: String,
    offsets: Vec<Range<usize>>,
}

impl<B: BufRead> CsvReader<B> {
    pub fn new(mut reader: B) -> io::Result<Self> {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            return Err(io::ErrorKind::UnexpectedEof.into());
        }

        let fields: HashMap<_, _> = line
            .split(',')
            .enumerate()
            .map(|(i, s)| (s.trim().to_owned(), i))
            .collect();
        let offsets = Vec::with_capacity(fields.len());

        line.clear();

        Ok(Self {
            reader,
            fields,
            line,
            offsets,
        })
    }

    pub fn index(&self, name: &str) -> Option<usize> {
        self.fields.get(name).copied()
    }

    pub fn next(&mut self) -> io::Result<Option<()>> {
        self.line.clear();
        if self.reader.read_line(&mut self.line)? == 0 {
            return Ok(None);
        }

        self.offsets.clear();
        let mut pos = 0;
        while let Some(delimiter_pos) = self.line[pos..].find(',') {
            self.offsets.push(pos..pos + delimiter_pos);
            pos = pos + delimiter_pos + 1;
        }

        self.offsets.push(pos..self.line.len());
        if self.offsets.len() != self.fields.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Field count mismatch",
            ));
        }

        Ok(Some(()))
    }

    pub fn field(&self, index: usize) -> &str {
        self.line[self.offsets[index].clone()].trim()
    }
}
