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

use std::{fmt, num::NonZeroUsize};

use super::error::AppError;

mod read;
mod write;

pub use read::{RecordAccessor, VotableReader};
pub use write::VotableWriter;

const VOTABLE_NS: &[u8] = b"http://www.ivoa.net/xml/VOTable/v1.3";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    Short,
    Int,
    Long,
    Float,
    Double,
    Char,
    String(Option<NonZeroUsize>),
}

impl DataType {
    fn parse_bytes(bstr: &[u8]) -> Result<Self, AppError> {
        match bstr {
            b"short" => Ok(Self::Short),
            b"int" => Ok(Self::Int),
            b"long" => Ok(Self::Long),
            b"float" => Ok(Self::Float),
            b"double" => Ok(Self::Double),
            b"char" => Ok(Self::Char),
            _ => Err(AppError::parse("Unsupported data type")),
        }
    }

    fn width(&self) -> Option<NonZeroUsize> {
        match self {
            Self::Short => NonZeroUsize::new(2),
            Self::Int => NonZeroUsize::new(4),
            Self::Long => NonZeroUsize::new(8),
            Self::Float => NonZeroUsize::new(4),
            Self::Double => NonZeroUsize::new(8),
            Self::Char => NonZeroUsize::new(1),
            Self::String(Some(s)) => Some(*s),
            Self::String(None) => None,
        }
    }
}

impl fmt::Display for DataType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Short => f.write_str("short"),
            Self::Int => f.write_str("int"),
            Self::Long => f.write_str("long"),
            Self::Float => f.write_str("float"),
            Self::Double => f.write_str("double"),
            Self::Char => f.write_str("char"),
            Self::String(_) => f.write_str("char array"),
        }
    }
}

enum FieldGetter<R> {
    Short(fn(&R) -> Option<i16>),
    Long(fn(&R) -> Option<i64>),
    Float(fn(&R) -> f32),
    Double(fn(&R) -> f64),
}

pub struct FieldInfo<R> {
    name: &'static str,
    unit: Option<&'static str>,
    ucd: Option<&'static str>,
    description: &'static str,
    getter: FieldGetter<R>,
}

impl<R> FieldInfo<R> {
    pub fn short(
        name: &'static str,
        unit: Option<&'static str>,
        ucd: Option<&'static str>,
        description: &'static str,
        getter: fn(&R) -> Option<i16>,
    ) -> Self {
        Self {
            name,
            unit,
            ucd,
            description,
            getter: FieldGetter::Short(getter),
        }
    }

    pub fn long(
        name: &'static str,
        unit: Option<&'static str>,
        ucd: Option<&'static str>,
        description: &'static str,
        getter: fn(&R) -> Option<i64>,
    ) -> Self {
        Self {
            name,
            unit,
            ucd,
            description,
            getter: FieldGetter::Long(getter),
        }
    }

    pub fn float(
        name: &'static str,
        unit: Option<&'static str>,
        ucd: Option<&'static str>,
        description: &'static str,
        getter: fn(&R) -> f32,
    ) -> Self {
        Self {
            name,
            unit,
            ucd,
            description,
            getter: FieldGetter::Float(getter),
        }
    }

    pub fn double(
        name: &'static str,
        unit: Option<&'static str>,
        ucd: Option<&'static str>,
        description: &'static str,
        getter: fn(&R) -> f64,
    ) -> Self {
        Self {
            name,
            unit,
            ucd,
            description,
            getter: FieldGetter::Double(getter),
        }
    }
}
