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

use std::{borrow::Cow, error, fmt, io};

use pyo3::{exceptions::PyRuntimeError, PyErr};

use super::votable::DataType;

#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    Parse(&'static str),
    FieldType(usize, DataType, DataType),
    MissingField(Cow<'static, [u8]>),
    MissingId(String),
    Io(io::Error),
    Xml(quick_xml::Error),
    Other(Box<dyn error::Error + Send + Sync>),
}

impl Error {
    pub fn parse(message: &'static str) -> Self {
        Self::Parse(message)
    }

    pub fn field_type(ordinal: usize, expected: DataType, actual: DataType) -> Self {
        Self::FieldType(ordinal, expected, actual)
    }

    pub fn new(error: impl Into<Box<dyn error::Error + Send + Sync>>) -> Self {
        Self::Other(error.into())
    }

    pub fn missing_id(id_type: &str) -> Self {
        Self::MissingId(id_type.into())
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(s) => write!(f, "Parser failure: {}", s),
            Self::FieldType(o, e, a) => write!(
                f,
                "Field {} type mismatch (expected {}, actual {})",
                o, e, a
            ),
            Self::MissingField(s) => write!(f, "Missing field {}", String::from_utf8_lossy(s)),
            Self::MissingId(s) => write!(f, "Missing ID ({})", s),
            Self::Io(e) => write!(f, "Io error: {}", e),
            Self::Xml(e) => write!(f, "XML error: {}", e),
            Self::Other(e) => write!(f, "Error: {}", e),
        }
    }
}

impl error::Error for Error {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::Xml(e) => Some(e),
            Self::Other(e) => Some(e.as_ref()),
            _ => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

impl From<io::ErrorKind> for Error {
    fn from(e: io::ErrorKind) -> Self {
        Self::Io(e.into())
    }
}

impl From<quick_xml::Error> for Error {
    fn from(e: quick_xml::Error) -> Self {
        Self::Xml(e)
    }
}

impl From<globset::Error> for Error {
    fn from(e: globset::Error) -> Self {
        Self::new(e)
    }
}

impl From<Error> for PyErr {
    fn from(e: Error) -> Self {
        match e {
            Error::Io(inner) => inner.into(),
            _ => PyRuntimeError::new_err(e.to_string()),
        }
    }
}
