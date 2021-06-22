use std::{error, fmt, io};

use pyo3::{exceptions::PyRuntimeError, PyErr};

use super::votable::DataType;

#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    Parse(&'static str),
    FieldType(usize, DataType, DataType),
    MissingField(Vec<u8>),
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

impl From<Error> for PyErr {
    fn from(e: Error) -> Self {
        match e {
            Error::Io(inner) => inner.into(),
            _ => PyRuntimeError::new_err(e.to_string()),
        }
    }
}
