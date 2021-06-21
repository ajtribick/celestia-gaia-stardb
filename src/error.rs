use std::{error, fmt, io};

use pyo3::{exceptions::PyRuntimeError, PyErr};

#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    Parse(&'static str),
    FieldType(&'static str),
    Io(io::Error),
    Xml(quick_xml::Error),
    Other(Box<dyn error::Error + Send + Sync>),
}

impl Error {
    pub(crate) fn parse(message: &'static str) -> Self {
        Self::Parse(message)
    }

    pub(crate) fn field_type(message: &'static str) -> Self {
        Self::FieldType(message)
    }

    pub fn new(error: impl Into<Box<dyn error::Error + Send + Sync>>) -> Self {
        Self::Other(error.into())
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(s) => write!(f, "Parser failure: {}", s),
            Self::FieldType(s) => write!(f, "Field type mismatch: {}", s),
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
            Error::Parse(s) => PyRuntimeError::new_err(s),
            Error::FieldType(s) => PyRuntimeError::new_err(s),
            Error::Io(e) => e.into(),
            Error::Xml(e) => PyRuntimeError::new_err(e.to_string()),
            Error::Other(e) => PyRuntimeError::new_err(e.to_string()),
        }
    }
}
