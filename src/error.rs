use std::{error, fmt, io};

use pyo3::{exceptions::PyRuntimeError, PyErr};

#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    GlobError(globset::Error),
    IoError(io::Error),
    XmlError(quick_xml::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GlobError(e) => write!(f, "GlobError({})", e),
            Self::IoError(e) => write!(f, "IoError({})", e),
            Self::XmlError(e) => write!(f, "XmlError({})", e),
        }
    }
}

impl error::Error for Error {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::GlobError(e) => Some(e),
            Self::IoError(e) => Some(e),
            Self::XmlError(e) => Some(e),
        }
    }
}

impl From<globset::Error> for Error {
    fn from(e: globset::Error) -> Self {
        Self::GlobError(e)
    }
}

impl From<io::Error> for Error {
    fn from(e: io::Error) -> Self {
        Self::IoError(e)
    }
}

impl From<io::ErrorKind> for Error {
    fn from(e: io::ErrorKind) -> Self {
        Self::IoError(e.into())
    }
}

impl From<quick_xml::Error> for Error {
    fn from(e: quick_xml::Error) -> Self {
        Self::XmlError(e)
    }
}

impl From<Error> for PyErr {
    fn from(e: Error) -> Self {
        match e {
            Error::GlobError(e) => PyRuntimeError::new_err(e.to_string()),
            Error::IoError(e) => e.into(),
            Error::XmlError(e) => PyRuntimeError::new_err(e.to_string()),
        }
    }
}
