use std::{fmt, hash::Hash, io::Read, num::NonZeroUsize};

use super::error::Error;

mod read;
mod write;

pub use read::{RecordAccessor, VotableReader};
pub use write::VotableWriter;

const VOTABLE_NS: &[u8] = b"http://www.ivoa.net/xml/VOTable/v1.3";

pub trait Ordinals: Sized {
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, Error>;
}

pub trait VotableRecord: Sized + 'static {
    type Ordinals: Ordinals;
    type Id: Eq + Hash + Copy;

    fn from_accessor(accessor: &RecordAccessor, ordinals: &Self::Ordinals) -> Result<Self, Error>;
    fn id(&self) -> Self::Id;
    fn fields() -> &'static [FieldInfo<Self>];
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    Short,
    Int,
    Long,
    Float,
    Double,
    Char,
}

impl DataType {
    fn parse_bytes(bstr: &[u8]) -> Result<Self, Error> {
        match bstr {
            b"short" => Ok(Self::Short),
            b"int" => Ok(Self::Int),
            b"long" => Ok(Self::Long),
            b"float" => Ok(Self::Float),
            b"double" => Ok(Self::Double),
            b"char" => Ok(Self::Char),
            _ => Err(Error::parse("Unsupported data type")),
        }
    }

    fn width(&self) -> Option<NonZeroUsize> {
        match self {
            Self::Short => NonZeroUsize::new(2),
            Self::Int => NonZeroUsize::new(4),
            Self::Long => NonZeroUsize::new(8),
            Self::Float => NonZeroUsize::new(4),
            Self::Double => NonZeroUsize::new(8),
            Self::Char => None,
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
        }
    }
}

enum FieldGetter<R> {
    Short(fn(&R) -> Option<i16>),
    Int(fn(&R) -> Option<i32>),
    Long(fn(&R) -> Option<i64>),
    Float(fn(&R) -> f32),
    Double(fn(&R) -> f64),
    Char(fn(&R) -> Option<&[u8]>),
}

pub struct FieldInfo<R: VotableRecord> {
    name: &'static str,
    unit: Option<&'static str>,
    ucd: Option<&'static str>,
    description: &'static str,
    getter: FieldGetter<R>,
}

impl<R: VotableRecord> FieldInfo<R> {
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

    pub fn int(
        name: &'static str,
        unit: Option<&'static str>,
        ucd: Option<&'static str>,
        description: &'static str,
        getter: fn(&R) -> Option<i32>,
    ) -> Self {
        Self {
            name,
            unit,
            ucd,
            description,
            getter: FieldGetter::Int(getter),
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

    pub fn char(
        name: &'static str,
        unit: Option<&'static str>,
        ucd: Option<&'static str>,
        description: &'static str,
        getter: fn(&R) -> Option<&[u8]>,
    ) -> Self {
        Self {
            name,
            unit,
            ucd,
            description,
            getter: FieldGetter::Char(getter),
        }
    }
}
