use std::{
    fs::{read_dir, File},
    io::{self, BufRead, BufReader, ErrorKind},
    path::Path,
};

use bitvec::prelude::*;
use flate2::read::GzDecoder;
use globset::Glob;
use pyo3::{prelude::*, wrap_pyfunction};

mod error;
mod votable;

use votable::VotableReader;

const HIP_SIZE: usize = 120404;

fn set_hip_ids(path: impl AsRef<Path>) -> io::Result<(BitVec, Vec<f32>)> {
    let mut hip_ids = bitvec![1; 120404];
    let mut magnitudes = vec![f32::NAN; HIP_SIZE];
    let file = File::open(path)?;
    let decoder = GzDecoder::new(file);
    let mut reader = BufReader::new(decoder);
    let mut line_buf = String::with_capacity(277);
    while reader.read_line(&mut line_buf)? != 0 {
        let hip: usize = line_buf[0..6]
            .trim()
            .parse()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse HIP id"))?;
        hip_ids.set(hip - 1, false);

        let hp_mag = line_buf[129..136]
            .trim()
            .parse()
            .map_err(|_| io::Error::new(ErrorKind::InvalidData, "Could not parse Hpmag"))?;

        magnitudes[hip - 1] = hp_mag;

        line_buf.clear();
    }

    Ok((hip_ids, magnitudes))
}

fn check_hip_ids_impl(
    hip_file: impl AsRef<Path>,
    gaia_dir: impl AsRef<Path>,
) -> Result<(), error::Error> {
    let (mut hip_ids, magnitudes) = set_hip_ids(hip_file)?;
    let pattern = Glob::new("**/gaiaedr3-hip2-*.vot.gz")?.compile_matcher();
    for entry_result in read_dir(gaia_dir)? {
        let entry = entry_result?;
        if !entry.metadata()?.is_file() {
            continue;
        }

        let entry_path = entry.path();
        if !pattern.is_match(&entry_path) {
            continue;
        }

        let file = File::open(entry_path)?;
        let mut reader = VotableReader::new(file)?;
        let hip_ordinal = reader
            .ordinal(b"hip")
            .ok_or(io::Error::new(ErrorKind::InvalidData, "Missing HIP field"))?;
        while let Some(accessor) = reader.read()? {
            let hip = accessor.read_i32(hip_ordinal)?.unwrap();
            hip_ids.set(hip as usize - 1, true);
        }
    }

    hip_ids
        .iter()
        .enumerate()
        .filter(|(_, b)| *b == false)
        .for_each(|(hip, _)| println!("HIP {} ({})", hip + 1, magnitudes[hip]));

    Ok(())
}

#[pyfunction]
fn check_hip_ids(hip_file: String, gaia_dir: String) -> PyResult<()> {
    check_hip_ids_impl(hip_file, gaia_dir).map_err(|e| e.into())
}

#[pymodule]
fn celestia_gaia(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(check_hip_ids, m)?)?;

    Ok(())
}
