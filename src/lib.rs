use std::{
    collections::HashMap,
    fs::{read_dir, File},
};

use globset::Glob;
use pyo3::{prelude::*, wrap_pyfunction};

mod astro;
mod error;
mod votable;

use error::Error;
use votable::VotableReader;

#[pyfunction]
fn check_hip_ids(gaia_dir: String) -> PyResult<()> {
    let pattern = Glob::new("**/gaiaedr3-hip2-*.vot.gz")
        .map_err(Error::new)?
        .compile_matcher();
    let mut hip_to_gaia: HashMap<i32, Vec<i64>> = HashMap::new();
    let mut gaia_to_hip: HashMap<i64, Vec<i32>> = HashMap::new();
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
            .ok_or(Error::parse("Missing hip field"))?;
        let gaia_ordinal = reader
            .ordinal(b"source_id")
            .ok_or(Error::parse("Missing source_id field"))?;
        while let Some(accessor) = reader.read()? {
            let hip = accessor.read_i32(hip_ordinal)?.unwrap();
            let source_id = accessor.read_i64(gaia_ordinal)?.unwrap();
            hip_to_gaia
                .entry(hip)
                .and_modify(|v| v.push(source_id))
                .or_insert(vec![source_id]);
            gaia_to_hip
                .entry(source_id)
                .and_modify(|v| v.push(hip))
                .or_insert(vec![hip]);
        }
    }

    let hip_stars = hip_to_gaia.len();
    let mut unique_gaia: usize = 0;
    let mut unique_gaia_unshared: usize = 0;

    hip_to_gaia.into_iter().for_each(|(_, source_ids)| {
        if source_ids.len() == 1 {
            unique_gaia += 1;
            if gaia_to_hip[&source_ids[0]].len() == 1 {
                unique_gaia_unshared += 1;
            }
        }
    });

    println!("HIP stars matched: {}", hip_stars);
    println!("Unique matches: {}", unique_gaia);
    println!("Unique unshared: {}", unique_gaia_unshared);

    Ok(())
}

#[pymodule]
fn celestia_gaia(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(check_hip_ids, m)?)?;

    Ok(())
}
