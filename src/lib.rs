use std::{
    collections::HashSet,
    fs::{read_dir, File},
};

use globset::Glob;
use pyo3::{prelude::*, wrap_pyfunction};

mod astro;
mod error;
mod votable;
mod xmatch;

use crate::{
    error::Error,
    votable::VotableReader,
    xmatch::{GaiaOrdinals, GaiaStar, HipOrdinals, HipStar, TycOrdinals, TycStar},
};

#[pyfunction]
fn check_hip_ids(gaia_dir: String) -> PyResult<()> {
    let pattern = Glob::new("**/gaiaedr3-hip2-*.vot.gz")
        .map_err(Error::new)?
        .compile_matcher();
    let mut hip_ids = HashSet::new();
    let mut gaia_ids = HashSet::new();
    let mut matches = Vec::new();
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
        let gaia_ordinals = GaiaOrdinals::from_reader(&reader)?;
        let ordinals = HipOrdinals::from_reader(&reader)?;
        while let Some(accessor) = reader.read()? {
            let hip_star = HipStar::from_accessor(&accessor, &ordinals)?;
            let gaia_star = GaiaStar::from_accessor(&accessor, &gaia_ordinals)?;
            let score = hip_star.score(&gaia_star);

            hip_ids.insert(hip_star.hip);
            gaia_ids.insert(gaia_star.source_id);
            matches.push((hip_star.hip, gaia_star.source_id, score));
        }
    }

    println!("Found {} distinct HIP stars", hip_ids.len());
    matches.sort_by(|(_, _, a), (_, _, b)| b.partial_cmp(a).unwrap());
    while let Some((hip, gaia, _)) = matches.pop() {
        if gaia_ids.contains(&gaia) && hip_ids.remove(&hip) {
            gaia_ids.remove(&gaia);
        }
    }

    println!("{} unmatched", hip_ids.len());

    Ok(())
}

#[pyfunction]
fn check_tyc_ids(gaia_dir: String) -> PyResult<()> {
    let pattern = Glob::new("**/gaiaedr3-tyctdsc-*.vot.gz")
        .map_err(Error::new)?
        .compile_matcher();
    let mut tyc_ids = HashSet::new();
    let mut gaia_ids = HashSet::new();
    let mut matches = Vec::new();
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
        let gaia_ordinals = GaiaOrdinals::from_reader(&reader)?;
        let ordinals = TycOrdinals::from_reader(&reader)?;
        while let Some(accessor) = reader.read()? {
            let tyc_star = TycStar::from_accessor(&accessor, &ordinals)?;
            let gaia_star = GaiaStar::from_accessor(&accessor, &gaia_ordinals)?;
            let score = tyc_star.score(&gaia_star);

            tyc_ids.insert(tyc_star.tyc);
            gaia_ids.insert(gaia_star.source_id);
            matches.push((tyc_star.tyc, gaia_star.source_id, score));
        }
    }

    println!("Found {} distinct TYC stars", tyc_ids.len());
    matches.sort_by(|(_, _, a), (_, _, b)| b.partial_cmp(a).unwrap());
    while let Some((tyc, gaia, _)) = matches.pop() {
        if gaia_ids.contains(&gaia) && tyc_ids.remove(&tyc) {
            gaia_ids.remove(&gaia);
        }
    }

    println!("{} unmatched", tyc_ids.len());

    Ok(())
}

#[pymodule]
fn celestia_gaia(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(check_hip_ids, m)?)?;
    m.add_function(wrap_pyfunction!(check_tyc_ids, m)?)?;

    Ok(())
}
