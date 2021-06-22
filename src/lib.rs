use std::{
    collections::HashMap,
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
    xmatch::{GaiaId, GaiaOrdinals, GaiaStar, HipId, HipOrdinals, HipStar},
};

#[pyfunction]
fn check_hip_ids(gaia_dir: String) -> PyResult<()> {
    let pattern = Glob::new("**/gaiaedr3-hip2-*.vot.gz")
        .map_err(Error::new)?
        .compile_matcher();
    let mut hip_to_gaia: HashMap<HipId, Vec<(GaiaId, f64)>> = HashMap::new();
    let mut gaia_to_hip: HashMap<GaiaId, Vec<HipId>> = HashMap::new();
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

            hip_to_gaia
                .entry(hip_star.hip)
                .and_modify(|v| v.push((gaia_star.source_id, score)))
                .or_insert(vec![(gaia_star.source_id, score)]);
            gaia_to_hip
                .entry(gaia_star.source_id)
                .and_modify(|v| v.push(hip_star.hip))
                .or_insert(vec![hip_star.hip]);
        }
    }

    let hip_stars = hip_to_gaia.len();
    let mut best_unshared: usize = 0;

    hip_to_gaia.into_iter().for_each(|(_hip, mut v)| {
        v.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let source_id = v.first().unwrap().0;
        if gaia_to_hip[&source_id].len() == 1 {
            best_unshared += 1;
        }
    });

    println!("HIP stars matched: {}", hip_stars);
    println!("Best match unshared: {}", best_unshared);

    Ok(())
}

#[pymodule]
fn celestia_gaia(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(check_hip_ids, m)?)?;

    Ok(())
}
