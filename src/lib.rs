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

use std::{
    borrow::Cow,
    collections::{hash_map::Entry, HashMap, HashSet},
    fs::{read_dir, File},
    iter::FromIterator,
    path::{Path, PathBuf},
};

use flate2::{write::GzEncoder, Compression};
use globset::Glob;
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

mod astro;
mod csv;
mod error;
mod hip2dist;
mod votable;
mod xmatch;

use crate::{
    astro::{HipId, TycId},
    error::AppError,
    hip2dist::estimate_distances,
    votable::VotableReader,
    xmatch::Crossmatcher,
};

const HIP_PATTERN: &str = "**/gaiaedr3-hip2-*.vot.gz";
const TYC_PATTERN: &str = "**/gaiaedr3-tyctdsc-*.vot.gz";
const XMATCH_PATTERN: &str = "**/xmatch-*.vot.gz";
const DISTANCE_PATTERN: &str = "**/gaiaedr3-distance-*.vot.gz";

fn load_tyc2hip(path: &Path) -> Result<HashMap<TycId, HipId>, AppError> {
    let mut path = PathBuf::from(path);
    path.push("tyc2tdsc_hip_xmatch.vot.gz");

    let file = File::open(path)?;
    let mut reader = VotableReader::new(file)?;

    let id_tycho_col = reader.ordinal(b"id_tycho")?;
    let hip_col = reader.ordinal(b"hip")?;
    let comp_col = reader.ordinal(b"cmp")?;

    let mut hip2tyc = HashMap::new();
    while let Some(accessor) = reader.read()? {
        let id_tycho = TycId(
            accessor
                .read_i64(id_tycho_col)?
                .ok_or(AppError::missing_id("id_tycho"))?,
        );
        let hip = HipId(
            accessor
                .read_i32(hip_col)?
                .ok_or(AppError::missing_id("hip"))?,
        );
        let cmp = accessor.read_char::<2>(comp_col)?;

        match hip2tyc.entry(hip) {
            Entry::Vacant(v) => {
                v.insert((id_tycho, cmp));
            }
            Entry::Occupied(mut o) => {
                if &cmp < &o.get().1 {
                    o.insert((id_tycho, cmp));
                }
            }
        }
    }

    Ok(hip2tyc.iter().map(|(h, (t, _))| (*t, *h)).collect())
}

fn full_crossmatch(path: &Path, output_name: &str) -> Result<(), AppError> {
    let tyc2hip = load_tyc2hip(path)?;
    let mut crossmatcher = Crossmatcher::new(tyc2hip);

    let hip_pattern = Glob::new(HIP_PATTERN)?.compile_matcher();
    let tyc_pattern = Glob::new(TYC_PATTERN)?.compile_matcher();
    for entry in read_dir(path)? {
        let entry = entry?;
        if !entry.metadata()?.is_file() {
            continue;
        }
        let entry_path = entry.path();
        if hip_pattern.is_match(&entry_path) {
            let file = File::open(entry_path)?;
            let reader = VotableReader::new(file)?;
            crossmatcher.add_hip(reader)?;
        } else if tyc_pattern.is_match(&entry_path) {
            let file = File::open(entry_path)?;
            let reader = VotableReader::new(file)?;
            crossmatcher.add_tyc(reader)?;
        }
    }

    let mut output_path = path.to_path_buf();
    output_path.push(output_name);

    let file = File::create(output_path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    crossmatcher.finalize(&mut encoder)?;
    encoder.finish()?.sync_all()?;

    Ok(())
}

fn get_required_dist_source_ids(path: &Path) -> Result<Vec<i64>, AppError> {
    let pattern = Glob::new(XMATCH_PATTERN)?.compile_matcher();
    let mut source_ids = HashSet::new();
    for entry_result in read_dir(path)? {
        let entry = entry_result?;
        let entry_path = entry.path();
        if !pattern.is_match(&entry_path) || !entry.metadata()?.is_file() {
            continue;
        }

        let file = File::open(entry_path)?;
        let mut reader = VotableReader::new(file)?;

        let ordinal = reader.ordinal(b"source_id")?;
        while let Some(accessor) = reader.read()? {
            let source_id = accessor
                .read_i64(ordinal)?
                .ok_or(AppError::MissingField(Cow::Borrowed(b"source_id")))?;
            source_ids.insert(source_id);
        }
    }

    let pattern = Glob::new(DISTANCE_PATTERN)?.compile_matcher();
    for entry_result in read_dir(path)? {
        let entry = entry_result?;
        let entry_path = entry.path();
        if !pattern.is_match(&entry_path) || !entry.metadata()?.is_file() {
            continue;
        }

        let file = File::open(entry_path)?;
        let mut reader = VotableReader::new(file)?;

        let ordinal = reader.ordinal(b"source_id")?;
        while let Some(accessor) = reader.read()? {
            if let Some(source_id) = accessor.read_i64(ordinal)? {
                source_ids.remove(&source_id);
            }
        }
    }

    let mut vec = Vec::from_iter(source_ids.into_iter());
    vec.sort_unstable();
    Ok(vec)
}

fn apply_distances(gaia_dir: &Path, source_ids: &[i64]) -> Result<Vec<f32>, AppError> {
    let mut geometric = vec![f32::NAN; source_ids.len()];
    let source_indices: HashMap<_, _> = source_ids
        .into_iter()
        .enumerate()
        .map(|(i, g)| (*g, i))
        .collect();

    let pattern = Glob::new(DISTANCE_PATTERN)?.compile_matcher();
    for entry_result in read_dir(gaia_dir)? {
        let entry = entry_result?;
        let entry_path = entry.path();
        if !pattern.is_match(&entry_path) || !entry.metadata()?.is_file() {
            continue;
        }

        let file = File::open(entry_path)?;
        let mut reader = VotableReader::new(file)?;
        let source_ordinal = reader.ordinal(b"source_id")?;
        let geometric_ordinal = reader.ordinal(b"r_med_geo")?;
        let photogeometric_ordinal = reader.ordinal(b"r_med_photogeo")?;
        while let Some(accessor) = reader.read()? {
            let source_id = accessor
                .read_i64(source_ordinal)?
                .ok_or(AppError::MissingField(Cow::Borrowed(b"source_id")))?;
            if let Some(&index) = source_indices.get(&source_id) {
                let mut distance = accessor.read_f32(photogeometric_ordinal)?;
                if distance.is_nan() {
                    distance = accessor.read_f32(geometric_ordinal)?;
                }

                geometric[index] = distance;
            }
        }
    }

    Ok(geometric)
}

#[pymodule]
fn celestia_gaia(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m)]
    #[pyo3(name = "build_xmatch", text_signature = "()")]
    fn build_xmatch_py<'py>(_py: Python<'py>, gaia_dir: &PyAny, output_name: &str) -> PyResult<()> {
        full_crossmatch(gaia_dir.str()?.to_str()?.as_ref(), output_name)?;
        Ok(())
    }

    #[pyfn(m)]
    #[pyo3(
        name = "get_required_dist_source_ids",
        text_signature = "(gaia_dir, /)"
    )]
    fn get_required_dist_source_ids_py<'py>(
        py: Python<'py>,
        gaia_dir: &PyAny,
    ) -> PyResult<&'py PyArray1<i64>> {
        Ok(get_required_dist_source_ids(gaia_dir.str()?.to_str()?.as_ref())?.into_pyarray(py))
    }

    #[pyfn(m)]
    #[pyo3(name = "apply_distances", text_signature = "(gaia_dir, source_ids, /)")]
    fn apply_distances_py<'py>(
        py: Python<'py>,
        gaia_dir: &'py PyAny,
        source_ids: &'py PyArray1<i64>,
    ) -> PyResult<&'py PyArray1<f32>> {
        let distances = apply_distances(
            gaia_dir.str()?.to_str()?.as_ref(),
            source_ids.readonly().as_slice()?,
        )?;
        Ok(distances.into_pyarray(py))
    }

    #[pyfn(m)]
    #[pyo3(
        name = "estimate_distances",
        text_signature = "(prior_file, hip2_file, output_file, /)"
    )]
    fn estimate_distances_py<'py>(
        _py: Python<'py>,
        prior_file: &'py PyAny,
        hip2_file: &'py PyAny,
        output_file: &'py PyAny,
    ) -> PyResult<()> {
        estimate_distances(
            prior_file.str()?.to_str()?.to_owned(),
            hip2_file.str()?.to_str()?.to_owned(),
            output_file.str()?.to_str()?.to_owned(),
        )
        .map_err(Into::into)
    }

    Ok(())
}
