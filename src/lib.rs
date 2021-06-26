use std::{
    borrow::Cow,
    collections::{HashMap, HashSet},
    fs::{read_dir, File},
    iter::FromIterator,
    path::Path,
};

use flate2::{write::GzEncoder, Compression};
use globset::Glob;
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

mod astro;
mod error;
mod votable;
mod xmatch;

use crate::{
    error::Error,
    votable::{VotableReader, VotableRecord},
    xmatch::{Crossmatchable, Crossmatcher, GaiaStar, HipStar, TycStar},
};

const HIP_PATTERN: &str = "**/gaiaedr3-hip2-*.vot.gz";
const TYC_PATTERN: &str = "**/gaiaedr3-tyctdsc-*.vot.gz";
const XMATCH_PATTERN: &str = "**/xmatch-*.vot.gz";
const DISTANCE_PATTERN: &str = "**/gaiaedr3-distance-*.vot.gz";

fn crossmatch_directory<A, B>(path: &Path, pattern: &str, output_name: &str) -> Result<(), Error>
where
    A: VotableRecord + Crossmatchable<B>,
    B: VotableRecord,
{
    let mut crossmatcher = Crossmatcher::<A, B>::new();

    let pattern = Glob::new(pattern)?.compile_matcher();
    for entry_result in read_dir(path)? {
        let entry = entry_result?;
        let entry_path = entry.path();
        if !pattern.is_match(&entry_path) || !entry.metadata()?.is_file() {
            continue;
        }

        let file = File::open(entry_path)?;
        let reader = VotableReader::new(file)?;
        crossmatcher.add_reader(reader)?;
    }

    let mut output_path = path.to_path_buf();
    output_path.push(output_name);

    let file = File::create(output_path)?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    crossmatcher.finalize(&mut encoder)?;
    encoder.finish()?.sync_all()?;

    Ok(())
}

fn get_xmatch_source_ids(path: &Path) -> Result<Vec<i64>, Error> {
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
                .ok_or(Error::MissingField(Cow::Borrowed(b"source_id")))?;
            source_ids.insert(source_id);
        }
    }

    let mut vec = Vec::from_iter(source_ids.into_iter());
    vec.sort_unstable();
    Ok(vec)
}

fn apply_distances(gaia_dir: &Path, source_ids: &[i64]) -> Result<Vec<f32>, Error> {
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
                .ok_or(Error::MissingField(Cow::Borrowed(b"source_id")))?;
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
    #[pyfn(m, "build_hip_xmatch")]
    #[text_signature = "(gaia_dir, output_name, /)"]
    fn build_hip_xmatch_py<'py>(
        _py: Python<'py>,
        gaia_dir: &PyAny,
        output_name: &str,
    ) -> PyResult<()> {
        crossmatch_directory::<HipStar, GaiaStar>(
            gaia_dir.str()?.to_str()?.as_ref(),
            HIP_PATTERN,
            output_name,
        )?;
        Ok(())
    }

    #[pyfn(m, "build_tyc_xmatch")]
    #[text_signature = "(gaia_dir, output_name, /)"]
    fn build_tyc_xmatch_py<'py>(
        _py: Python<'py>,
        gaia_dir: &PyAny,
        output_name: &str,
    ) -> PyResult<()> {
        crossmatch_directory::<TycStar, GaiaStar>(
            gaia_dir.str()?.to_str()?.as_ref(),
            TYC_PATTERN,
            output_name,
        )?;
        Ok(())
    }

    #[pyfn(m, "get_source_ids")]
    #[text_signature = "(gaia_dir, /)"]
    fn get_source_ids_py<'py>(py: Python<'py>, gaia_dir: &PyAny) -> PyResult<&'py PyArray1<i64>> {
        Ok(get_xmatch_source_ids(gaia_dir.str()?.to_str()?.as_ref())?.into_pyarray(py))
    }

    #[pyfn(m, "apply_distances")]
    #[text_signature = "(gaia_dir, source_ids, /)"]
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

    Ok(())
}
