use std::{
    fs::{read_dir, File},
    path::Path,
};

use globset::Glob;
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

fn crossmatch_directory<A, B>(path: &Path, pattern: &str) -> Result<(), Error>
where
    A: VotableRecord + Crossmatchable<B>,
    B: VotableRecord,
{
    let mut crossmatcher = Crossmatcher::<A, B>::new();

    let pattern = Glob::new(pattern).map_err(Error::new)?.compile_matcher();
    for entry_result in read_dir(path)? {
        let entry = entry_result?;
        if !entry.metadata()?.is_file() {
            continue;
        }

        let entry_path = entry.path();
        if !pattern.is_match(&entry_path) {
            continue;
        }

        let file = File::open(entry_path)?;
        let reader = VotableReader::new(file)?;
        crossmatcher.add_reader(reader)?;
    }

    crossmatcher.finalize();

    Ok(())
}

#[pymodule]
fn celestia_gaia(_py: Python, m: &PyModule) -> PyResult<()> {
    #[pyfn(m, "build_hip_xmatch")]
    #[text_signature = "(gaia_dir, /)"]
    fn build_hip_xmatch_py<'py>(_py: Python<'py>, gaia_dir: &PyAny) -> PyResult<()> {
        crossmatch_directory::<HipStar, GaiaStar>(gaia_dir.str()?.to_str()?.as_ref(), HIP_PATTERN)?;
        Ok(())
    }

    #[pyfn(m, "build_tyc_xmatch")]
    #[text_signature = "(gaia_dir, /)"]
    fn build_tyc_xmatch_py<'py>(_py: Python<'py>, gaia_dir: &PyAny) -> PyResult<()> {
        crossmatch_directory::<TycStar, GaiaStar>(gaia_dir.str()?.to_str()?.as_ref(), TYC_PATTERN)?;
        Ok(())
    }

    Ok(())
}
