use std::{
    collections::HashMap,
    fs::{read_dir, File},
    io::Read,
};

use globset::Glob;
use pyo3::{prelude::*, wrap_pyfunction};

mod astro;
mod error;
mod votable;

use astro::{apply_pm, SkyCoords, ProperMotion, MAS_TO_DEG, Squarable};
use error::Error;
use votable::{RecordAccessor, VotableReader};

const HIP_EPOCH: f64 = 1991.25;

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
struct HipId(pub i32);

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
struct GaiaId(pub i64);

#[derive(Debug)]
struct HipStar {
    hip: HipId,
    coords: SkyCoords,
    hp_mag: f64,
}

impl HipStar {
    pub fn from_accessor(accessor: &RecordAccessor, ordinals: &Ordinals) -> Result<Self, Error> {
        Ok(Self {
            hip: HipId(accessor.read_i32(ordinals.hip)?.unwrap()),
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.hip_ra)?,
                dec: accessor.read_f64(ordinals.hip_dec)?,
            },
            hp_mag: accessor.read_f64(ordinals.hp_mag)?,
        })
    }
}

#[derive(Debug)]
struct GaiaStar {
    source_id: GaiaId,
    coords: SkyCoords,
    parallax: f64,
    parallax_error: f64,
    rv: f64,
    epoch: f64,
    pm: ProperMotion,
    e_pm: ProperMotion,
    phot_g_mean_mag: f32,
    phot_bp_mean_mag: f32,
    phot_rp_mean_mag: f32,
    astrometric_params_solved: i16,
}

impl GaiaStar {
    pub fn from_accessor(accessor: &RecordAccessor, ordinals: &Ordinals) -> Result<Self, Error> {
        Ok(Self {
            source_id: GaiaId(accessor.read_i64(ordinals.source_id)?.unwrap()),
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.ra)?,
                dec:accessor.read_f64(ordinals.dec)?,
            },
            parallax: accessor.read_f64(ordinals.parallax)?,
            parallax_error: accessor.read_f32(ordinals.parallax_error)? as f64,
            rv: accessor.read_f32(ordinals.dr2_radial_velocity)? as f64,
            epoch: accessor.read_f64(ordinals.ref_epoch)?,
            pm: ProperMotion {
                pm_ra: accessor.read_f64(ordinals.pm_ra)?,
                pm_dec: accessor.read_f64(ordinals.pm_dec)?,
            },
            e_pm: ProperMotion {
                pm_ra: accessor.read_f32(ordinals.pm_ra_error)? as f64,
                pm_dec: accessor.read_f32(ordinals.pm_dec_error)? as f64,
            },
            phot_g_mean_mag: accessor.read_f32(ordinals.phot_g_mean_mag)?,
            phot_bp_mean_mag: accessor.read_f32(ordinals.phot_bp_mean_mag)?,
            phot_rp_mean_mag: accessor.read_f32(ordinals.phot_rp_mean_mag)?,
            astrometric_params_solved: accessor.read_i16(ordinals.astrometric_params_solved)?.unwrap_or(0),
        })
    }
}

fn crossmatch_score(hip_star: &HipStar, gaia_star: &GaiaStar) -> f64 {
    let epoch_diff = gaia_star.epoch - HIP_EPOCH;
    let rv = if gaia_star.rv.is_nan() { 0.0 } else { gaia_star.rv };
    let parallax = f64::min(gaia_star.parallax, 0.01);
    let pm = ProperMotion {
        pm_ra: if gaia_star.pm.pm_ra.is_nan() { 0.0 } else { gaia_star.pm.pm_ra },
        pm_dec: if gaia_star.pm.pm_dec.is_nan() { 0.0 } else { gaia_star.pm.pm_dec },
    };
    let propagated = apply_pm(&gaia_star.coords, &pm, rv, parallax, gaia_star.epoch, HIP_EPOCH);
    let dist = hip_star.coords.ang_dist(&propagated);
    assert!(!dist.is_nan(), "HIP {}, Gaia {}\n{:?}\n{:?}", hip_star.hip.0, gaia_star.source_id.0, hip_star.coords, propagated);

    let calc_pm = ProperMotion {
        pm_ra: (gaia_star.coords.ra - hip_star.coords.ra) * gaia_star.coords.dec.cos() / epoch_diff * MAS_TO_DEG,
        pm_dec: (gaia_star.coords.dec - hip_star.coords.dec) / epoch_diff * MAS_TO_DEG,
    };

    let e_pm = ProperMotion {
        pm_ra: f64::min(gaia_star.e_pm.pm_ra, 0.001),
        pm_dec: f64::min(gaia_star.e_pm.pm_dec, 0.001),
    };
    let pm_diff = ((calc_pm.pm_ra - pm.pm_ra) / e_pm.pm_ra).sqr() + ((calc_pm.pm_dec - pm.pm_dec) / e_pm.pm_dec).sqr();
    assert!(!pm_diff.is_nan());

    let bp_mag = if gaia_star.phot_bp_mean_mag.is_nan() { gaia_star.phot_g_mean_mag } else { gaia_star.phot_bp_mean_mag };
    let rp_mag = if gaia_star.phot_rp_mean_mag.is_nan() { gaia_star.phot_g_mean_mag } else { gaia_star.phot_rp_mean_mag };
    let mag_diff = hip_star.hp_mag - (0.91 * bp_mag + 0.09 * rp_mag) as f64;
    assert!(!mag_diff.is_nan());

    pm_diff + (mag_diff / 0.1).sqr() + (dist / MAS_TO_DEG).sqr()
}

#[derive(Debug)]
struct Ordinals {
    hip: usize,
    hip_ra: usize,
    hip_dec: usize,
    hp_mag: usize,
    source_id: usize,
    ra: usize,
    dec: usize,
    parallax: usize,
    parallax_error: usize,
    dr2_radial_velocity: usize,
    ref_epoch: usize,
    pm_ra: usize,
    pm_ra_error: usize,
    pm_dec: usize,
    pm_dec_error: usize,
    phot_g_mean_mag: usize,
    phot_bp_mean_mag: usize,
    phot_rp_mean_mag: usize,
    astrometric_params_solved: usize,
}

impl Ordinals {
    pub fn from_reader<R: Read>(reader: &VotableReader<R>) -> Result<Self, Error> {
        Ok(Self {
            hip: reader.ordinal(b"hip")?,
            hip_ra: reader.ordinal(b"hip_ra")?,
            hip_dec: reader.ordinal(b"hip_dec")?,
            hp_mag: reader.ordinal(b"hp_mag")?,
            source_id: reader.ordinal(b"source_id")?,
            ra: reader.ordinal(b"ra")?,
            dec: reader.ordinal(b"dec")?,
            parallax: reader.ordinal(b"parallax")?,
            parallax_error: reader.ordinal(b"parallax_error")?,
            dr2_radial_velocity: reader.ordinal(b"dr2_radial_velocity")?,
            ref_epoch: reader.ordinal(b"ref_epoch")?,
            pm_ra: reader.ordinal(b"pmra")?,
            pm_ra_error: reader.ordinal(b"pmra_error")?,
            pm_dec: reader.ordinal(b"pmdec")?,
            pm_dec_error: reader.ordinal(b"pmdec_error")?,
            phot_g_mean_mag: reader.ordinal(b"phot_g_mean_mag")?,
            phot_bp_mean_mag: reader.ordinal(b"phot_bp_mean_mag")?,
            phot_rp_mean_mag: reader.ordinal(b"phot_rp_mean_mag")?,
            astrometric_params_solved: reader.ordinal(b"astrometric_params_solved")?,
        })
    }
}

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
        let ordinals = Ordinals::from_reader(&reader)?;
        while let Some(accessor) = reader.read()? {
            let hip_star = HipStar::from_accessor(&accessor, &ordinals)?;
            let gaia_star = GaiaStar::from_accessor(&accessor, &ordinals)?;
            let score = crossmatch_score(&hip_star, &gaia_star);

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
