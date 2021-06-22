use std::io::Read;

use super::{
    astro::{ProperMotion, SkyCoords, Squarable, MAS_TO_DEG},
    error::Error,
    votable::{RecordAccessor, VotableReader},
};

const HIP_EPOCH: f64 = 1991.25;

#[derive(Debug)]
pub struct GaiaOrdinals {
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

impl GaiaOrdinals {
    pub fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, Error> {
        Ok(Self {
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

#[derive(Debug)]
pub struct HipOrdinals {
    hip: usize,
    hip_ra: usize,
    hip_dec: usize,
    hp_mag: usize,
}

impl HipOrdinals {
    pub fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, Error> {
        Ok(Self {
            hip: reader.ordinal(b"hip")?,
            hip_ra: reader.ordinal(b"hip_ra")?,
            hip_dec: reader.ordinal(b"hip_dec")?,
            hp_mag: reader.ordinal(b"hp_mag")?,
        })
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub struct GaiaId(pub i64);

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub struct HipId(pub i32);

#[derive(Debug)]
pub struct GaiaStar {
    pub source_id: GaiaId,
    pub coords: SkyCoords,
    pub parallax: f64,
    pub parallax_error: f64,
    pub rv: f64,
    pub epoch: f64,
    pub pm: ProperMotion,
    pub e_pm: ProperMotion,
    pub phot_g_mean_mag: f32,
    pub phot_bp_mean_mag: f32,
    pub phot_rp_mean_mag: f32,
    pub astrometric_params_solved: i16,
}

impl GaiaStar {
    pub fn from_accessor(
        accessor: &RecordAccessor,
        ordinals: &GaiaOrdinals,
    ) -> Result<Self, Error> {
        Ok(Self {
            source_id: GaiaId(accessor.read_i64(ordinals.source_id)?.unwrap()),
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.ra)?,
                dec: accessor.read_f64(ordinals.dec)?,
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
            astrometric_params_solved: accessor
                .read_i16(ordinals.astrometric_params_solved)?
                .unwrap_or(0),
        })
    }
}

#[derive(Debug)]
pub struct HipStar {
    pub hip: HipId,
    pub coords: SkyCoords,
    pub hp_mag: f64,
}

impl HipStar {
    pub fn from_accessor(accessor: &RecordAccessor, ordinals: &HipOrdinals) -> Result<Self, Error> {
        Ok(Self {
            hip: HipId(accessor.read_i32(ordinals.hip)?.unwrap()),
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.hip_ra)?,
                dec: accessor.read_f64(ordinals.hip_dec)?,
            },
            hp_mag: accessor.read_f64(ordinals.hp_mag)?,
        })
    }

    pub fn score(&self, gaia_star: &GaiaStar) -> f64 {
        let epoch_diff = gaia_star.epoch - HIP_EPOCH;
        let rv = if gaia_star.rv.is_nan() {
            0.0
        } else {
            gaia_star.rv
        };
        let parallax = f64::min(gaia_star.parallax, 0.01);
        let pm = ProperMotion {
            pm_ra: if gaia_star.pm.pm_ra.is_nan() {
                0.0
            } else {
                gaia_star.pm.pm_ra
            },
            pm_dec: if gaia_star.pm.pm_dec.is_nan() {
                0.0
            } else {
                gaia_star.pm.pm_dec
            },
        };
        let propagated = gaia_star
            .coords
            .apply_pm(&pm, rv, parallax, gaia_star.epoch, HIP_EPOCH);
        let dist = self.coords.ang_dist(&propagated);
        assert!(
            !dist.is_nan(),
            "HIP {}, Gaia {}\n{:?}\n{:?}",
            self.hip.0,
            gaia_star.source_id.0,
            self.coords,
            propagated
        );

        let calc_pm = ProperMotion {
            pm_ra: (gaia_star.coords.ra - self.coords.ra) * gaia_star.coords.dec.cos() / epoch_diff
                * MAS_TO_DEG,
            pm_dec: (gaia_star.coords.dec - self.coords.dec) / epoch_diff * MAS_TO_DEG,
        };

        let e_pm = ProperMotion {
            pm_ra: f64::min(gaia_star.e_pm.pm_ra, 0.001),
            pm_dec: f64::min(gaia_star.e_pm.pm_dec, 0.001),
        };
        let pm_diff = ((calc_pm.pm_ra - pm.pm_ra) / e_pm.pm_ra).sqr()
            + ((calc_pm.pm_dec - pm.pm_dec) / e_pm.pm_dec).sqr();
        assert!(!pm_diff.is_nan());

        let bp_mag = if gaia_star.phot_bp_mean_mag.is_nan() {
            gaia_star.phot_g_mean_mag
        } else {
            gaia_star.phot_bp_mean_mag
        };
        let rp_mag = if gaia_star.phot_rp_mean_mag.is_nan() {
            gaia_star.phot_g_mean_mag
        } else {
            gaia_star.phot_rp_mean_mag
        };
        let mag_diff = self.hp_mag - (0.91 * bp_mag + 0.09 * rp_mag) as f64;
        assert!(!mag_diff.is_nan());

        pm_diff + (mag_diff / 0.1).sqr() + (dist / MAS_TO_DEG).sqr()
    }
}