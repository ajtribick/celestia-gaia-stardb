use std::io::Read;

use crate::{
    astro::{ProperMotion, SkyCoords, Squarable, MAS_TO_DEG},
    error::Error,
    votable::{Ordinals, RecordAccessor, VotableReader, VotableRecord},
};

use super::{Crossmatchable, GaiaStar};

const HIP_EPOCH: f64 = 1991.25;

#[derive(Debug)]
pub struct HipOrdinals {
    hip: usize,
    hip_ra: usize,
    hip_dec: usize,
    hp_mag: usize,
}

impl Ordinals for HipOrdinals {
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, Error> {
        Ok(Self {
            hip: reader.ordinal(b"hip")?,
            hip_ra: reader.ordinal(b"hip_ra")?,
            hip_dec: reader.ordinal(b"hip_dec")?,
            hp_mag: reader.ordinal(b"hp_mag")?,
        })
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub struct HipId(pub i32);

#[derive(Debug)]
pub struct HipStar {
    pub hip: HipId,
    pub coords: SkyCoords,
    pub hp_mag: f64,
}

impl VotableRecord for HipStar {
    type Ordinals = HipOrdinals;
    type Id = HipId;

    fn from_accessor(accessor: &RecordAccessor, ordinals: &HipOrdinals) -> Result<Self, Error> {
        Ok(Self {
            hip: HipId(
                accessor
                    .read_i32(ordinals.hip)?
                    .ok_or(Error::missing_id("hip"))?,
            ),
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.hip_ra)?,
                dec: accessor.read_f64(ordinals.hip_dec)?,
            },
            hp_mag: accessor.read_f64(ordinals.hp_mag)?,
        })
    }

    fn id(&self) -> HipId {
        self.hip
    }
}

impl Crossmatchable<GaiaStar> for HipStar {
    fn score(&self, gaia_star: &GaiaStar) -> f64 {
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
