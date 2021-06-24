use std::io::Read;

use crate::{
    astro::{ProperMotion, SkyCoords, Squarable, MAS_TO_DEG},
    error::Error,
    votable::{Ordinals, RecordAccessor, VotableReader, VotableRecord},
};

use super::{Crossmatchable, GaiaStar};

#[derive(Debug)]
pub struct TycOrdinals {
    id_tycho: usize,
    ra_deg: usize,
    de_deg: usize,
    bt_mag: usize,
    vt_mag: usize,
    ep_ra1990: usize,
    ep_de1990: usize,
}

impl Ordinals for TycOrdinals {
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, Error> {
        Ok(Self {
            id_tycho: reader.ordinal(b"id_tycho")?,
            ra_deg: reader.ordinal(b"tyc_ra")?,
            de_deg: reader.ordinal(b"tyc_dec")?,
            bt_mag: reader.ordinal(b"bt_mag")?,
            vt_mag: reader.ordinal(b"vt_mag")?,
            ep_ra1990: reader.ordinal(b"ep_ra1990")?,
            ep_de1990: reader.ordinal(b"ep_de1990")?,
        })
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub struct TycId(pub i64);

#[derive(Debug)]
pub struct TycStar {
    pub tyc: TycId,
    pub coords: SkyCoords,
    pub bt_mag: f32,
    pub vt_mag: f32,
    pub epoch_ra: f64,
    pub epoch_dec: f64,
}

fn g_vt(bp_rp: f32) -> f32 {
    -0.01077 + bp_rp * (-0.0682 + bp_rp * (-0.2387 + bp_rp * 0.02342))
}

fn g_bt(bp_rp: f32) -> f32 {
    -0.004288
        + bp_rp
            * (-0.8547 + bp_rp * (0.1244 + bp_rp * (-0.9085 + bp_rp * (0.4843 + bp_rp * -0.06814))))
}

impl VotableRecord for TycStar {
    type Ordinals = TycOrdinals;
    type Id = TycId;

    fn from_accessor(accessor: &RecordAccessor, ordinals: &TycOrdinals) -> Result<Self, Error> {
        Ok(Self {
            tyc: TycId(
                accessor
                    .read_i64(ordinals.id_tycho)?
                    .ok_or(Error::missing_id("id_tycho"))?,
            ),
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.ra_deg)?,
                dec: accessor.read_f64(ordinals.de_deg)?,
            },
            bt_mag: accessor.read_f32(ordinals.bt_mag)?,
            vt_mag: accessor.read_f32(ordinals.vt_mag)?,
            epoch_ra: accessor.read_f32(ordinals.ep_ra1990)? as f64 + 1990.0,
            epoch_dec: accessor.read_f32(ordinals.ep_de1990)? as f64 + 1990.0,
        })
    }

    fn id(&self) -> TycId {
        self.tyc
    }
}

impl Crossmatchable<GaiaStar> for TycStar {
    fn score(&self, gaia_star: &GaiaStar) -> f64 {
        let mean_epoch = self.epoch_ra * 0.5 + self.epoch_dec * 0.5;
        let epoch_diff = gaia_star.epoch - mean_epoch;
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
            .apply_pm(&pm, rv, parallax, gaia_star.epoch, mean_epoch);
        let dist = self.coords.ang_dist(&propagated);
        assert!(
            !dist.is_nan(),
            "TYC {}, Gaia {}\n{:?}\n{:?}",
            self.tyc.0,
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

        let bp_rp = bp_mag - rp_mag;

        let idx_diff = if self.bt_mag.is_nan() {
            gaia_star.phot_g_mean_mag - self.vt_mag - g_vt(bp_rp)
        } else if self.vt_mag.is_nan() {
            gaia_star.phot_g_mean_mag - self.bt_mag - g_bt(bp_rp)
        } else {
            gaia_star.phot_g_mean_mag
                - (self.vt_mag + self.bt_mag + g_vt(bp_rp) + g_bt(bp_rp)) / 2.0
        };

        assert!(!idx_diff.is_nan());

        pm_diff + (idx_diff as f64 / 0.06).sqr() + (dist / MAS_TO_DEG).sqr()
    }
}
