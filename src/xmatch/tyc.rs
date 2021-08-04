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

use std::io::{Read, Write};

use lazy_static::lazy_static;

use super::{Crossmatchable, GaiaStar};
use crate::{
    astro::{HipId, ProperMotion, SkyCoords, Squarable, TycId, MAS_TO_DEG},
    error::AppError,
    votable::{FieldInfo, Ordinals, RecordAccessor, VotableReader, VotableRecord},
};

#[derive(Debug)]
pub struct TycOrdinals {
    id_tycho: usize,
    hip: usize,
    cmp: usize,
    ra_deg: usize,
    de_deg: usize,
    bt_mag: usize,
    vt_mag: usize,
    ep_ra1990: usize,
    ep_de1990: usize,
}

impl Ordinals for TycOrdinals {
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, AppError> {
        Ok(Self {
            id_tycho: reader.ordinal(b"id_tycho")?,
            hip: reader.ordinal(b"hip")?,
            cmp: reader.ordinal(b"cmp")?,
            ra_deg: reader.ordinal(b"tyc_ra")?,
            de_deg: reader.ordinal(b"tyc_dec")?,
            bt_mag: reader.ordinal(b"bt_mag")?,
            vt_mag: reader.ordinal(b"vt_mag")?,
            ep_ra1990: reader.ordinal(b"ep_ra1990")?,
            ep_de1990: reader.ordinal(b"ep_de1990")?,
        })
    }
}

#[derive(Debug)]
pub struct TycStar {
    pub tyc: TycId,
    pub hip: Option<HipId>,
    pub cmp: Option<Vec<u8>>,
    pub coords: SkyCoords,
    pub bt_mag: f32,
    pub vt_mag: f32,
    pub epoch_ra: f32,
    pub epoch_dec: f32,
}

fn g_vt(bp_rp: f32) -> f32 {
    -0.01077 + bp_rp * (-0.0682 + bp_rp * (-0.2387 + bp_rp * 0.02342))
}

fn g_bt(bp_rp: f32) -> f32 {
    -0.004288
        + bp_rp
            * (-0.8547 + bp_rp * (0.1244 + bp_rp * (-0.9085 + bp_rp * (0.4843 + bp_rp * -0.06814))))
}

lazy_static! {
    static ref TYC_FIELDS: [FieldInfo<TycStar>; 9] = [
        FieldInfo::long(
            "id_tycho",
            None,
            Some("meta.id;meta.dataset"),
            "Numeric Tycho-2 identifier",
            |t| Some(t.tyc.0)
        ),
        FieldInfo::int(
            "hip",
            None,
            Some("meta.id.cross"),
            "Hipparcos number",
            |t| t.hip.map(|h| h.0)
        ),
        FieldInfo::char(
            "cmp",
            None,
            Some("meta.id.part"),
            "Component designation",
            |t| t.cmp.as_deref(),
        ),
        FieldInfo::double(
            "tyc_ra",
            Some("Angle[deg]"),
            Some("pos.eq.ra"),
            "Observed Tycho-2 Right Ascension, ICRS",
            |t| t.coords.ra
        ),
        FieldInfo::double(
            "tyc_dec",
            Some("Angle[deg]"),
            Some("pos.eq.dec"),
            "Observed Tycho-2 Declination, ICRS",
            |t| t.coords.ra
        ),
        FieldInfo::float(
            "bt_mag",
            Some("Magnitude[Mag]"),
            Some("em.opt.B;phot.mag"),
            "Tycho-2 BT magnitude",
            |t| t.bt_mag
        ),
        FieldInfo::float(
            "vt_mag",
            Some("Magnitude[Mag]"),
            Some("em.opt.V;phot.mag"),
            "Tycho-2 VT magnitude",
            |t| t.vt_mag
        ),
        FieldInfo::float(
            "ep_ra1990",
            Some("Time[year]"),
            Some("time.epoch"),
            "Epoch--1990 of raDeg",
            |g| g.epoch_ra
        ),
        FieldInfo::float(
            "ep_de1990",
            Some("Time[year]"),
            Some("time.epoch"),
            "Epoch--1990 of deDeg",
            |g| g.epoch_dec
        ),
    ];
}

impl VotableRecord for TycStar {
    type Ordinals = TycOrdinals;
    type Id = TycId;

    fn from_accessor(accessor: &RecordAccessor, ordinals: &TycOrdinals) -> Result<Self, AppError> {
        Ok(Self {
            tyc: TycId(
                accessor
                    .read_i64(ordinals.id_tycho)?
                    .ok_or(AppError::missing_id("id_tycho"))?,
            ),
            hip: accessor.read_i32(ordinals.hip)?.map(HipId),
            cmp: accessor.read_char(ordinals.cmp)?,
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.ra_deg)?,
                dec: accessor.read_f64(ordinals.de_deg)?,
            },
            bt_mag: accessor.read_f32(ordinals.bt_mag)?,
            vt_mag: accessor.read_f32(ordinals.vt_mag)?,
            epoch_ra: accessor.read_f32(ordinals.ep_ra1990)?,
            epoch_dec: accessor.read_f32(ordinals.ep_de1990)?,
        })
    }

    fn id(&self) -> TycId {
        self.tyc
    }

    fn fields() -> &'static [FieldInfo<Self>] {
        TYC_FIELDS.as_ref()
    }
}

impl Crossmatchable<GaiaStar> for TycStar {
    fn score(&self, gaia_star: &GaiaStar) -> f64 {
        let mean_epoch = (self.epoch_ra * 0.5 + self.epoch_dec * 0.5) as f64 + 1990.0;
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
        let propagated =
            gaia_star
                .coords
                .apply_pm(&pm, rv as f64, parallax, gaia_star.epoch, mean_epoch);
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
            pm_ra: f64::min(gaia_star.e_pm.pm_ra as f64, 0.001),
            pm_dec: f64::min(gaia_star.e_pm.pm_dec as f64, 0.001),
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

pub fn tyc_csv_crossmatch(
    mut reader: VotableReader<impl Read>,
    writer: &mut impl Write,
) -> Result<(), AppError> {
    let mut results = Vec::with_capacity(2561887);
    let tyc_ordinal = reader.ordinal(b"id_tycho")?;
    let source_id_ordinal = reader.ordinal(b"source_id")?;
    writeln!(writer, "tyc,source_id")?;
    while let Some(record) = reader.read()? {
        let tyc = record
            .read_i64(tyc_ordinal)?
            .ok_or(AppError::missing_id("id_tycho"))?;
        let source_id = record
            .read_i64(source_id_ordinal)?
            .ok_or(AppError::missing_id("source_id"))?;
        results.push((tyc, source_id));
    }

    results.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));

    for (tyc, source_id) in results {
        let tyc1 = tyc / 1000000;
        let tyc2 = (tyc / 10) % 100000;
        let tyc3 = tyc % 10;
        writeln!(writer, "\"{}-{}-{}\",{}", tyc1, tyc2, tyc3, source_id)?;
    }

    Ok(())
}
