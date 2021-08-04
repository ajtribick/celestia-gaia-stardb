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

use crate::{
    astro::{HipId, ProperMotion, SkyCoords, Squarable, MAS_TO_DEG},
    error::AppError,
    votable::{FieldInfo, Ordinals, RecordAccessor, VotableReader, VotableRecord},
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
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, AppError> {
        Ok(Self {
            hip: reader.ordinal(b"hip")?,
            hip_ra: reader.ordinal(b"hip_ra")?,
            hip_dec: reader.ordinal(b"hip_dec")?,
            hp_mag: reader.ordinal(b"hp_mag")?,
        })
    }
}

#[derive(Debug)]
pub struct HipStar {
    pub hip: HipId,
    pub coords: SkyCoords,
    pub hp_mag: f64,
}

lazy_static! {
    static ref HIP_FIELDS: [FieldInfo<HipStar>; 4] = [
        FieldInfo::int("hip", None, None, "Hipparcos identifier", |h| Some(h.hip.0)),
        FieldInfo::double(
            "hip_ra",
            Some("deg"),
            Some("pos.eq.ra;meta.main"),
            "Right Ascension in ICRS, Ep=1991.25",
            |g| g.coords.ra
        ),
        FieldInfo::double(
            "hip_dec",
            Some("deg"),
            Some("pos.eq.dec;meta.main"),
            "Declination in ICRS, Ep=1991.25",
            |g| g.coords.dec
        ),
        FieldInfo::double("hp_mag", Some("mag"), None, "Hipparcos magnitude", |g| g
            .hp_mag),
    ];
}

impl VotableRecord for HipStar {
    type Ordinals = HipOrdinals;
    type Id = HipId;

    fn from_accessor(accessor: &RecordAccessor, ordinals: &HipOrdinals) -> Result<Self, AppError> {
        Ok(Self {
            hip: HipId(
                accessor
                    .read_i32(ordinals.hip)?
                    .ok_or(AppError::missing_id("hip"))?,
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

    fn fields() -> &'static [FieldInfo<Self>] {
        HIP_FIELDS.as_ref()
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
        let propagated =
            gaia_star
                .coords
                .apply_pm(&pm, rv as f64, parallax, gaia_star.epoch, HIP_EPOCH);
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
        let mag_diff = self.hp_mag - (0.91 * bp_mag + 0.09 * rp_mag) as f64;
        assert!(!mag_diff.is_nan());

        pm_diff + (mag_diff / 0.1).sqr() + (dist / MAS_TO_DEG).sqr()
    }
}

pub fn hip_csv_crossmatch(
    mut reader: VotableReader<impl Read>,
    writer: &mut impl Write,
) -> Result<(), AppError> {
    let mut records = Vec::with_capacity(117955);
    let hip_ordinal = reader.ordinal(b"hip")?;
    let source_id_ordinal = reader.ordinal(b"source_id")?;
    writeln!(writer, "hip,source_id")?;
    while let Some(record) = reader.read()? {
        let hip = record
            .read_i32(hip_ordinal)?
            .ok_or(AppError::missing_id("hip"))?;
        let source_id = record
            .read_i64(source_id_ordinal)?
            .ok_or(AppError::missing_id("source_id"))?;
        records.push((hip, source_id));
    }

    records.sort_unstable_by(|(a, _), (b, _)| a.cmp(b));
    for (hip, source_id) in records {
        writeln!(writer, "{},{}", hip, source_id)?;
    }

    Ok(())
}
