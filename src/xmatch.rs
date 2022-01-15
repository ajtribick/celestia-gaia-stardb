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

use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::io::{Read, Write};

use lazy_static::lazy_static;

use crate::astro::{GaiaId, HipId, ProperMotion, SkyCoords, Squarable, TycId, MAS_TO_DEG};
use crate::error::AppError;
use crate::votable::{FieldInfo, RecordAccessor, VotableReader, VotableWriter};

#[derive(Debug)]
struct GaiaOrdinals {
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
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, AppError> {
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
struct HipOrdinals {
    pub hip: usize,
    pub hip_ra: usize,
    pub hip_dec: usize,
    pub hp_mag: usize,
    pub hip_parallax: Option<usize>,
}

impl HipOrdinals {
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, AppError> {
        Ok(Self {
            hip: reader.ordinal(b"hip")?,
            hip_ra: reader.ordinal(b"hip_ra")?,
            hip_dec: reader.ordinal(b"hip_dec")?,
            hp_mag: reader.ordinal(b"hp_mag")?,
            hip_parallax: reader.ordinal(b"hip_parallax").ok(),
        })
    }
}

#[derive(Debug)]
struct TycOrdinals {
    pub id_tycho: usize,
    pub ra_deg: usize,
    pub de_deg: usize,
    pub bt_mag: usize,
    pub vt_mag: usize,
    pub ep_ra1990: Option<usize>,
    pub ep_de1990: Option<usize>,
}

impl TycOrdinals {
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, AppError> {
        Ok(Self {
            id_tycho: reader.ordinal(b"id_tycho")?,
            ra_deg: reader.ordinal(b"tyc_ra")?,
            de_deg: reader.ordinal(b"tyc_dec")?,
            bt_mag: reader.ordinal(b"bt_mag")?,
            vt_mag: reader.ordinal(b"vt_mag")?,
            ep_ra1990: reader.ordinal(b"ep_ra1990").ok(),
            ep_de1990: reader.ordinal(b"ep_de1990").ok(),
        })
    }
}

#[derive(Debug)]
struct GaiaStar {
    pub source_id: GaiaId,
    pub coords: SkyCoords,
    pub parallax: f64,
    pub parallax_error: f32,
    pub rv: f32,
    pub epoch: f64,
    pub pm: ProperMotion<f64>,
    pub e_pm: ProperMotion<f32>,
    pub phot_g_mean_mag: f32,
    pub phot_bp_mean_mag: f32,
    pub phot_rp_mean_mag: f32,
    pub astrometric_params_solved: i16,
}

#[derive(Debug)]
struct CrossmatchStar {
    pub id: i64,
    pub coords: SkyCoords,
    pub hp_mag: f64,
    pub parallax: f64,
    pub bt_mag: f32,
    pub vt_mag: f32,
    pub epoch_ra: f32,
    pub epoch_dec: f32,
}

lazy_static! {
    static ref GAIA_FIELDS: [FieldInfo<GaiaStar>; 15] = [
        FieldInfo::long(
            "source_id",
            None,
            Some("meta.id"),
            "Unique source identifier (unique within a particular Data Release",
            |g| Some(g.source_id.0)
        ),
        FieldInfo::double(
            "ra",
            Some("deg"),
            Some("pos.eq.ra;meta.main"),
            "Right ascension",
            |g| g.coords.ra
        ),
        FieldInfo::double(
            "dec",
            Some("deg"),
            Some("pos.eq.dec;meta.main"),
            "Declination",
            |g| g.coords.dec
        ),
        FieldInfo::double(
            "parallax",
            Some("mas"),
            Some("pos.parallax.trig"),
            "Parallax",
            |g| g.parallax
        ),
        FieldInfo::float(
            "parallax_error",
            Some("mas"),
            Some("stat.error;pos.parallax.trig"),
            "Standard error of parallax",
            |g| g.parallax_error
        ),
        FieldInfo::float(
            "dr2_radial_velocity",
            Some("km.s**-1"),
            Some("spect.dopplerVeloc.opt;em.opt.I"),
            "Radial velocity from Gaia DR2",
            |g| g.rv
        ),
        FieldInfo::double(
            "ref_epoch",
            Some("yr"),
            Some("meta.ref;time.epoch"),
            "Reference epoch",
            |g| g.epoch
        ),
        FieldInfo::double(
            "pmra",
            Some("mas.yr**-1"),
            Some("pos.pm;pos.eq.ra"),
            "Proper motion in right ascension direction",
            |g| g.pm.pm_ra
        ),
        FieldInfo::float(
            "pmra_error",
            Some("mas.yr**-1"),
            Some("stat.error;pos.pm;pos.eq.ra"),
            "Standard error of proper motion in right ascension direction",
            |g| g.e_pm.pm_ra
        ),
        FieldInfo::double(
            "pmdec",
            Some("mas.yr**-1"),
            Some("pos.pm;pos.eq.dec"),
            "Proper motion in declination direction",
            |g| g.pm.pm_dec
        ),
        FieldInfo::float(
            "pmdec_error",
            Some("mas.yr**-1"),
            Some("stat.error;pos.pm;pos.eq.dec"),
            "Standard error of proper motion in declination direction",
            |g| g.e_pm.pm_dec
        ),
        FieldInfo::float(
            "phot_g_mean_mag",
            Some("mag"),
            Some("phot.mag;em.opt"),
            "G-band mean magnitude",
            |g| g.phot_g_mean_mag
        ),
        FieldInfo::float(
            "phot_bp_mean_mag",
            Some("mag"),
            Some("phot.mag;em.opt.B"),
            "Integrated BP mean magnitude",
            |g| g.phot_bp_mean_mag
        ),
        FieldInfo::float(
            "phot_rp_mean_mag",
            Some("mag"),
            Some("phot.mag;em.opt.R"),
            "Integrated RP mean magnitude",
            |g| g.phot_rp_mean_mag
        ),
        FieldInfo::short(
            "astrometric_params_solved",
            None,
            Some("meta.number"),
            "Which parameters have been solved for?",
            |g| Some(g.astrometric_params_solved)
        ),
    ];
    static ref CROSSMATCH_FIELDS: [FieldInfo<CrossmatchStar>; 8] = [
        FieldInfo::long("id", None, None, "Celestia identifier", |s| Some(s.id)),
        FieldInfo::double(
            "src_ra",
            Some("Angle[deg]"),
            Some("pos.eq.ra;meta.main"),
            "Right Ascension in ICRS, Ep=1991.25",
            |s| s.coords.ra
        ),
        FieldInfo::double(
            "src_dec",
            Some("Angle[deg]"),
            Some("pos.eq.dec;meta.main"),
            "Declination in ICRS, Ep=1991.25",
            |s| s.coords.dec
        ),
        FieldInfo::double(
            "hp_mag",
            Some("Magnitude[Mag]"),
            None,
            "Hipparcos magnitude",
            |s| s.hp_mag
        ),
        FieldInfo::float(
            "bt_mag",
            Some("Magnitude[Mag]"),
            Some("em.opt.B;phot.mag"),
            "Tycho-2 BT magnitude",
            |s| s.bt_mag
        ),
        FieldInfo::float(
            "vt_mag",
            Some("Magnitude[Mag]"),
            Some("em.opt.V;phot.mag"),
            "Tycho-2 VT magnitude",
            |s| s.bt_mag
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

impl GaiaStar {
    fn from_accessor(accessor: &RecordAccessor, ordinals: &GaiaOrdinals) -> Result<Self, AppError> {
        let mut gaia_star = Self {
            source_id: GaiaId(
                accessor
                    .read_i64(ordinals.source_id)?
                    .ok_or_else(|| AppError::missing_id("source_id"))?,
            ),
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.ra)?,
                dec: accessor.read_f64(ordinals.dec)?,
            },
            parallax: accessor.read_f64(ordinals.parallax)?,
            parallax_error: accessor.read_f32(ordinals.parallax_error)?,
            rv: accessor.read_f32(ordinals.dr2_radial_velocity)?,
            epoch: accessor.read_f64(ordinals.ref_epoch)?,
            pm: ProperMotion {
                pm_ra: accessor.read_f64(ordinals.pm_ra)?,
                pm_dec: accessor.read_f64(ordinals.pm_dec)?,
            },
            e_pm: ProperMotion {
                pm_ra: accessor.read_f32(ordinals.pm_ra_error)?,
                pm_dec: accessor.read_f32(ordinals.pm_dec_error)?,
            },
            phot_g_mean_mag: accessor.read_f32(ordinals.phot_g_mean_mag)?,
            phot_bp_mean_mag: accessor.read_f32(ordinals.phot_bp_mean_mag)?,
            phot_rp_mean_mag: accessor.read_f32(ordinals.phot_rp_mean_mag)?,
            astrometric_params_solved: accessor
                .read_i16(ordinals.astrometric_params_solved)?
                .unwrap_or(0),
        };

        // magnitude correction from "Gaia Early Data Release 3: Summary of the contents and survey propreties" (erratum)
        if gaia_star.astrometric_params_solved == 31 || gaia_star.phot_g_mean_mag < 13.0 {
            return Ok(gaia_star);
        }

        let mut bp_rp = (gaia_star.phot_bp_mean_mag - gaia_star.phot_rp_mean_mag).clamp(0.25, 3.0);
        if bp_rp.is_nan() {
            bp_rp = 0.25;
        }

        if gaia_star.phot_g_mean_mag < 16.0 {
            gaia_star.phot_g_mean_mag -=
                2.5 * (1.00876 + bp_rp * (-0.02540 + bp_rp * (0.01747 - bp_rp * 0.00277))).log10();
        } else {
            gaia_star.phot_g_mean_mag -=
                2.5 * (1.00525 + bp_rp * (-0.02323 + bp_rp * (0.01740 - bp_rp * 0.00253))).log10();
        }

        Ok(gaia_star)
    }
}

fn g_vt(bp_rp: f32) -> f32 {
    -0.01077 + bp_rp * (-0.0682 + bp_rp * (-0.2387 + bp_rp * 0.02342))
}

fn g_bt(bp_rp: f32) -> f32 {
    -0.004288
        + bp_rp
            * (-0.8547 + bp_rp * (0.1244 + bp_rp * (-0.9085 + bp_rp * (0.4843 + bp_rp * -0.06814))))
}

impl CrossmatchStar {
    pub fn from_hip_accessor(
        accessor: &RecordAccessor,
        ordinals: &HipOrdinals,
    ) -> Result<Self, AppError> {
        Ok(Self {
            id: accessor
                .read_i32(ordinals.hip)?
                .ok_or_else(|| AppError::missing_id("hip"))? as i64,
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.hip_ra)?,
                dec: accessor.read_f64(ordinals.hip_dec)?,
            },
            hp_mag: accessor
                .read_f64(ordinals.hp_mag)
                .or_else(|_| accessor.read_f32(ordinals.hp_mag).map(|h| h as f64))?,
            parallax: ordinals.hip_parallax.map_or(Ok(f64::NAN), |n| accessor.read_f64(n))?,
            bt_mag: f32::NAN,
            vt_mag: f32::NAN,
            epoch_ra: 1.25,
            epoch_dec: 1.25,
        })
    }

    pub fn from_tyc_accessor(
        accessor: &RecordAccessor,
        ordinals: &TycOrdinals,
    ) -> Result<Self, AppError> {
        Ok(Self {
            id: accessor
                .read_i64(ordinals.id_tycho)?
                .ok_or_else(|| AppError::missing_id("id_tycho"))?,
            coords: SkyCoords {
                ra: accessor.read_f64(ordinals.ra_deg)?,
                dec: accessor.read_f64(ordinals.de_deg)?,
            },
            hp_mag: f64::NAN,
            parallax: f64::NAN,
            bt_mag: accessor
                .read_f32(ordinals.bt_mag)
                .or_else(|_| accessor.read_f64(ordinals.bt_mag).map(|x| x as f32))?,
            vt_mag: accessor
                .read_f32(ordinals.vt_mag)
                .or_else(|_| accessor.read_f64(ordinals.vt_mag).map(|x| x as f32))?,
            epoch_ra: ordinals
                .ep_ra1990
                .map_or(Ok(1.25), |ord| accessor.read_f32(ord))?,
            epoch_dec: ordinals
                .ep_de1990
                .map_or(Ok(1.25), |ord| accessor.read_f32(ord))?,
        })
    }

    pub fn merge(&mut self, other: &CrossmatchStar) {
        if self.id <= other.id {
            if f64::is_nan(self.hp_mag) {
                self.hp_mag = other.hp_mag;
            }
            if !f32::is_nan(other.bt_mag) {
                self.bt_mag = other.bt_mag;
            }
            if !f32::is_nan(other.vt_mag) {
                self.vt_mag = other.vt_mag;
            }
        } else {
            self.id = other.id;
            self.coords = other.coords;
            if !f64::is_nan(other.hp_mag) {
                self.hp_mag = other.hp_mag;
            }
            if f32::is_nan(self.bt_mag) {
                self.bt_mag = other.bt_mag;
            }
            if f32::is_nan(self.vt_mag) {
                self.vt_mag = other.vt_mag;
            }
            self.epoch_ra = other.epoch_ra;
            self.epoch_dec = other.epoch_dec;
        }
    }

    pub fn score(&self, gaia_star: &GaiaStar) -> f64 {
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
            "ID {}, Gaia {}\n{:?}\n{:?}",
            self.id,
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

        let plx_diff = if self.parallax.is_nan() || gaia_star.parallax.is_nan() {
            100.0
        } else {
            (self.parallax - gaia_star.parallax) / gaia_star.parallax_error as f64
        };

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

        let mag_diff = if self.hp_mag.is_nan() {
            let mag_diff = if self.bt_mag.is_nan() {
                gaia_star.phot_g_mean_mag - self.vt_mag - g_vt(bp_rp)
            } else if self.vt_mag.is_nan() {
                gaia_star.phot_g_mean_mag - self.bt_mag - g_bt(bp_rp)
            } else {
                gaia_star.phot_g_mean_mag
                    - (self.vt_mag + self.bt_mag + g_vt(bp_rp) + g_bt(bp_rp)) / 2.0
            };

            mag_diff as f64 / 0.06
        } else {
            let mag_diff = self.hp_mag - (0.91 * bp_mag + 0.09 * rp_mag) as f64;
            mag_diff / 0.1
        };

        assert!(!mag_diff.is_nan());

        pm_diff.recip() + plx_diff.sqr().recip() + mag_diff.sqr().recip() + (dist / MAS_TO_DEG).sqr().recip()
    }
}

pub struct Crossmatcher {
    tyc2hip: HashMap<TycId, HipId>,
    crossmatch_stars: HashMap<i64, CrossmatchStar>,
    gaia_stars: HashMap<GaiaId, GaiaStar>,
    matches: Vec<(i64, GaiaId, f64)>,
}

impl Crossmatcher {
    pub fn new(tyc2hip: HashMap<TycId, HipId>) -> Self {
        Self {
            tyc2hip,
            crossmatch_stars: HashMap::new(),
            gaia_stars: HashMap::new(),
            matches: Vec::new(),
        }
    }

    pub fn add_hip(&mut self, mut reader: VotableReader<impl Read>) -> Result<(), AppError> {
        let hip_ordinals = HipOrdinals::from_reader(&reader)?;
        let gaia_ordinals = GaiaOrdinals::from_reader(&reader)?;

        while let Some(accessor) = reader.read()? {
            let hip_star = CrossmatchStar::from_hip_accessor(&accessor, &hip_ordinals)?;
            let hip_id = hip_star.id;
            let gaia_star = GaiaStar::from_accessor(&accessor, &gaia_ordinals)?;
            let gaia_id = gaia_star.source_id;
            let score = hip_star.score(&gaia_star);

            self.crossmatch_stars
                .entry(hip_id)
                .and_modify(|e| e.merge(&hip_star))
                .or_insert(hip_star);
            self.gaia_stars.insert(gaia_id, gaia_star);
            self.matches.push((hip_id, gaia_id, score));
        }

        Ok(())
    }

    pub fn add_tyc(&mut self, mut reader: VotableReader<impl Read>) -> Result<(), AppError> {
        let tyc_ordinals = TycOrdinals::from_reader(&reader)?;
        let gaia_ordinals = GaiaOrdinals::from_reader(&reader)?;

        while let Some(accessor) = reader.read()? {
            let tyc_star = CrossmatchStar::from_tyc_accessor(&accessor, &tyc_ordinals)?;
            let tyc_id = self
                .tyc2hip
                .get(&TycId(tyc_star.id))
                .map_or(tyc_star.id, |h| h.0 as i64);
            let gaia_star = GaiaStar::from_accessor(&accessor, &gaia_ordinals)?;
            let gaia_id = gaia_star.source_id;
            let score = tyc_star.score(&gaia_star);

            self.crossmatch_stars
                .entry(tyc_id)
                .and_modify(|e| e.merge(&tyc_star))
                .or_insert(tyc_star);
            self.gaia_stars.insert(gaia_id, gaia_star);
            self.matches.push((tyc_id, gaia_id, score));
        }

        Ok(())
    }

    pub fn finalize(mut self, writer: impl Write) -> Result<(), AppError> {
        println!(
            "Found {} distinct source stars",
            self.crossmatch_stars.len()
        );
        self.matches
            .sort_unstable_by(|(_, _, a), (_, _, b)| a.partial_cmp(b).unwrap());

        let mut writer = VotableWriter::new(writer)?;
        writer.write_fields(CROSSMATCH_FIELDS.as_ref())?;
        writer.write_fields(GAIA_FIELDS.as_ref())?;

        while let Some((source_id, crossmatch_id, _)) = self.matches.pop() {
            if let Entry::Occupied(entry) = self.gaia_stars.entry(crossmatch_id) {
                if let Some(source_star) = self.crossmatch_stars.remove(&source_id) {
                    let gaia_star = entry.remove();
                    writer.add_data(&source_star, CROSSMATCH_FIELDS.as_ref())?;
                    writer.add_data(&gaia_star, GAIA_FIELDS.as_ref())?;
                    writer.write_row()?;
                }
            }
        }

        writer.finish()?;

        println!("{} unmatched", self.crossmatch_stars.len());
        Ok(())
    }
}
