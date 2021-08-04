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
    collections::{hash_map::Entry, HashMap},
    io::{Read, Write},
};

use lazy_static::lazy_static;

use crate::votable::FieldInfo;

use super::{
    astro::{GaiaId, ProperMotion, SkyCoords},
    error::AppError,
    votable::{Ordinals, RecordAccessor, VotableReader, VotableRecord, VotableWriter},
};

mod hip;
mod tyc;

pub use hip::{hip_csv_crossmatch, HipStar};
pub use tyc::{tyc_csv_crossmatch, TycStar};

pub trait Crossmatchable<C> {
    fn score(&self, gaia_star: &C) -> f64;
}

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

impl Ordinals for GaiaOrdinals {
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
pub struct GaiaStar {
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
}

impl VotableRecord for GaiaStar {
    type Ordinals = GaiaOrdinals;
    type Id = GaiaId;

    fn from_accessor(accessor: &RecordAccessor, ordinals: &GaiaOrdinals) -> Result<Self, AppError> {
        let mut gaia_star = Self {
            source_id: GaiaId(
                accessor
                    .read_i64(ordinals.source_id)?
                    .ok_or(AppError::missing_id("source_id"))?,
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

    fn id(&self) -> Self::Id {
        self.source_id
    }

    fn fields() -> &'static [FieldInfo<Self>] {
        GAIA_FIELDS.as_ref()
    }
}

pub struct Crossmatcher<A, B>
where
    A: VotableRecord + Crossmatchable<B>,
    B: VotableRecord,
{
    source_stars: HashMap<A::Id, A>,
    crossmatch_stars: HashMap<B::Id, B>,
    matches: Vec<(A::Id, B::Id, f64)>,
}

impl<A, B> Crossmatcher<A, B>
where
    A: VotableRecord + Crossmatchable<B>,
    B: VotableRecord,
{
    pub fn new() -> Self {
        Self {
            source_stars: HashMap::new(),
            crossmatch_stars: HashMap::new(),
            matches: Vec::new(),
        }
    }

    pub fn add_reader(&mut self, mut reader: VotableReader<impl Read>) -> Result<(), AppError> {
        let source_ordinals = <A as VotableRecord>::Ordinals::from_reader(&reader)?;
        let crossmatch_ordinals = <B as VotableRecord>::Ordinals::from_reader(&reader)?;

        while let Some(accessor) = reader.read()? {
            let source_star = A::from_accessor(&accessor, &source_ordinals)?;
            let source_id = source_star.id();
            let crossmatch_star = B::from_accessor(&accessor, &crossmatch_ordinals)?;
            let crossmatch_id = crossmatch_star.id();
            let score = source_star.score(&crossmatch_star);

            self.source_stars.insert(source_id, source_star);
            self.crossmatch_stars.insert(crossmatch_id, crossmatch_star);
            self.matches.push((source_id, crossmatch_id, score));
        }

        Ok(())
    }

    pub fn finalize(mut self, writer: impl Write) -> Result<(), AppError> {
        println!("Found {} distinct source stars", self.source_stars.len());
        self.matches
            .sort_unstable_by(|(_, _, a), (_, _, b)| b.partial_cmp(a).unwrap());

        let mut writer = VotableWriter::new(writer)?;
        writer.write_fields::<A>()?;
        writer.write_fields::<B>()?;

        while let Some((source_id, crossmatch_id, _)) = self.matches.pop() {
            if let Entry::Occupied(entry) = self.crossmatch_stars.entry(crossmatch_id) {
                if let Some(source_star) = self.source_stars.remove(&source_id) {
                    let crossmatch_star = entry.remove();
                    writer.add_data(&source_star)?;
                    writer.add_data(&crossmatch_star)?;
                    writer.write_row()?;
                }
            }
        }

        writer.finish()?;

        println!("{} unmatched", self.source_stars.len());
        Ok(())
    }
}
