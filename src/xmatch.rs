use std::{collections::HashSet, hash::Hash, io::Read};

use super::{
    astro::{ProperMotion, SkyCoords},
    error::Error,
    votable::{Ordinals, RecordAccessor, VotableReader, VotableRecord},
};

mod hip;
mod tyc;

pub use hip::HipStar;
pub use tyc::TycStar;

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
    fn from_reader(reader: &VotableReader<impl Read>) -> Result<Self, Error> {
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

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
pub struct GaiaId(pub i64);

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

impl VotableRecord for GaiaStar {
    type Ordinals = GaiaOrdinals;
    type Id = GaiaId;

    fn from_accessor(accessor: &RecordAccessor, ordinals: &GaiaOrdinals) -> Result<Self, Error> {
        Ok(Self {
            source_id: GaiaId(
                accessor
                    .read_i64(ordinals.source_id)?
                    .ok_or(Error::missing_id("source_id"))?,
            ),
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

    fn id(&self) -> Self::Id {
        self.source_id
    }
}

pub struct Crossmatcher<A, B>
where
    A: VotableRecord + Crossmatchable<B>,
    B: VotableRecord,
{
    source_ids: HashSet<A::Id>,
    crossmatch_ids: HashSet<B::Id>,
    matches: Vec<(A::Id, B::Id, f64)>,
}

impl<A, B> Crossmatcher<A, B>
where
    A: VotableRecord + Crossmatchable<B>,
    B: VotableRecord,
{
    pub fn new() -> Self {
        Self {
            source_ids: HashSet::new(),
            crossmatch_ids: HashSet::new(),
            matches: Vec::new(),
        }
    }

    pub fn add_reader(&mut self, mut reader: VotableReader<impl Read>) -> Result<(), Error> {
        let source_ordinals = <A as VotableRecord>::Ordinals::from_reader(&reader)?;
        let crossmatch_ordinals = <B as VotableRecord>::Ordinals::from_reader(&reader)?;

        while let Some(accessor) = reader.read()? {
            let source_star = A::from_accessor(&accessor, &source_ordinals)?;
            let crossmatch_star = B::from_accessor(&accessor, &crossmatch_ordinals)?;
            let score = source_star.score(&crossmatch_star);

            self.source_ids.insert(source_star.id());
            self.crossmatch_ids.insert(crossmatch_star.id());
            self.matches
                .push((source_star.id(), crossmatch_star.id(), score));
        }

        Ok(())
    }

    pub fn finalize(&mut self) {
        println!("Found {} distinct source stars", self.source_ids.len());
        self.matches
            .sort_by(|(_, _, a), (_, _, b)| b.partial_cmp(a).unwrap());
        while let Some((hip, gaia, _)) = self.matches.pop() {
            if self.crossmatch_ids.contains(&gaia) && self.source_ids.remove(&hip) {
                self.crossmatch_ids.remove(&gaia);
            }
        }

        println!("{} unmatched", self.source_ids.len());
    }
}
