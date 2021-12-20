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

use std::f64;
use std::hash::Hash;
use std::ops;

use nalgebra::vector;

pub const MAS_TO_RADIANS: f64 = f64::consts::PI / (3600000.0 * 180.0);
pub const MAS_TO_DEG: f64 = 1.0 / 3600000.0;
pub const AU_IN_KM_YEAR_PER_S: f64 = 149597870.7 / (365.25 * 86400.0);

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
#[repr(transparent)]
pub struct GaiaId(pub i64);

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy, PartialOrd, Ord)]
#[repr(transparent)]
pub struct HipId(pub i32);

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
#[repr(transparent)]
pub struct TycId(pub i64);

pub trait Squarable: ops::Mul<Output = Self> + Sized + Copy {
    fn sqr(&self) -> Self {
        *self * *self
    }
}

impl Squarable for f32 {}
impl Squarable for f64 {}

#[derive(Debug, Clone, Copy)]
pub struct SkyCoords {
    pub ra: f64,
    pub dec: f64,
}

impl SkyCoords {
    pub fn ang_dist(&self, other: &Self) -> f64 {
        let (dec_sin_a, dec_cos_a) = self.dec.to_radians().sin_cos();
        let (dec_sin_b, dec_cos_b) = other.dec.to_radians().sin_cos();
        (dec_sin_a * dec_sin_b + dec_cos_a * dec_cos_b * (self.ra - other.ra).to_radians().cos())
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    }

    // basic linear motion, not accounting for light travel time
    pub fn apply_pm(
        &self,
        pm: &ProperMotion<f64>,
        rv: f64,
        parallax: f64,
        epoch1: f64,
        epoch2: f64,
    ) -> Self {
        let (ra_sin, ra_cos) = self.ra.to_radians().sin_cos();
        let (dec_sin, dec_cos) = self.dec.to_radians().sin_cos();
        let start = vector![ra_cos * dec_cos, ra_sin * dec_cos, dec_sin];

        let x = vector![-ra_sin, ra_cos, 0.0];
        let y = vector![-ra_cos * dec_sin, -ra_sin * dec_sin, dec_cos];

        let velocity = (pm.pm_ra * x + pm.pm_dec * y + rv * parallax / AU_IN_KM_YEAR_PER_S * start)
            * MAS_TO_RADIANS;
        let end = start + velocity * (epoch2 - epoch1);

        SkyCoords {
            ra: end.y.atan2(end.x).to_degrees(),
            dec: (end.z / end.norm()).asin().to_degrees(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ProperMotion<F> {
    pub pm_ra: F,
    pub pm_dec: F,
}

#[cfg(test)]
mod test {
    use super::*;

    // test cases from ADQL EPOCH_PROP_POS in Gaia archive

    #[test]
    fn apply_pm_zero_rv() {
        let coords = SkyCoords {
            ra: 132.0,
            dec: -38.0,
        };
        let pm = ProperMotion {
            pm_ra: 4200.0,
            pm_dec: 3200.0,
        };
        let (parallax, rv, epoch1, epoch2) = (850.0, 0.0, 1990.0, 2020.0);
        let expected = SkyCoords {
            ra: 2.30460952981089_f64.to_degrees(),
            dec: -0.662759549026617_f64.to_degrees(),
        };
        let actual = coords.apply_pm(&pm, rv, parallax, epoch1, epoch2);
        let ra_diff = (actual.ra - expected.ra).abs();
        let dec_diff = (actual.dec - expected.dec).abs();
        let max_diff = f64::max(ra_diff, dec_diff);
        assert!(
            max_diff < 1e-9,
            "Coordinate mismatch: expected {:?}, actual {:?}, diff ({}, {})",
            expected,
            actual,
            ra_diff,
            dec_diff,
        );
    }

    #[test]
    fn apply_pm_with_rv() {
        let coords = SkyCoords {
            ra: 132.0,
            dec: -38.0,
        };
        let pm = ProperMotion {
            pm_ra: 4200.0,
            pm_dec: 3200.0,
        };
        let (parallax, rv, epoch1, epoch2) = (850.0, 80.0, 1990.0, 2020.0);
        let expected = SkyCoords {
            ra: 2.30460791702764_f64.to_degrees(),
            dec: -0.662760518633619_f64.to_degrees(),
        };
        let actual = coords.apply_pm(&pm, rv, parallax, epoch1, epoch2);
        let ra_diff = (actual.ra - expected.ra).abs();
        let dec_diff = (actual.dec - expected.dec).abs();
        let max_diff = f64::max(ra_diff, dec_diff);
        assert!(
            max_diff < 1e-9,
            "Coordinate mismatch: expected {:?}, actual {:?}, diff ({}, {})",
            expected,
            actual,
            ra_diff,
            dec_diff
        );
    }
}
