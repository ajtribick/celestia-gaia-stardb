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

// Solve mode of exponentially decaying space distribution
// details given in Bailer-Jones (2015)

use std::f64::consts::{FRAC_PI_3, FRAC_PI_6};

use super::{HipInfo, PriorInfo};

const FRAC_2_SQRT_3: f64 = 1.1547005383792515; // 2/sqrt(3)
const FRAC_3_SQRT_3_2: f64 = 2.598076211353316; // 3*sqrt(3)/2

#[allow(clippy::many_single_char_names)]
pub(super) fn edsd_mode(edsd_length: f64, plx: f64, e_plx: f64) -> f64 {
    let e_plx2 = (e_plx * e_plx).recip();
    let a = -2.0 * edsd_length;
    let b = plx * e_plx2 * edsd_length * 1000.0;
    let c = -1000000.0 * edsd_length * e_plx2;

    // solve cubic
    let p = b - a * a / 3.0;
    let q = 2.0 * a * a * a / 27.0 - a * b / 3.0 + c;
    let delta = 0.25 * q * q + p * p * p / 27.0;
    assert!(delta.is_finite(), "Delta is not finite");
    assert!(delta != 0.0, "Delta is zero");

    let a_3 = a / 3.0;
    let mut result = if delta > 0.0 {
        let r = -q * 0.5;
        let root_delta = delta.sqrt();
        (r + root_delta).cbrt() + (r - root_delta).cbrt() - a_3
    } else {
        let p_sqrt = (-p).sqrt();
        let x = FRAC_2_SQRT_3 * p_sqrt;
        let y = (FRAC_3_SQRT_3_2 * q / (p_sqrt * p_sqrt * p_sqrt)).asin() / 3.0;
        IntoIterator::into_iter([
            x * y.sin() - a_3,
            -x * (y + FRAC_PI_3).sin() - a_3,
            x * (y + FRAC_PI_6).cos() - a_3,
        ])
        .filter(|r| *r >= 0.0)
        .reduce(f64::min)
        .unwrap()
    };

    // use Newton-Raphson to get better approximation
    let a2 = 2.0 * a;
    loop {
        let offset =
            (c + result * (b + result * (a + result))) / (b + result * (a2 + result * 3.0));
        result -= offset;
        if offset.abs() < 1e-6 {
            return result;
        }
    }
}

pub(super) fn geometric_posterior(hip_info: &HipInfo, prior: &PriorInfo) -> impl Fn(f64) -> f64 {
    // Create geometric posterior, prior from Bailer-Jones et al. (2021)
    let plx = hip_info.plx / 1000.0;
    let l = prior.ggd_l;
    let alpha = prior.ggd_alpha;
    let beta = prior.ggd_beta;
    let x = 500000.0 / (hip_info.e_plx * hip_info.e_plx);
    move |r| {
        let diff = plx - r.recip();
        beta * r.ln() - (r / l).powf(alpha) - x * diff * diff
    }
}
