use std::f64;

use nalgebra::Vector3;

const MAS_TO_RADIANS: f64 = f64::consts::PI / (3600000.0 * 180.0);
const AU_IN_KM_YEAR_PER_S: f64 = 149597870.7 / (365.25 * 86400.0);

// basic linear motion, not accounting for light travel time
pub fn apply_pm(
    ra: f64,
    dec: f64,
    pm_ra: f64,
    pm_dec: f64,
    rv: f64,
    parallax: f64,
    epoch1: f64,
    epoch2: f64,
) -> (f64, f64) {
    let (ra_sin, ra_cos) = ra.to_radians().sin_cos();
    let (dec_sin, dec_cos) = dec.to_radians().sin_cos();
    let start = Vector3::new(ra_cos * dec_cos, ra_sin * dec_cos, dec_sin);

    // get coordinate vectors
    let z = start.normalize();
    let x = Vector3::new(-ra_sin, ra_cos, 0.0);
    let y = z.cross(&x);

    let velocity = pm_ra * MAS_TO_RADIANS * x
        + pm_dec * MAS_TO_RADIANS * y
        + rv * parallax * MAS_TO_RADIANS / AU_IN_KM_YEAR_PER_S * z;

    let end = start + velocity * (epoch2 - epoch1);

    (
        end.y.atan2(end.x).to_degrees(),
        (end.z / end.norm()).asin().to_degrees(),
    )
}

#[cfg(test)]
mod test {
    use super::*;

    // test cases from ADQL EPOCH_PROP_POS in Gaia archive

    #[test]
    fn apply_pm_zero_rv() {
        let (ra, dec, parallax, pm_ra, pm_dec, rv, epoch1, epoch2) =
            (132.0, -38.0, 850.0, 4200.0, 3200.0, 0.0, 1990.0, 2020.0);
        let ra_expected = 2.30460952981089_f64.to_degrees();
        let dec_expected = -0.662759549026617_f64.to_degrees();
        let (ra_actual, dec_actual) =
            apply_pm(ra, dec, pm_ra, pm_dec, rv, parallax, epoch1, epoch2);
        let ra_diff = (ra_actual - ra_expected).abs();
        let dec_diff = (dec_actual - dec_expected).abs();
        let max_diff = f64::max(ra_diff, dec_diff);
        assert!(
            max_diff < 1e-9,
            "Coordinate mismatch: expected ({}, {}), actual ({}, {}), diff ({}, {})",
            ra_expected,
            dec_expected,
            ra_actual,
            dec_actual,
            ra_diff,
            dec_diff
        );
    }

    #[test]
    fn apply_pm_with_rv() {
        let (ra, dec, parallax, pm_ra, pm_dec, rv, epoch1, epoch2) =
            (132.0, -38.0, 850.0, 4200.0, 3200.0, 80.0, 1990.0, 2020.0);
        let ra_expected = 2.30460791702764_f64.to_degrees();
        let dec_expected = -0.662760518633619_f64.to_degrees();
        let (ra_actual, dec_actual) =
            apply_pm(ra, dec, pm_ra, pm_dec, rv, parallax, epoch1, epoch2);
        let ra_diff = (ra_actual - ra_expected).abs();
        let dec_diff = (dec_actual - dec_expected).abs();
        let max_diff = f64::max(ra_diff, dec_diff);
        assert!(
            max_diff < 1e-9,
            "Coordinate mismatch: expected ({}, {}), actual ({}, {}), diff ({}, {})",
            ra_expected,
            dec_expected,
            ra_actual,
            dec_actual,
            ra_diff,
            dec_diff
        );
    }
}
