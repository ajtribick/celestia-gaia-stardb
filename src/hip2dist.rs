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
    fs::File,
    io::{self, BufReader, BufWriter, ErrorKind, Write},
    path::Path,
};

use crate::{astro::HipId, csv::CsvReader, error::AppError};

mod distributions;
mod estimate;

#[derive(Debug)]
struct PriorInfo {
    ggd_l: f64,
    ggd_alpha: f64,
    ggd_beta: f64,
    edsd_length: f64,
}

#[derive(Debug)]
struct HipInfo {
    hip: HipId,
    plx: f64,
    e_plx: f64,
    healpix: usize,
}

#[derive(Debug)]
struct DistanceInfo {
    hip: HipId,
    lower: f64,
    median: f64,
    upper: f64,
}

fn load_priors(path: impl AsRef<Path>) -> io::Result<Vec<PriorInfo>> {
    let file = File::open(path)?;
    let mut reader = CsvReader::new(BufReader::new(file))?;
    let healpix_col = reader.index("healpix").ok_or_else(|| io::Error::new(
        ErrorKind::InvalidData,
        "Missing healpix field",
    ))?;
    let ggd_l_col = reader.index("GGDrlen").ok_or_else(|| io::Error::new(
        ErrorKind::InvalidData,
        "Missing GGDrlen field",
    ))?;
    let ggd_alpha_col = reader.index("GGDalpha").ok_or_else(|| io::Error::new(
        ErrorKind::InvalidData,
        "Missing GGDalpha field",
    ))?;
    let ggd_beta_col = reader.index("GGDbeta").ok_or_else(|| io::Error::new(
        ErrorKind::InvalidData,
        "Missing field GGDbeta",
    ))?;
    let edsd_length_col = reader.index("EDSDrlen").ok_or_else(|| io::Error::new(
        ErrorKind::InvalidData,
        "Missing EDSDrlen field",
    ))?;

    let mut result = Vec::with_capacity(12288);
    while reader.next()?.is_some() {
        let healpix: usize = reader
            .field(healpix_col)
            .parse()
            .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?;
        if healpix != result.len() {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                "Prior file is not sequential",
            ));
        }
        let prior_info = PriorInfo {
            ggd_l: reader
                .field(ggd_l_col)
                .parse()
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
            ggd_alpha: reader
                .field(ggd_alpha_col)
                .parse()
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
            ggd_beta: reader
                .field(ggd_beta_col)
                .parse()
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
            edsd_length: reader
                .field(edsd_length_col)
                .parse()
                .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
        };
        result.push(prior_info);
    }

    Ok(result)
}

pub fn estimate_distances(
    prior_path: impl AsRef<Path>,
    hip_healpix_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
) -> Result<(), AppError> {
    let priors = load_priors(prior_path)?;
    let result = estimate::estimate(priors.into(), hip_healpix_path)?;

    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    writeln!(writer, "HIP,dist_low,dist_med,dist_high")?;
    for distance_info in result {
        writeln!(
            writer,
            "{},{},{},{}",
            distance_info.hip.0, distance_info.lower, distance_info.median, distance_info.upper
        )?;
    }

    Ok(())
}
