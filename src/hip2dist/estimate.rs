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

use std::f64::consts::FRAC_1_SQRT_2;
use std::fs::File;
use std::io::{self, BufRead, BufReader, ErrorKind};
use std::path::Path;
use std::sync::Arc;
use std::thread::{self, JoinHandle};

use crossbeam_channel::{Receiver, Sender};
use lazy_static::lazy_static;
use rand::distributions::{Distribution, Uniform};
use statrs::{distribution::Normal, function::erf};

use super::distributions::{edsd_mode, geometric_posterior};
use super::{DistanceInfo, HipInfo, PriorInfo};

use crate::astro::HipId;
use crate::csv::CsvReader;
use crate::error::AppError;

const MCMC_SAMPLES: usize = 50000;
const BURN_IN_SAMPLES: usize = MCMC_SAMPLES / 10;

lazy_static! {
    static ref LOWER_POS: f64 = 0.5 * (1.0 + erf::erf(-FRAC_1_SQRT_2));
    static ref UPPER_POS: f64 = 0.5 * (1.0 + erf::erf(FRAC_1_SQRT_2));
}

fn percentile(samples: &[f64], position: f64) -> f64 {
    // Find percentile in sorted set, interpolating if necessary
    assert!(!samples.is_empty());
    assert!((0.0..=1.0).contains(&position), "Position out of range");
    let position = (samples.len() - 1) as f64 * position;
    let index = position as usize;
    let weight = position.fract();
    if weight == 0.0 {
        samples[index]
    } else {
        (1.0 - weight) * samples[index] + weight * samples[index + 1]
    }
}

struct Parser<B: BufRead + Send + 'static> {
    reader: CsvReader<B>,
    hip_col: usize,
    plx_col: usize,
    e_plx_col: usize,
    healpix_col: usize,
    sender: Sender<HipInfo>,
}

impl<B: BufRead + Send> Parser<B> {
    pub fn new(reader: CsvReader<B>, sender: Sender<HipInfo>) -> io::Result<Self> {
        let hip_col = reader
            .index("HIP")
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Missing HIP field"))?;
        let plx_col = reader
            .index("Plx")
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Missing Plx field"))?;
        let e_plx_col = reader
            .index("e_Plx")
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Missing e_Plx field"))?;
        let healpix_col = reader
            .index("healpix")
            .ok_or_else(|| io::Error::new(ErrorKind::InvalidData, "Missing healpix field"))?;
        Ok(Self {
            reader,
            hip_col,
            plx_col,
            e_plx_col,
            healpix_col,
            sender,
        })
    }

    fn process(&mut self) -> Result<usize, AppError> {
        let mut processed = 0;
        while self.reader.next()?.is_some() {
            let hip_info = HipInfo {
                hip: HipId(
                    self.reader
                        .field(self.hip_col)
                        .parse()
                        .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
                ),
                plx: self
                    .reader
                    .field(self.plx_col)
                    .parse()
                    .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
                e_plx: self
                    .reader
                    .field(self.e_plx_col)
                    .parse()
                    .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
                healpix: self
                    .reader
                    .field(self.healpix_col)
                    .parse()
                    .map_err(|e| io::Error::new(ErrorKind::InvalidData, e))?,
            };

            self.sender
                .send(hip_info)
                .map_err(|e| AppError::Other(Box::new(e)))?;
            processed += 1;
        }

        Ok(processed)
    }

    pub fn launch(mut self) -> JoinHandle<Result<usize, AppError>> {
        thread::spawn(move || self.process())
    }
}

struct Calculator {
    priors: Arc<[PriorInfo]>,
    receiver: Receiver<HipInfo>,
    sender: Sender<DistanceInfo>,
    burn_in_samples: usize,
    mcmc_samples: usize,
}

impl Calculator {
    pub fn new(
        priors: Arc<[PriorInfo]>,
        receiver: Receiver<HipInfo>,
        sender: Sender<DistanceInfo>,
        burn_in_samples: usize,
        mcmc_samples: usize,
    ) -> Self {
        Self {
            priors,
            receiver,
            sender,
            burn_in_samples,
            mcmc_samples,
        }
    }

    fn process(&self) -> Result<(), AppError> {
        let mut rng = rand::thread_rng();
        let udist: Uniform<f64> = Uniform::new(0.0, 1.0);
        let mut samples = vec![0.0; self.mcmc_samples];

        while let Ok(hip_info) = self.receiver.recv() {
            let prior = &self.priors[hip_info.healpix];

            // Initialize MCMC, setup given in Bailer-Jones et al. (2021)
            let mode = edsd_mode(prior.edsd_length, hip_info.plx, hip_info.e_plx);
            let p_dist = geometric_posterior(&hip_info, prior);
            let mut r = mode;
            let mut log_p = p_dist(r);

            let step_size = 0.75 * r * f64::min((hip_info.e_plx / hip_info.plx).abs(), 1.0 / 3.0);
            let step_dist = Normal::new(0.0, step_size).unwrap();

            // Burn-in
            for _ in 0..self.burn_in_samples {
                let r_new = r + step_dist.sample(&mut rng);
                if r_new <= 0.0 {
                    continue;
                }

                let log_p_new = p_dist(r_new);
                if log_p_new >= log_p || udist.sample(&mut rng).ln() <= log_p_new - log_p {
                    r = r_new;
                    log_p = log_p_new;
                }
            }

            // Collect samples
            for sample in samples.iter_mut() {
                let r_new = r + step_dist.sample(&mut rng);
                if r_new > 0.0 {
                    let log_p_new = p_dist(r_new);
                    if log_p_new >= log_p || udist.sample(&mut rng).ln() <= log_p_new - log_p {
                        r = r_new;
                        log_p = log_p_new;
                    }
                }
                *sample = r;
            }

            samples.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            let distance_info = DistanceInfo {
                hip: hip_info.hip,
                lower: percentile(&samples, *LOWER_POS),
                median: percentile(&samples, 0.5),
                upper: percentile(&samples, *UPPER_POS),
            };

            self.sender
                .send(distance_info)
                .map_err(|e| AppError::Other(Box::new(e)))?;
        }
        Ok(())
    }

    pub fn launch(self) -> JoinHandle<Result<(), AppError>> {
        thread::spawn(move || self.process())
    }
}

pub(super) fn estimate(
    priors: Arc<[PriorInfo]>,
    path: impl AsRef<Path>,
) -> Result<Vec<DistanceInfo>, AppError> {
    let file = File::open(path)?;
    let reader = CsvReader::new(BufReader::new(file))?;

    let num_calculators = usize::max(num_cpus::get().saturating_sub(2), 1);
    let (info_sender, info_receiver) = crossbeam_channel::bounded(num_calculators * 2);
    let (dist_sender, dist_receiver) = crossbeam_channel::bounded(num_calculators * 2);

    let parser = Parser::new(reader, info_sender)?.launch();
    let mut calculators = Vec::with_capacity(num_calculators);
    for _ in 0..num_calculators {
        let calculator = Calculator::new(
            priors.clone(),
            info_receiver.clone(),
            dist_sender.clone(),
            MCMC_SAMPLES,
            BURN_IN_SAMPLES,
        );
        calculators.push(calculator.launch());
    }

    drop(info_receiver);
    drop(dist_sender);

    let mut distances = Vec::with_capacity(120000);
    while let Ok(distance_info) = dist_receiver.recv() {
        distances.push(distance_info);
    }

    let num_processed = parser.join().map_err(AppError::Thread)??;
    let result = if distances.len() == num_processed {
        distances.sort_unstable_by(|a, b| a.hip.cmp(&b.hip));
        Ok(distances)
    } else {
        Err(AppError::Thread(Box::new(
            "Mismatch between processed count and results",
        )))
    };

    calculators.into_iter().fold(result, |res, t| {
        t.join()
            .unwrap_or_else(|e| Err(AppError::Thread(e)))
            .and(res)
    })
}
