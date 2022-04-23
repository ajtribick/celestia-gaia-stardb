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

use std::borrow::Cow;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use tinyvec::ArrayVec;

use crate::astro::{HipId, TycId};
use crate::error::AppError;
use crate::votable::VotableReader;

#[derive(Debug)]
struct TycComponent(TycId, ArrayVec<[u8; 2]>);

#[derive(Debug)]
struct TycHipMap(HashMap<HipId, TycComponent>);

impl TycHipMap {
    fn new() -> Self {
        Self(HashMap::new())
    }

    fn add(&mut self, hip: HipId, tyc: TycId, cmp: ArrayVec<[u8; 2]>) {
        match self.0.entry(hip) {
            Entry::Vacant(v) => {
                v.insert(TycComponent(tyc, cmp));
            }
            Entry::Occupied(mut e) => {
                if cmp < e.get().1 {
                    e.insert(TycComponent(tyc, cmp));
                }
            }
        }
    }
}

pub fn load_tyc2hip(
    gaia_path: &Path,
    vizier_path: &Path,
) -> Result<HashMap<TycId, HipId>, AppError> {
    let mut tdsc_path = PathBuf::from(gaia_path);
    tdsc_path.push("tyc2tdsc_hip_xmatch.vot.gz");

    let mut hip2tyc = TycHipMap::new();

    load_tyc2tdsc_hip(&tdsc_path, &mut hip2tyc)?;
    load_tyc_suppl1_hip(&vizier_path, &mut hip2tyc)?;

    Ok(hip2tyc.0.into_iter().map(|(h, tc)| (tc.0, h)).collect())
}

fn load_tyc2tdsc_hip(path: &Path, hip2tyc: &mut TycHipMap) -> Result<(), AppError> {
    let file = File::open(path)?;
    let mut reader = VotableReader::new(file)?;

    let id_tycho_col = reader.ordinal(b"id_tycho")?;
    let hip_col = reader.ordinal(b"hip")?;
    let comp_col = reader.ordinal(b"cmp")?;

    while let Some(accessor) = reader.read()? {
        let id_tycho = TycId(
            accessor
                .read_i64(id_tycho_col)?
                .ok_or_else(|| AppError::missing_id("id_tycho"))?,
        );
        let hip = HipId(
            accessor
                .read_i32(hip_col)?
                .ok_or_else(|| AppError::missing_id("hip"))?,
        );
        let cmp = accessor.read_string::<2>(comp_col)?;

        hip2tyc.add(hip, id_tycho, cmp);
    }

    Ok(())
}

struct Suppl1Fields {
    tyc1: (usize, usize),
    tyc2: (usize, usize),
    tyc3: (usize, usize),
    hip: (usize, usize),
    cmp: (usize, usize),
}

fn load_tyc_suppl1_hip(path: &Path, hip2tyc: &mut TycHipMap) -> Result<(), AppError> {
    let mut readme_path = path.to_path_buf();
    readme_path.push("tyc2.readme");
    let fields = read_tyc2_suppl1_readme(&readme_path)?;

    let mut suppl1_path = path.to_path_buf();
    suppl1_path.push("tyc2suppl_1.dat.gz");

    let file = File::open(suppl1_path)?;
    let decoder = GzDecoder::new(file);
    let reader = BufReader::new(decoder);
    for line_result in reader.lines() {
        let line = line_result?;
        if line.len() < fields.hip.1 {
            continue;
        }
        let hip_str = line[fields.hip.0..fields.hip.1].trim();
        if hip_str.is_empty() {
            continue;
        }
        let hip = HipId(hip_str.parse()?);
        let tyc1: i64 = line[fields.tyc1.0..fields.tyc1.1].parse()?;
        let tyc2: i64 = line[fields.tyc2.0..fields.tyc2.1].parse()?;
        let tyc3: i64 = line[fields.tyc3.0..fields.tyc3.1].parse()?;
        let id_tycho = TycId(tyc1 * 1000000 + tyc2 * 10 + tyc3);
        let cmp = if line.len() < fields.cmp.1 {
            ArrayVec::new()
        } else {
            line[fields.cmp.0..fields.cmp.1]
                .as_bytes()
                .iter()
                .copied()
                .collect()
        };
        hip2tyc.add(hip, id_tycho, cmp);
    }

    Ok(())
}

fn read_tyc2_suppl1_readme(path: &Path) -> Result<Suppl1Fields, AppError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    loop {
        let line = lines
            .next()
            .ok_or(AppError::Parse("No information for suppl1.dat in file"))??;
        if line.starts_with("Byte-by-byte Description of file:") && line.contains("suppl_1.dat") {
            break;
        }
    }

    let mut tyc1 = None;
    let mut tyc2 = None;
    let mut tyc3 = None;
    let mut hip = None;
    let mut cmp = None;

    let mut separator_count = 0;
    while separator_count < 3 {
        let line = lines
            .next()
            .ok_or(AppError::Parse("Unexpected end of file"))??;
        if line.starts_with("----") {
            separator_count += 1;
            continue;
        }
        if separator_count < 2 {
            continue;
        }
        let range_str = line[..8].trim();
        if range_str.is_empty() {
            continue;
        }
        let range = match range_str.split_once('-') {
            Some((a, b)) => (a.trim().parse::<usize>()? - 1, b.trim().parse()?),
            None => {
                let start = range_str.trim().parse()?;
                (start - 1, start)
            }
        };

        let name = line[8..]
            .trim()
            .split_ascii_whitespace()
            .nth(2)
            .ok_or(AppError::Parse("Missing label"))?;
        match name {
            "TYC1" => tyc1 = Some(range),
            "TYC2" => tyc2 = Some(range),
            "TYC3" => tyc3 = Some(range),
            "HIP" => hip = Some(range),
            "CCDM" => cmp = Some(range),
            _ => (),
        }
    }

    Ok(Suppl1Fields {
        tyc1: tyc1.ok_or(AppError::MissingField(Cow::Borrowed(b"TYC1")))?,
        tyc2: tyc2.ok_or(AppError::MissingField(Cow::Borrowed(b"TYC2")))?,
        tyc3: tyc3.ok_or(AppError::MissingField(Cow::Borrowed(b"TYC3")))?,
        hip: hip.ok_or(AppError::MissingField(Cow::Borrowed(b"HIP")))?,
        cmp: cmp.ok_or(AppError::MissingField(Cow::Borrowed(b"CCDM")))?,
    })
}
