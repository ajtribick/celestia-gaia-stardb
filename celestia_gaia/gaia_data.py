# gaia-stardb: Processing Gaia DR2 for celestia.Sci/Celestia
# Copyright (C) 2019–2021  Andrew Tribick
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""Routines for downloading and cross-referencing Gaia data."""

from pathlib import Path
import re
import time

from astropy.table import Table
from astroquery.gaia import Gaia
import numpy as np

from .directories import GAIA_EDR3_DIR
from .ranges import MultiRange
from .utils import confirm_action
from .celestia_gaia import build_hip_xmatch, build_tyc_xmatch, get_source_ids


_HIP_MAX = 120404
_TYC_MAX = 9537


# --- DOWNLOAD CROSSMATCH CANDIDATES ---

def _hip_query(start: int, end: int) -> str:
    return f"""SELECT
	h.hip, h.plx AS hip_parallax, h.e_plx AS hip_parallax_error, h.ra AS hip_ra, h.dec AS hip_dec,
    h.hp_mag,
    g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.dr2_radial_velocity, g.ref_epoch,
    g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
    g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
    g.astrometric_params_solved
FROM
	public.hipparcos_newreduction h
	JOIN gaiaedr3.gaia_source g ON 1=CONTAINS(
		POINT(
			'ICRS',
			h.ra+COALESCE(h.pm_ra, 0)/COS(RADIANS(h.dec))*(2016.0-1991.25)/3600000.0,
			h.dec+COALESCE(h.pm_de, 0)*(2016.0-1991.25)/3600000.0
		),
		CIRCLE('ICRS', g.ra, g.dec, 2.0/60.0)
	)
WHERE
    h.hip BETWEEN {start} and {end}
	AND ABS(h.hp_mag-g.phot_g_mean_mag) <= 3
    """


def _tyc_query(start: int, end: int) -> str:
    id_start = start * 1000000
    id_end = (end+1) * 1000000 - 1
    return f"""SELECT
	t.id_tycho, t.hip, t.cmp, t.ra_deg AS tyc_ra, t.de_deg AS tyc_dec,
    t.bt_mag, t.vt_mag, t.ep_ra1990, t.ep_de1990,
    g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.dr2_radial_velocity, g.ref_epoch,
    g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
    g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
    g.astrometric_params_solved
FROM
	gaiaedr3.tycho2tdsc_merge t
	JOIN gaiaedr3.gaia_source g ON 1=CONTAINS(
		POINT(
			'ICRS',
			t.ra_deg+COALESCE(t.pm_ra, 0)/COS(RADIANS(t.de_deg))*(2016.0-t.ep_ra1990-1990.0)/3600000.0,
			t.de_deg+COALESCE(t.pm_de, 0)*(2016.0-t.ep_de1990-1990.0)/3600000.0
		),
		CIRCLE('ICRS', g.ra, g.dec, 2.0/60.0)
	)
WHERE
    t.id_tycho BETWEEN {id_start} AND {id_end}
	AND ABS(COALESCE(t.vt_mag, t.bt_mag)-g.phot_g_mean_mag) <= 3
    """


def _run_query(
    query: str, output_file: Path, upload_table: Table = None, upload_name: str = None
) -> None:
    job = Gaia.launch_job_async(
        query,
        dump_to_file=True,
        output_file=output_file,
        output_format='votable',
        upload_resource=upload_table,
        upload_table_name=upload_name,
        verbose=False,
        background=True,
    )

    print(f'  Launched job id {job.jobid}')
    delay=15
    while True:
        phase = job.get_phase(update=True)
        if job.is_finished():
            break
        print(f'  {phase}, waiting {delay} seconds')
        time.sleep(delay)
        delay = min(delay*2, 120)

    print(f'  {phase}')
    if phase != 'COMPLETED':
        print(f'  {phase}: {job.get_error()}')
        raise RuntimeError('Failed to download Gaia data')

    job.save_results()
    Gaia.remove_jobs([job.jobid])


def download_gaia_hip(ranges: MultiRange, chunk_size: int = 50000) -> None:
    """Download HIP data from the Gaia archive."""
    for section in ranges.chunk_ranges(chunk_size):
        hip_file = GAIA_EDR3_DIR/f'gaiaedr3-hip2-{section.begin:06}-{section.end:06}.vot.gz'

        query = _hip_query(section.begin, section.end)
        print(f'Querying HIP stars in range {section.begin} to {section.end}')
        _run_query(query, hip_file)


def download_gaia_tyc(ranges: MultiRange, chunk_size: int = 200) -> None:
    """Download TYC/TDSC data from the Gaia archive."""
    for section in ranges.chunk_ranges(chunk_size):
        tyc_file = (
            GAIA_EDR3_DIR/f'gaiaedr3-tyctdsc-part{section.begin:04}-{section.end:04}.vot.gz'
        )

        query = _tyc_query(section.begin, section.end)
        print(f'Querying TYC/TDSC stars in regions {section.begin} to {section.end}')
        _run_query(query, tyc_file)


_RANGE_PATTERN = re.compile(r'-([0-9]+)-([0-9]+)$')


def _getranges(start: int, end: int, path: Path, pattern: str) -> MultiRange:
    required_ranges = MultiRange(start, end)
    for existing in path.glob(pattern):
        match = _RANGE_PATTERN.search(existing.stem)
        if match:
            groups = match.groups()
            required_ranges.remove(int(groups[0]), int(groups[1]))
    return required_ranges


def download_gaia() -> None:
    """Download data from the Gaia archive."""
    GAIA_EDR3_DIR.mkdir(parents=True, exist_ok=True)

    hip_ranges = _getranges(1, _HIP_MAX, GAIA_EDR3_DIR, 'gaiaedr3-hip2-*.vot.gz')
    if not hip_ranges:
        if confirm_action('Hipparcos cross-match data already downloaded, replace?'):
            hip_ranges = MultiRange(1, _HIP_MAX)
    download_gaia_hip(hip_ranges)

    tyc_ranges = _getranges(1, _TYC_MAX, GAIA_EDR3_DIR, 'gaiaedr3-tyctdsc-*.vot.gz')
    if not tyc_ranges:
        if confirm_action('Tycho cross-match data already downloaded, replace?'):
            tyc_ranges = MultiRange(1, _TYC_MAX)
    download_gaia_tyc(tyc_ranges)


HIP_XMATCH = 'xmatch-gaia-hip.vot.gz'
TYC_XMATCH = 'xmatch-gaia-tyc.vot.gz'


def build_xmatches() -> None:
    """Build the cross-matches"""
    if (
        not (GAIA_EDR3_DIR / HIP_XMATCH).exists()
        or confirm_action('Re-generate Hipparcos cross-match?')
    ):
        build_hip_xmatch(GAIA_EDR3_DIR, HIP_XMATCH)

    if (
        not (GAIA_EDR3_DIR / TYC_XMATCH).exists()
        or confirm_action('Re-generate Tycho cross-match?')
    ):
        build_tyc_xmatch(GAIA_EDR3_DIR, TYC_XMATCH)


def download_gaia_distances(chunk_size: int = 250000) -> None:
    """Downloads the distances from the Gaia archive"""
    query = """SELECT
    d.source_id, d.r_med_geo, d.r_med_photogeo
FROM
    tap_upload.source_ids s
    JOIN external.gaiaedr3_distance d ON s.source_id = d.source_id
    """

    source_ids: np.ndarray = get_source_ids(GAIA_EDR3_DIR)
    source_ids = source_ids.astype('int64')  # https://github.com/numpy/numpy/issues/12264
    position = 0
    part = 1
    num_parts = (len(source_ids)+chunk_size-1) // chunk_size
    while position < len(source_ids):
        next_position = min(position + chunk_size, len(source_ids))
        print(f'Querying distances, part {part} of {num_parts}')

        table = Table([source_ids[position:next_position]], names=("source_id",))
        section_path = GAIA_EDR3_DIR/f'gaiaedr3-distance-{part:04}.vot.gz'
        _run_query(query, section_path, table, 'source_ids')

        position = next_position
        part += 1
