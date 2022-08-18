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
import gzip
import re
import time

import astropy.io.ascii as io_ascii
from astropy.table import Table
from astroquery.gaia import Gaia
import numpy as np

from .directories import GAIA_EDR3_DIR, VIZIER_DIR
from .ranges import MultiRange
from .utils import confirm_action
from .celestia_gaia import build_xmatch, get_required_dist_source_ids


_HIP_MAX = 120404
_TYC_MAX = 9537
_MAX_MAG_DIFF = 4


# --- DOWNLOAD CROSSMATCH CANDIDATES ---

def _hip1_query(upload_name: str) -> str:
    return f"""SELECT
    h.hip, h.ra AS hip_ra, h.dec AS hip_dec, h.hp_mag,
    g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.dr2_radial_velocity, g.ref_epoch,
    g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
    g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
    g.astrometric_params_solved
FROM
    tap_upload.{upload_name} h
    JOIN gaiaedr3.gaia_source g ON 1=CONTAINS(
        POINT('ICRS', h.ra, h.dec),
        CIRCLE('ICRS', g.ra, g.dec, 2.0/60.0)
    )
WHERE ABS(h.hp_mag-g.phot_g_mean_mag) <= {_MAX_MAG_DIFF}
    """


def _hip2_query(start: int, end: int) -> str:
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
    AND ABS(h.hp_mag-g.phot_g_mean_mag) <= {_MAX_MAG_DIFF}
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
    AND ABS(COALESCE(t.vt_mag, t.bt_mag)-g.phot_g_mean_mag) <= {_MAX_MAG_DIFF}
    """


def _tyc2_supplement1_query(upload_name: str) -> str:
    return f"""SELECT
    t.tyc1, t.tyc2, t.tyc3, t.hip, t.cmp, t.ra_deg AS tyc_ra, t.de_deg AS tyc_dec,
    t.bt_mag, t.vt_mag,
    g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.dr2_radial_velocity, g.ref_epoch,
    g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
    g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
    g.astrometric_params_solved
FROM
    tap_upload.{upload_name} t
    JOIN gaiaedr3.gaia_source g ON 1=CONTAINS(
        POINT(
            'ICRS',
            t.ra_deg+COALESCE(t.pm_ra, 0)/COS(RADIANS(t.de_deg))*24.75/3600000.0,
            t.de_deg+COALESCE(t.pm_de, 0)*24.75/3600000.0
        ),
        CIRCLE('ICRS', g.ra, g.dec, 2.0/60.0)
    )
WHERE
    ABS(COALESCE(t.vt_mag, t.bt_mag)-g.phot_g_mean_mag) <= {_MAX_MAG_DIFF}
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


def download_gaia_hip2(ranges: MultiRange, chunk_size: int = 50000) -> None:
    """Download HIP2 data from the Gaia archive."""
    for section in ranges.chunk_ranges(chunk_size):
        hip_file = GAIA_EDR3_DIR/f'gaiaedr3-hip2-{section.begin:06}-{section.end:06}.vot.gz'

        query = _hip2_query(section.begin, section.end)
        print(f'Querying HIP stars in range {section.begin} to {section.end}')
        _run_query(query, hip_file)


def download_gaia_hip1() -> None:
    """Donwload HIP1 data from the Gaia archive."""
    hip1_file = GAIA_EDR3_DIR/'gaiaedr3-hip1.vot.gz'
    if (
        hip1_file.exists()
        and not confirm_action('Re-download HIP1 crossmatch?')
    ):
        return

    table = Table.read(GAIA_EDR3_DIR/'hip1_subset.vot.gz', format='votable')
    table.remove_columns(['hd', 'b_v', 'e_b_v', 'v_i', 'e_v_i', 'sptype'])

    # Parse sexagesimal RA and Dec coordinates
    table = table[np.logical_not(np.logical_or(table['rahms'].mask, table['dedms'].mask))]
    table['rahms'] = table['rahms'].filled('')
    table['dedms'] = table['dedms'].filled('')
    pcoord = np.vectorize(_parse_coord)
    table['ra'] = pcoord(table['rahms']) * 15.0
    table['dec'] = pcoord(table['dedms'])
    table.remove_columns(['rahms', 'dedms'])

    # Use Vmag as a proxy for Hpmag when Hpmag is missing
    table = table[np.logical_not(np.logical_and(table['hpmag'].mask, table['vmag'].mask))]
    table['hpmag'] = table['hpmag'].filled(table['vmag'].filled(np.nan))
    table.remove_column('vmag')
    table.rename_column('hpmag', 'hp_mag')

    query = _hip1_query('hip1_subset')
    _run_query(query, hip1_file, table, 'hip1_subset')


def download_gaia_tyctdsc(ranges: MultiRange, chunk_size: int = 200) -> None:
    """Download TYC/TDSC data from the Gaia archive."""
    for section in ranges.chunk_ranges(chunk_size):
        tyc_file = (
            GAIA_EDR3_DIR/f'gaiaedr3-tyctdsc-{section.begin:04}-{section.end:04}.vot.gz'
        )

        query = _tyc_query(section.begin, section.end)
        print(f'Querying TYC/TDSC stars in regions {section.begin} to {section.end}')
        _run_query(query, tyc_file)


def download_gaia_tyc2_supplement1() -> None:
    """Download TYC2 Supplement 1 data from the Gaia archive."""
    tyc2_supplement1_file = GAIA_EDR3_DIR/'gaiaedr3-tyc2suppl1.vot.gz'
    if (
        tyc2_supplement1_file.exists()
        and not confirm_action('Re-download TYC2 supplement 1 crossmatch?')
    ):
        return

    reader = io_ascii.get_reader(
        io_ascii.Cds,
        readme=str(VIZIER_DIR/'tyc2.readme'),
        include_names=[
            'TYC1', 'TYC2', 'TYC3', 'RAdeg', 'DEdeg', 'pmRA', 'pmDE', 'BTmag', 'VTmag',
            'HIP', 'CCDM',
        ],
    )
    reader.data.table_name = 'suppl_1.dat'
    with gzip.open(VIZIER_DIR/'tyc2suppl_1.dat.gz', 'rb') as gzf:
        table = reader.read(gzf)

    table.rename_columns(
        [
            'TYC1', 'TYC2', 'TYC3',
            'RAdeg', 'DEdeg', 'pmRA', 'pmDE',
            'BTmag', 'VTmag',
            'HIP', 'CCDM',
        ],
        [
            'tyc1', 'tyc2', 'tyc3',
            'ra_deg', 'de_deg', 'pm_ra', 'pm_de',
            'bt_mag', 'vt_mag',
            'hip', 'cmp',
        ],
    )

    query = _tyc2_supplement1_query('tyc2suppl1')
    _run_query(query, tyc2_supplement1_file, table, 'tyc2suppl1')


def download_tyc2tdsc_xmatch():
    """Download the TYC2-HIP cross-index from the Gaia archive."""
    tyc2tdsc_xmatch_file = GAIA_EDR3_DIR/'tyc2tdsc_hip_xmatch.vot.gz'
    if (
        tyc2tdsc_xmatch_file.exists()
        and not confirm_action('Re-download TYC2TDSC-HIP identifier map?')
    ):
        return

    query = 'SELECT id_tycho, hip, cmp FROM tycho2tdsc_merge WHERE hip IS NOT NULL'
    print('Querying for HIP-TYC crossmatches from tycho2tdsc_merge')
    _run_query(query, tyc2tdsc_xmatch_file)

_RANGE_PATTERN = re.compile(r'-([0-9]+)-([0-9]+)\.')


def _parse_coord(coord_str: str) -> float:
    coord_str = coord_str.strip()
    coord_parts = [float(p) for p in coord_str.split()]
    if coord_str.startswith('-'):
        coord = coord_parts[0] - coord_parts[1]/60 - coord_parts[2]/3600
    else:
        coord = coord_parts[0] + coord_parts[1]/60 + coord_parts[2]/3600
    return coord


def download_hip1_subset() -> None:
    """Downloads details of HIP1 stars not in HIP2."""
    query = """SELECT
    h.hip, h.hd, h.rahms, h.dedms, h.hpmag, h.vmag, h.b_v, h.e_b_v, h.v_i, h.e_v_i, h.sptype
FROM
    public.hipparcos h
    LEFT JOIN public.hipparcos_newreduction h2 ON h2.hip = h.hip
WHERE
    h2.hip IS NULL
    """
    hip1_subset_file = GAIA_EDR3_DIR/'hip1_subset.vot.gz'
    if (
        not hip1_subset_file.exists()
        or confirm_action('Re-download HIP1 subset?')
    ):
        _run_query(query, hip1_subset_file)


def _getranges(start: int, end: int, path: Path, pattern: str) -> MultiRange:
    required_ranges = MultiRange(start, end)
    for existing in path.glob(pattern):
        match = _RANGE_PATTERN.search(str(existing))
        if match:
            groups = match.groups()
            required_ranges.remove(int(groups[0]), int(groups[1]))
    return required_ranges


def download_gaia() -> None:
    """Download data from the Gaia archive."""
    GAIA_EDR3_DIR.mkdir(parents=True, exist_ok=True)

    download_hip1_subset()
    download_tyc2tdsc_xmatch()

    hip_ranges = _getranges(1, _HIP_MAX, GAIA_EDR3_DIR, 'gaiaedr3-hip2-*.vot.gz')
    if not hip_ranges:
        if confirm_action('HIP2-Gaia cross-match data already downloaded, replace?'):
            hip_ranges = MultiRange(1, _HIP_MAX)
    download_gaia_hip2(hip_ranges)

    download_gaia_hip1()

    tyc_ranges = _getranges(1, _TYC_MAX, GAIA_EDR3_DIR, 'gaiaedr3-tyctdsc-*.vot.gz')
    if not tyc_ranges:
        if confirm_action('TYC2TDSC-Gaia cross-match data already downloaded, replace?'):
            tyc_ranges = MultiRange(1, _TYC_MAX)
    download_gaia_tyctdsc(tyc_ranges)

    download_gaia_tyc2_supplement1()


def build_xmatches() -> None:
    """Build the cross-matches"""
    if (
        not (GAIA_EDR3_DIR/'xmatch-gaia-hiptyc.vot.gz').exists()
        or confirm_action('Re-generate HIP/TYC cross-match?')
    ):
        build_xmatch(GAIA_EDR3_DIR, VIZIER_DIR, 'xmatch-gaia-hiptyc.vot.gz')


def download_gaia_distances(chunk_size: int = 250000) -> None:
    """Downloads the distances from the Gaia archive"""
    query = """SELECT
    s.source_id, d.r_med_geo, d.r_med_photogeo
FROM
    tap_upload.source_ids s
    LEFT JOIN external.gaiaedr3_distance d ON s.source_id = d.source_id
    """

    source_ids: np.ndarray = get_required_dist_source_ids(GAIA_EDR3_DIR)
    if len(source_ids) == 0 and confirm_action('Re-download distances?'):
        for f in GAIA_EDR3_DIR.glob('gaiaedr3-distance-*.vot.gz'):
            f.unlink()
        source_ids = get_required_dist_source_ids(GAIA_EDR3_DIR)
    source_ids = source_ids.astype('int64')  # https://github.com/numpy/numpy/issues/12264
    position = 0
    part = 1
    for f in GAIA_EDR3_DIR.glob('gaiaedr3-distance-*.vot.gz'):
        fpart = int(f.name.removeprefix('gaiaedr3-distance-').removesuffix('.vot.gz'))
        if fpart >= part:
            part = fpart + 1

    p = 1
    num_parts = (len(source_ids)+chunk_size-1) // chunk_size
    while position < len(source_ids):
        next_position = min(position + chunk_size, len(source_ids))
        print(f'Querying distances, part {p} of {num_parts}')

        table = Table([source_ids[position:next_position]], names=("source_id",))
        section_path = GAIA_EDR3_DIR/f'gaiaedr3-distance-{part:04}.vot.gz'
        _run_query(query, section_path, table, 'source_ids')

        position = next_position
        part += 1
        p += 1
