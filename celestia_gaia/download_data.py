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

"""Routines for downloading the data files."""

import re
import time
from pathlib import Path

import astropy.io.ascii as io_ascii
import requests
from astropy import units
from astroquery.gaia import Gaia
from astroquery.xmatch import XMatch

from .directories import GAIA_EDR3_DIR, VIZIER_DIR, XMATCH_DIR
from .ranges import MultiRange


def _yesno(prompt: str, default: bool=False) -> bool:
    """Prompt the user for yes/no input."""
    if default:
        new_prompt = f'{prompt} (Y/n): '
    else:
        new_prompt = f'{prompt} (y/N): '

    while True:
        answer = input(new_prompt)
        if answer == '':
            return default
        if answer in ('y', 'Y'):
            return True
        if answer in ('n', 'N'):
            return False


def _proceed_checkfile(path: Path) -> bool:
    """Check if a file exists, if so prompt the user if they want to replace it."""
    if path.exists():
        if _yesno(f'{path} already exists, replace?'):
            path.unlink()
        else:
            return False
    return True


def _download_file(path: Path, url: str) -> bool:
    """Download a file using requests."""
    if not _proceed_checkfile(path):
        return

    print(f'Downloading {url}')
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with path.open('wb') as f:
            f.write(response.raw.read())
    else:
        print('Failed to download')


# --- GAIA DATA DOWNLOAD ---

_HIP_MAX = 120404
_TYC_MAX = 9537


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
			h.ra+COALESCE(h.pm_ra, 0)*COS(RADIANS(h.dec))*(2016.0-1991.25)/3600000.0,
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
	t.tyc1, t.tyc2, t.tyc3, t.hip, t.cmp, t.ra_deg AS tyc_ra, t.de_deg AS tyc_dec,
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
			t.ra_deg+COALESCE(t.pm_ra, 0)*COS(RADIANS(t.de_deg))*(2016.0-t.ep_ra1990-1990.0)/3600000.0,
			t.de_deg+COALESCE(t.pm_de, 0)*(2016.0-t.ep_de1990-1990.0)/3600000.0
		),
		CIRCLE('ICRS', g.ra, g.dec, 2.0/60.0)
	)
WHERE
    t.id_tycho BETWEEN {id_start} AND {id_end}
	AND ABS(COALESCE(t.vt_mag, t.bt_mag)-g.phot_g_mean_mag) <= 3
    """


def _run_query(query: str, output_file: Path) -> None:
    job = Gaia.launch_job_async(
        query,
        dump_to_file=True,
        output_file=output_file,
        output_format='votable',
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
        hip_file = GAIA_EDR3_DIR/f'gaiaedr3-hip2-{section.begin:06}-{section.end:06}.votable'

        query = _hip_query(section.begin, section.end)
        print(f'Querying HIP stars in range {section.begin} to {section.end}')
        _run_query(query, hip_file)


def download_gaia_tyc(ranges: MultiRange, chunk_size: int = 200) -> None:
    """Download TYC/TDSC data from the Gaia archive."""
    for section in ranges.chunk_ranges(chunk_size):
        tyc_file = (
            GAIA_EDR3_DIR/f'gaiaedr3-tyctdsc-part{section.begin:04}-{section.end:04}.votable'
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

    hip_ranges = _getranges(1, _HIP_MAX, GAIA_EDR3_DIR, 'gaiaedr3-hip2-*.votable')
    if not hip_ranges:
        if _yesno('Hipparcos cross-match data already downloaded, replace?'):
            hip_ranges = MultiRange(1, _HIP_MAX)
    download_gaia_hip(hip_ranges)

    tyc_ranges = _getranges(1, _TYC_MAX, GAIA_EDR3_DIR, 'gaiaedr3-tyctdsc-*.votable')
    if not tyc_ranges:
        if _yesno('Tycho cross-match data already downloaded, replace?'):
            tyc_ranges = MultiRange(1, _TYC_MAX)
    download_gaia_tyc(tyc_ranges)


# --- SAO XMATCH DOWNLOAD ---

def download_xmatch(cat1: str, cat2: str, path: Path) -> None:
    """Download a cross-match from VizieR."""
    if not _proceed_checkfile(path):
        return

    result = XMatch.query(
        cat1=cat1, cat2=cat2, max_distance=5 * units.arcsec,
    )

    io_ascii.write(result, path, format='csv')


def download_sao_xmatch() -> None:
    """Download cross-matches to the SAO catalogue."""
    XMATCH_DIR.mkdir(parents=True, exist_ok=True)

    cross_matches = [
        ('vizier:I/131A/sao', 'vizier:I/311/hip2', 'sao_hip_xmatch.csv'),
        ('vizier:I/131A/sao', 'vizier:I/259/tyc2', 'sao_tyc2_xmatch.csv'),
        ('vizier:I/131A/sao', 'vizier:I/259/suppl_1', 'sao_tyc2_suppl1_xmatch.csv'),
        ('vizier:I/131A/sao', 'vizier:I/259/suppl_2', 'sao_tyc2_suppl2_xmatch.csv'),
    ]

    for cat1, cat2, filename in cross_matches:
        print(f'Downloading {cat1}-{cat2} crossmatch')
        download_xmatch(cat1, cat2, XMATCH_DIR/filename)


# --- VIZIER DOWNLOAD ---

def download_vizier() -> None:
    """Download catalogue archive files from VizieR."""
    VIZIER_DIR.mkdir(parents=True, exist_ok=True)

    files_urls = [
        ('ascc.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?I/280B'),
        ('hipgpma.tar.gz', 'https://cdsarc.unistra.fr/viz-bin/nph-Cat/tar.gz?J/A+A/623/A72'),
        # for some reason, the SAO archive at VizieR does not work, so download files individually
        ('sao.dat.gz', 'https://cdsarc.unistra.fr/ftp/I/131A/sao.dat.gz'),
        ('sao.readme', 'https://cdsarc.unistra.fr/ftp/I/131A/ReadMe'),
        ('tyc2hd.tar.gz', 'https://cdsarc.unistra.fr/viz-bin/nph-Cat/tar.gz?IV/25'),
        ('tyc2spec.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?III/231'),
        ('tyc2specnew.tar.gz', 'https://cdsarc.unistra.fr/viz-bin/nph-Cat/tar.gz?J/PAZh/34/21'),
        ('tyc2teff.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?V/136'),
        ('ubvriteff.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J/ApJS/193/1'),
        ('xhip.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?V/137D'),
    ]

    for file_name, url in files_urls:
        _download_file(VIZIER_DIR/file_name, url)
