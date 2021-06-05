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

import time
from pathlib import Path

import astropy.io.ascii as io_ascii
import requests
from astropy import units
from astroquery.gaia import Gaia
from astroquery.xmatch import XMatch

from .directories import GAIA_DR2_DIR, GAIA_EDR3_DIR, VIZIER_DIR, XMATCH_DIR


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
    hip, hip2_plx, hip2_e_plx, hip2_ra, hip2_dec, hp_mag,
    source_id, ra, dec, parallax, parallax_error,
    pmra, pmra_error, pmdec, pmdec_error,
    IF_THEN_ELSE(
        bp_rp > -20,
        TO_REAL(CASE_CONDITION(
            phot_g_mean_mag - 2.5*LOG10(
                1.00525 - 0.02323*GREATEST(0.25, LEAST(bp_rp, 3))
                + 0.01740*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 2)
                - 0.00253*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 3)
            ),
            astrometric_params_solved = 31,
            phot_g_mean_mag,
            phot_g_mean_mag < 13,
            phot_g_mean_mag,
            phot_g_mean_mag < 16,
            phot_g_mean_mag - 2.5*LOG10(
                1.00876 - 0.02540*GREATEST(0.25, LEAST(bp_rp, 3))
                + 0.01747*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 2)
                - 0.00277*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 3)
            )
        )),
        phot_g_mean_mag
    ) AS phot_g_mean_mag,
    bp_rp, ref_epoch,
    r_med_geo, r_med_photogeo,
    delta_mag,
    DISTANCE(hip2_pos, gaia_prop_pos)*3600 AS dist,
    CASE_CONDITION(
        delta_ra,
        delta_ra < -180, delta_ra + 360,
        delta_ra > 180, delta_ra - 360
    )*3600000/(ref_epoch-1991.25) AS calc_pmra,
    delta_dec*3600000/(ref_epoch-1991.25) AS calc_pmdec
FROM
    (
        SELECT
            hip2.hip, hip2.plx AS hip2_plx, hip2.e_plx AS hip2_e_plx,
            hip2.ra AS hip2_ra, hip2.dec AS hip2_dec, hip2.hp_mag,
            gaia.source_id, gaia.ra, gaia.dec, gaia.parallax, gaia.parallax_error,
            gaia.pmra, gaia.pmra_error, gaia.pmdec, gaia.pmdec_error,
            gaia.phot_g_mean_mag, gaia.bp_rp,
            gaia.ref_epoch, gaia.astrometric_params_solved,
            dist.r_med_geo, dist.r_med_photogeo,
            hip2.hp_mag-(
                0.91*COALESCE(gaia.phot_bp_mean_mag, gaia.phot_g_mean_mag)
                +0.09*COALESCE(gaia.phot_rp_mean_mag, gaia.phot_g_mean_mag)
            ) AS delta_mag,
            POINT('ICRS', hip2.ra, hip2.dec) AS hip2_pos,
            EPOCH_PROP_POS(
                gaia.ra, gaia.dec, gaia.parallax, gaia.pmra, gaia.pmdec,
                COALESCE(gaia.dr2_radial_velocity, 0),
                gaia.ref_epoch, 1991.25
            ) AS gaia_prop_pos,
            gaia.ra - hip2.ra AS delta_ra,
            gaia.dec - hip2.dec AS delta_dec
        FROM
            public.hipparcos_newreduction hip2
            JOIN gaiaedr3.gaia_source gaia ON 1=CONTAINS(
                POINT('ICRS', hip2.ra, hip2.dec),
                CIRCLE('ICRS', gaia.ra, gaia.dec, 0.05)
            )
            LEFT JOIN external.gaiaedr3_distance dist ON dist.source_id = gaia.source_id
        WHERE
            hip2.hip BETWEEN {start} AND {end}
    ) m
WHERE
    phot_g_mean_mag <= 13.7
    OR delta_mag BETWEEN -1 AND 0.5
    OR DISTANCE(hip2_pos, gaia_prop_pos) < 0.00027777777777777778
    """


def _tyc_query(start: int, end: int) -> str:
    id_start = start * 1000000
    id_end = (end+1) * 1000000 - 1
    return f"""SELECT
    tyc1, tyc2, tyc3, hip, cmp,
    tyc2_ra, tyc2_dec, tyc2_epoch, vt_mag, bt_mag,
    source_id, ra, dec, parallax, parallax_error,
    pmra, pmra_error, pmdec, pmdec_error,
    IF_THEN_ELSE(
        bp_rp > -20,
        TO_REAL(CASE_CONDITION(
            phot_g_mean_mag - 2.5*LOG10(
                1.00525 - 0.02323*GREATEST(0.25, LEAST(bp_rp, 3))
                + 0.01740*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 2)
                - 0.00253*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 3)
            ),
            astrometric_params_solved = 31,
            phot_g_mean_mag,
            phot_g_mean_mag < 13,
            phot_g_mean_mag,
            phot_g_mean_mag < 16,
            phot_g_mean_mag - 2.5*LOG10(
                1.00876 - 0.02540*GREATEST(0.25, LEAST(bp_rp, 3))
                + 0.01747*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 2)
                - 0.00277*POWER(GREATEST(0.25, LEAST(bp_rp, 3)), 3)
            )
        )),
        phot_g_mean_mag
    ) AS phot_g_mean_mag,
    bp_rp, ref_epoch,
    r_med_geo, r_med_photogeo,
    delta_mag,
    DISTANCE(tyc2_pos, gaia_prop_pos)*3600 AS dist,
    CASE_CONDITION(
        delta_ra,
        delta_ra < -180, delta_ra + 360,
        delta_ra > 180, delta_ra - 360
    )*3600000/(ref_epoch-tyc2_epoch) AS calc_pmra,
    delta_dec*3600000/(ref_epoch-tyc2_epoch) AS calc_pmdec
FROM
    (
        SELECT
            TO_INTEGER(tyc2.tyc1) AS tyc1,
            TO_INTEGER(tyc2.tyc2) AS tyc2,
            TO_INTEGER(tyc2.tyc3) AS tyc3,
            tyc2.ra_deg AS tyc2_ra, tyc2.de_deg AS tyc2_dec,
            0.5*(tyc2.ep_ra1990 + tyc2.ep_de1990)+1990 AS tyc2_epoch,
            tyc2.vt_mag, tyc2.bt_mag, tyc2.hip, tyc2.cmp,
            gaia.source_id, gaia.ra, gaia.dec, gaia.parallax, gaia.parallax_error,
            gaia.pmra, gaia.pmra_error, gaia.pmdec, gaia.pmdec_error,
            gaia.phot_g_mean_mag, gaia.bp_rp,
            gaia.ref_epoch, gaia.astrometric_params_solved,
            dist.r_med_geo, dist.r_med_photogeo,
            (
                0.99*COALESCE(tyc2.vt_mag, tyc2.bt_mag)
                +0.01*COALESCE(tyc2.bt_mag, tyc2.vt_mag)
            )-(
                0.91*COALESCE(gaia.phot_bp_mean_mag, gaia.phot_g_mean_mag)
                +0.09*COALESCE(gaia.phot_rp_mean_mag, gaia.phot_g_mean_mag)
            ) AS delta_mag,
            POINT('ICRS', tyc2.ra, tyc2.dec) AS tyc2_pos,
            EPOCH_PROP_POS(
                gaia.ra, gaia.dec, gaia.parallax, gaia.pmra, gaia.pmdec,
                COALESCE(gaia.dr2_radial_velocity, 0),
                gaia.ref_epoch, 0.5*(tyc2.ep_ra1990 + tyc2.ep_de1990)+1990
            ) AS gaia_prop_pos,
            gaia.ra - tyc2.ra AS delta_ra,
            gaia.dec - tyc2.dec AS delta_dec
        FROM
            gaiaedr3.tycho2tdsc_merge tyc2
            JOIN gaiaedr3.gaia_source gaia ON 1=CONTAINS(
                POINT('ICRS', tyc2.ra, tyc2.dec),
                CIRCLE('ICRS', gaia.ra, gaia.dec, 0.05)
            )
            LEFT JOIN external.gaiaedr3_distance dist ON dist.source_id = gaia.source_id
        WHERE
            tyc2.id_tycho BETWEEN {id_start} AND {id_end}
    ) m
WHERE
    delta_mag BETWEEN -1 AND 0.5
    OR DISTANCE(tyc2_pos, gaia_prop_pos) < 0.00027777777777777778
    """

def download_gaia_hip(
    chunk_size: int = 5000, *, begin_section: int = 0, begin_at: int = 1
) -> None:
    """Download HIP data from the Gaia archive."""

    section = begin_section
    start = begin_at
    end = begin_at+chunk_size-1
    while start <= _HIP_MAX:
        hip_file = GAIA_EDR3_DIR/f'gaiaedr3-hip2-part{section:02}.votable'

        query = _hip_query(start, end)
        print(f'Querying HIP stars in range {start} to {end}')
        job = Gaia.launch_job_async(
            query,
            dump_to_file=True,
            output_file=hip_file,
            output_format='votable',
            verbose=False,
            background=True,
        )

        print(f'  Launched job id {job.jobid}')
        delay=10
        while True:
            phase = job.get_phase(update=True)
            if job.is_finished():
                break
            print(f'  {phase}, waiting {delay} seconds')
            time.sleep(delay)
            delay = min(delay+10, 60)

        print(f'  {phase}')
        if phase != 'COMPLETED':
            raise RuntimeError('Failed to download Gaia data')

        job.save_results()
        Gaia.remove_jobs([job.jobid])

        section += 1
        start = end+1
        end = min(end+chunk_size, _HIP_MAX)


def download_gaia_tyc(
    chunk_size = 20, *, begin_section: int = 0, begin_at = 1
) -> None:
    """Download TYC/TDSC data from the Gaia archive."""

    section = begin_section
    start = begin_at
    end = begin_at+chunk_size-1
    while start <= _TYC_MAX:
        hip_file = GAIA_EDR3_DIR/f'gaiaedr3-tyctdsc-part{section:02}.votable'

        query = _tyc_query(start, end)
        print(f'Querying TYC/TDSC stars in regions {start} to {end}')
        job = Gaia.launch_job_async(
            query,
            dump_to_file=True,
            output_file=hip_file,
            output_format='votable',
            verbose=False,
            background=True,
        )

        print(f'  Launched job id {job.jobid}')
        delay=10
        while True:
            phase = job.get_phase(update=True)
            if job.is_finished():
                break
            print(f'  {phase}, waiting {delay} seconds')
            time.sleep(delay)
            delay = min(delay+10, 60)

        print(f'  {phase}')
        if phase != 'COMPLETED':
            raise RuntimeError('Failed to download Gaia data')

        job.save_results()
        Gaia.remove_jobs([job.jobid])

        section += 1
        start = end+1
        end = min(end+chunk_size, _HIP_MAX)


def download_gaia() -> None:
    """Download data from the Gaia archive."""
    GAIA_DR2_DIR.mkdir(parents=True, exist_ok=True)
    GAIA_EDR3_DIR.mkdir(parents=True, exist_ok=True)

    download_gaia_hip()
    download_gaia_tyc()


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
