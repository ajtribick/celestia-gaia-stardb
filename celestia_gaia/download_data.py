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

from pathlib import Path, PurePath
from typing import Union, cast
from zipfile import ZipFile

import astropy.io.ascii as io_ascii
import astropy.io.votable as votable
import numpy as np
import requests
from astropy import units
from astropy.table import Table, join, unique, vstack
from astroquery.gaia import Gaia
from astroquery.utils.tap import Tap
from astroquery.xmatch import XMatch

from .parse_utils import open_cds_tarfile

GAIA_DIR = Path('data', 'gaiadr2')
SIMBAD_DIR = Path('data', 'simbad')
VIZIER_DIR = Path('data', 'vizier')
XMATCH_DIR = Path('data', 'xmatch')


def yesno(prompt: str, default: bool=False) -> bool:
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


def proceed_checkfile(path: Path) -> bool:
    """Check if a file exists, if so prompt the user if they want to replace it."""
    if path.exists():
        if yesno(f'{path} already exists, replace?'):
            path.unlink()
        else:
            return False
    return True


def download_file(path: Path, url: str) -> bool:
    """Download a file using requests."""
    if not proceed_checkfile(path):
        return

    print(f'Downloading {url}')
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with path.open('wb') as f:
            f.write(response.raw.read())
    else:
        print('Failed to download')


# --- GAIA DATA DOWNLOAD ---

def download_gaia_data(colname: str, xindex: Union[str, Table], outfile_path: Path) -> None:
    """Query and download Gaia data."""

    if isinstance(xindex, str):
        xindex_name = cast(str, xindex)
        upload_resource = None
        upload_table_name = None
    else:
        xindex_name = "TAP_UPLOAD.cel_xindex"
        upload_resource = cast(Table, xindex)
        upload_table_name = 'cel_xindex'

    query = f"""SELECT
    x.source_id, x.original_ext_source_id AS {colname},
    g.ra, g.dec, g.parallax, g.parallax_error, g.pmra,
    g.pmdec, g.phot_g_mean_mag, g.bp_rp, g.teff_val,
    d.r_est, d.r_lo, d.r_hi
FROM
    {xindex_name} x
    JOIN gaiadr2.gaia_source g ON g.source_id = x.source_id
    LEFT JOIN external.gaiadr2_geometric_distance d ON d.source_id = x.source_id"""

    print(query)
    job = Gaia.launch_job_async(
        query,
        upload_resource=upload_resource,
        upload_table_name=upload_table_name,
        dump_to_file=True,
        output_file=outfile_path,
        output_format='csv',
    )
    try:
        job.save_results()
    finally:
        Gaia.remove_jobs(job.jobid)


CONESEARCH_URL = (
    'https://www.cosmos.esa.int/documents/29201/1769576/Hipparcos2GaiaDR2coneSearch.zip'
)


def download_gaia_hip() -> None:
    """Download HIP data from the Gaia archive."""
    hip_file = GAIA_DIR/'gaiadr2_hip-result.csv'
    if not proceed_checkfile(hip_file):
        return

    conesearch_file = GAIA_DIR/'hip2conesearch.zip'
    if proceed_checkfile(conesearch_file):
        download_file(conesearch_file, CONESEARCH_URL)

    # the gaiadr2.hipparcos2_best_neighbour table misses a large number of HIP stars that are
    # actually present, so use the mapping from Kervella et al. (2019) "Binarity of Hipparcos
    # stars from Gaia pm anomaly" instead.

    with open_cds_tarfile(VIZIER_DIR/'hipgpma.tar.gz') as tf:
        hip_map = unique(tf.read_gzip('hipgpma.dat', ['HIP', 'GDR2']))

    with ZipFile(conesearch_file, 'r') as csz:
        with csz.open('Hipparcos2GaiaDR2coneSearch.csv', 'r') as f:
            cone_map = io_ascii.read(
                f, format='csv', names=['HIP', 'GDR2', 'dist'],
                include_names=['HIP', 'GDR2'],
            )

    cone_map = unique(cone_map)

    hip_map = join(hip_map, cone_map, join_type='outer', keys='HIP', table_names=['pm', 'cone'])
    hip_map['GDR2'] = hip_map['GDR2_pm'].filled(hip_map['GDR2_cone'])
    hip_map.remove_columns(['GDR2_pm', 'GDR2_cone'])
    hip_map.rename_column('HIP', 'original_ext_source_id')
    hip_map.rename_column('GDR2', 'source_id')

    download_gaia_data('hip_id', hip_map, hip_file)


def _load_gaia_tyc_ids(path: Path) -> Table:
    with path.open('r') as f:
        header = f.readline().split(',')
        col_idx = header.index('tyc2_id')
        tyc1 = []
        tyc2 = []
        tyc3 = []
        for line in f:
            try:
                tyc2_id = line.split(',')[col_idx]
            except IndexError:
                continue

            tyc = tyc2_id.split('-')
            tyc1.append(int(tyc[0]))
            tyc2.append(int(tyc[1]))
            tyc3.append(int(tyc[2]))

        return Table([tyc1, tyc2, tyc3], names=['TYC1','TYC2','TYC3'], dtype=('i4', 'i4', 'i4'))


def _load_ascc_tyc_ids(filename: Path) -> Table:
    data = None
    with open_cds_tarfile(filename) as tf:
        for data_file in tf.tf:
            path = PurePath(data_file.name)
            if path.parent != PurePath('.') or not path.stem.startswith('cc'):
                continue

            section_data = tf.read_gzip(
                path.stem, ['TYC1', 'TYC2', 'TYC3'], readme_name='cc*.dat',
            )

            if data is None:
                data = section_data
            else:
                data = vstack([data, section_data], join_type='exact')

    return data


def get_missing_tyc_ids(tyc_file: Path, ascc_file: Path) -> Table:
    """Finds the ASCC TYC ids that are not present in Gaia cross-match."""
    print("Finding missing TYC ids in ASCC")
    t_asc = unique(_load_ascc_tyc_ids(ascc_file))
    t_gai = _load_gaia_tyc_ids(tyc_file)

    t_gai['in_gaia'] = True

    t_mgd = join(t_asc, t_gai, join_type='left')
    t_mgd['in_gaia'] = t_mgd['in_gaia'].filled(False)

    t_missing = t_mgd[np.logical_not(t_mgd['in_gaia'])]
    t_missing = t_missing[t_missing['TYC1'] != 0] # remove invalid entries

    return Table([[f"TYC {t['TYC1']}-{t['TYC2']}-{t['TYC3']}" for t in t_missing]], names=['id'])


def download_gaia_tyc() -> None:
    """Download TYC data from the Gaia archive."""

    tyc_file = GAIA_DIR/'gaiadr2_tyc-result.csv'
    if proceed_checkfile(tyc_file):
        download_gaia_data('tyc2_id', 'gaiadr2.tycho2_best_neighbour', tyc_file)

    # Use SIMBAD to fill in some of the missing entries
    SIMBAD_DIR.mkdir(parents=True, exist_ok=True)

    simbad_file = SIMBAD_DIR/'tyc-gaia.votable'
    if proceed_checkfile(simbad_file):
        ascc_file = VIZIER_DIR/'ascc.tar.gz'
        missing_ids = get_missing_tyc_ids(tyc_file, ascc_file)
        print("Querying SIMBAD for Gaia DR2 identifiers")
        simbad = Tap(url='http://simbad.u-strasbg.fr:80/simbad/sim-tap')
        query = """SELECT
    id1.id tyc_id, id2.id gaia_id
FROM
    TAP_UPLOAD.missing_tyc src
    JOIN IDENT id1 ON id1.id = src.id
    JOIN IDENT id2 ON id2.oidref = id1.oidref
WHERE
    id2.id LIKE 'Gaia DR2 %'"""
        print(query)
        job = simbad.launch_job_async(
            query,
            upload_resource=missing_ids,
            upload_table_name='missing_tyc',
            output_file=simbad_file,
            output_format='votable',
            dump_to_file=True,
        )
        job.save_results()

    tyc2_file = GAIA_DIR/'gaiadr2_tyc-result-extra.csv'
    if proceed_checkfile(tyc2_file):
        missing_ids = votable.parse(simbad_file).resources[0].tables[0].to_table()
        missing_ids['tyc_id'] = [
            m[m.rfind(' ')+1:] for m in missing_ids['tyc_id'].astype('U')
        ]
        missing_ids.rename_column('tyc_id', 'original_ext_source_id')
        missing_ids['gaia_id'] = [
            int(m[m.rfind(' ')+1:]) for m in missing_ids['gaia_id'].astype('U')
        ]
        missing_ids.rename_column('gaia_id', 'source_id')

        download_gaia_data('tyc2_id', missing_ids, tyc2_file)


def download_gaia() -> None:
    """Download data from the Gaia archive."""
    GAIA_DIR.mkdir(parents=True, exist_ok=True)

    download_gaia_hip()
    download_gaia_tyc()


# --- SAO XMATCH DOWNLOAD ---

def download_xmatch(cat1: str, cat2: str, path: Path) -> None:
    """Download a cross-match from VizieR."""
    if not proceed_checkfile(path):
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
        download_file(VIZIER_DIR/file_name, url)
