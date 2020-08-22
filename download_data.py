#!/usr/bin/env python3

# gaia-stardb: Processing Gaia DR2 for celestia.Sci/Celestia
# Copyright (C) 2019  Andrew Tribick
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

import contextlib
import getpass
import gzip
import os
import tarfile

from zipfile import ZipFile

import requests
import astropy.io.ascii as io_ascii

from astropy import units
from astropy.table import join, unique
from astroquery.gaia import Gaia
from astroquery.xmatch import XMatch

def yesno(prompt, default=False):
    """Prompt the user for yes/no input."""
    if default:
        new_prompt = prompt + ' (Y/n): '
    else:
        new_prompt = prompt + ' (y/N): '

    while True:
        answer = input(new_prompt)
        if answer == '':
            return default
        elif answer == 'y' or answer == 'Y':
            return True
        elif answer == 'n' or answer == 'N':
            return False

def proceed_checkfile(filename):
    """Check if a file exists, if so prompt the user if they want to replace it."""
    if os.path.exists(filename):
        if yesno(filename + ' already exists, replace?'):
            with contextlib.suppress(FileNotFoundError):
                os.remove(filename)
        else:
            return False
    return True

def download_file(outfile_name, url):
    """Download a file using requests."""
    if not proceed_checkfile(outfile_name):
        return

    print('Downloading ' + url)
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(outfile_name, 'wb') as f:
            f.write(response.raw.read())
    else:
        print('Failed to download')

# --- GAIA DATA DOWNLOAD ---

def download_gaia_data(colname, xindex_table, outfile_name):
    """Query and download Gaia data."""
    query = """SELECT
    x.source_id, x.original_ext_source_id AS """ + colname + """,
    g.ra, g.dec, g.parallax, g.parallax_error, g.pmra,
    g.pmdec, g.phot_g_mean_mag, g.bp_rp, g.teff_val,
    d.r_est, d.r_lo, d.r_hi
FROM
    """ + xindex_table + """ x
    JOIN gaiadr2.gaia_source g ON g.source_id = x.source_id
    LEFT JOIN external.gaiadr2_geometric_distance d ON d.source_id = x.source_id"""

    print(query)
    job = Gaia.launch_job_async(query,
                                dump_to_file=True,
                                output_file=outfile_name,
                                output_format='csv')
    try:
        job.save_results()
    finally:
        Gaia.remove_jobs(job.jobid)

CONESEARCH_URL = \
    'https://www.cosmos.esa.int/documents/29201/1769576/Hipparcos2GaiaDR2coneSearch.zip'

def download_gaia_hip(username):
    """Download HIP data from the Gaia archive."""
    hip_file = os.path.join('gaia', 'gaiadr2_hip-result.csv')
    if not proceed_checkfile(hip_file):
        return

    conesearch_file = os.path.join('gaia', 'hip2conesearch.zip')
    if proceed_checkfile(conesearch_file):
        download_file(conesearch_file, CONESEARCH_URL)

    # the gaiadr2.hipparcos2_best_neighbour table misses a large number of HIP stars that are
    # actually present, so use the mapping from Kervella et al. (2019) "Binarity of Hipparcos
    # stars from Gaia pm anomaly" instead.
    with tarfile.open(os.path.join('vizier', 'hipgpma.tar.gz'), 'r:gz') as tf:
        with tf.extractfile('./ReadMe') as readme:
            reader = io_ascii.get_reader(io_ascii.Cds,
                                            readme=readme,
                                            include_names=['HIP', 'GDR2'])
            reader.data.table_name = 'hipgpma.dat'
            with tf.extractfile('./hipgpma.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                hip_map = reader.read(f)

    hip_map = unique(hip_map)

    with ZipFile(conesearch_file, 'r') as csz:
        with csz.open('Hipparcos2GaiaDR2coneSearch.csv', 'r') as f:
            cone_map = io_ascii.read(f,
                                     format='csv',
                                     names=['HIP', 'GDR2', 'dist'],
                                     include_names=['HIP', 'GDR2'])

    cone_map = unique(cone_map)

    hip_map = join(hip_map, cone_map, join_type='outer', keys='HIP', table_names=['pm', 'cone'])
    hip_map['GDR2'] = hip_map['GDR2_pm'].filled(hip_map['GDR2_cone'])
    hip_map.remove_columns(['GDR2_pm', 'GDR2_cone'])
    hip_map.rename_column('HIP', 'original_ext_source_id')
    hip_map.rename_column('GDR2', 'source_id')

    Gaia.upload_table(upload_resource=hip_map, table_name='hipgpma')
    try:
        download_gaia_data('hip_id', 'user_'+username+'.hipgpma', hip_file)
    finally:
        Gaia.delete_user_table('hipgpma')

def download_gaia():
    """Download data from the Gaia archive."""
    with contextlib.suppress(FileExistsError):
        os.mkdir('gaia')

    print('Login to Gaia Archive')
    username = input('Username: ')
    if not username:
        print('Login aborted')
        return
    password = getpass.getpass('Password: ')
    if not password:
        print('Login aborted')
        return

    Gaia.login(user=username, password=password)
    try:
        download_gaia_hip(username)

        # download TYC data

        tyc_file = os.path.join('gaia', 'gaiadr2_tyc-result.csv')
        if proceed_checkfile(tyc_file):
            download_gaia_data('tyc2_id', 'gaiadr2.tycho2_best_neighbour', tyc_file)

    finally:
        Gaia.logout()

# --- SAO XMATCH DOWNLOAD ---

def download_xmatch(cat1, cat2, outfile_name):
    """Download a cross-match from VizieR."""
    if not proceed_checkfile(outfile_name):
        return

    result = XMatch.query(cat1=cat1,
                          cat2=cat2,
                          max_distance=5 * units.arcsec)

    io_ascii.write(result, outfile_name, format='csv')

def download_sao_xmatch():
    """Download cross-matches to the SAO catalogue."""
    with contextlib.suppress(FileExistsError):
        os.mkdir('xmatch')

    cross_matches = [
        ('vizier:I/131A/sao', 'vizier:I/311/hip2', 'sao_hip_xmatch.csv'),
        ('vizier:I/131A/sao', 'vizier:I/259/tyc2', 'sao_tyc2_xmatch.csv')
    ]

    for cat1, cat2, filename in cross_matches:
        print('Downloading '+cat1+'-'+cat2+' crossmatch')
        download_xmatch(cat1, cat2, os.path.join('xmatch', filename))

# --- VIZIER DOWNLOAD ---
def download_vizier():
    """Download catalogue archive files from VizieR."""
    with contextlib.suppress(FileExistsError):
        os.mkdir('vizier')

    files_urls = [
        ('ascc.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?I/280B'),
        ('tyc2spec.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?III/231'),
        ('tyc2teff.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?V/136'),
        ('ubvriteff.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J/ApJS/193/1'),
        ('xhip.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?V/137D'),
        ('hipgpma.tar.gz', 'https://cdsarc.unistra.fr/viz-bin/nph-Cat/tar.gz?J/A+A/623/A72'),
    ]

    for file_name, url in files_urls:
        download_file(os.path.join('vizier', file_name), url)

if __name__ == "__main__":
    download_vizier()
    download_gaia()
    download_sao_xmatch()
