#!/usr/bin/env python3

import contextlib, getpass, os
import astropy, requests
from zipfile import ZipFile
from astropy import units
from astropy.io import ascii
from astroquery.gaia import Gaia
from astroquery.xmatch import XMatch
from requests.exceptions import HTTPError

def yesno(prompt, default=False):
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
    if os.path.exists(filename):
        if yesno(filename + ' already exists, replace?'):
            with contextlib.suppress(FileNotFoundError):
                os.remove(filename)
        else:
            return False
    return True

def download_file(outfile_name, url):
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
    if not proceed_checkfile(outfile_name):
        return

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

def download_gaia():
    with contextlib.suppress(FileExistsError):
        os.mkdir("gaia")

    print("Login to Gaia Archive")
    username = input("Username: ")
    if not username:
        print("Login aborted")
        return
    password = getpass.getpass("Password: ")
    if not password:
        print("Login aborted")
        return

    Gaia.login(user=username, password=password)
    try:
        # the gaiadr2.hipparcos2_best_neighbour table misses a large number of HIP stars that are actually present
        # so use the cone search file
        conesearch_file = os.path.join('gaia', 'hip2conesearch.zip')
        download_file(
            conesearch_file,
            'https://www.cosmos.esa.int/documents/29201/1769576/Hipparcos2GaiaDR2coneSearch.zip')

        with ZipFile(conesearch_file, 'r') as csz:
            with csz.open('Hipparcos2GaiaDR2coneSearch.csv', 'r') as f:
                hip_map = ascii.read(f, names=['original_ext_source_id', 'source_id', 'dist'])

        Gaia.upload_table(upload_resource=hip_map, table_name='hip_cone')
        try:
            download_gaia_data('hip_id', 'user_' + username + '.hip_cone', os.path.join('gaia', 'gaiadr2_hip-result.csv'))
            download_gaia_data('tyc2_id', 'gaiadr2.tycho2_best_neighbour', os.path.join('gaia', 'gaiadr2_tyc-result.csv'))
        finally:
            Gaia.delete_user_table('hip_cone')

    finally:
        Gaia.logout()

# --- SAO XMATCH DOWNLOAD ---

def download_xmatch(cat1, cat2, outfile_name):
    if not proceed_checkfile(outfile_name):
        return

    result = XMatch.query(cat1=cat1,
                          cat2=cat2,
                          max_distance=5 * units.arcsec)

    astropy.io.ascii.write(result, outfile_name, format='csv')

def download_sao_xmatch():
    with contextlib.suppress(FileExistsError):
        os.mkdir("xmatch")

    print("SAO-HIP cross match")
    download_xmatch('vizier:I/131A/sao', 'vizier:I/311/hip2', os.path.join('xmatch', 'sao_hip_xmatch.csv'))
    print("SAO-TYC2 cross match")
    download_xmatch('vizier:I/131A/sao', 'vizier:I/259/tyc2', os.path.join('xmatch', 'sao_tyc2_xmatch.csv'))

# --- VIZIER DOWNLOAD ---
def download_vizier():
    with contextlib.suppress(FileExistsError):
        os.mkdir('vizier')

    files_urls = [
        ('ascc.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?I/280B'),
        ('tyc2hd.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?IV/25'),
        ('tyc2spec.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?III/231'),
        ('tyc2teff.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?V/136'),
        ('ubvriteff.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?J/ApJS/193/1'),
        ('xhip.tar.gz', 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Cat/tar.gz?V/137D')]

    for file_name, url in files_urls:
        download_file(os.path.join('vizier', file_name), url)

if __name__ == "__main__":
    download_gaia()
    download_vizier()
    download_sao_xmatch()
