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

from pathlib import Path

import astropy.io.ascii as io_ascii
import requests
from astropy import units
from astroquery.xmatch import XMatch

from .directories import VIZIER_DIR, XMATCH_DIR
from .utils import confirm_action


def _proceed_checkfile(path: Path) -> bool:
    """Check if a file exists, if so prompt the user if they want to replace it."""
    if path.exists():
        if confirm_action(f'{path} already exists, replace?'):
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
