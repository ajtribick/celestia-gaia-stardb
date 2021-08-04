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

"""Routines for downloading data for the HIP2 distance estimation."""

import astropy.units as u
from astropy_healpix import HEALPix
import numpy as np

from .directories import HIP2_DIST_DIR, VIZIER_DIR
from .utils import confirm_action, download_file, open_cds_tarfile
from .celestia_gaia import estimate_distances


_PRIORS_FILE = HIP2_DIST_DIR/'prior_summary.csv'
_HEALPIX_FILE = HIP2_DIST_DIR/'hip2healpix.csv'
_DIST_FILE = HIP2_DIST_DIR/'hip2dist.csv'


def download_dist_prior() -> None:
    """Download the Gaia EDR3 distance prior."""
    HIP2_DIST_DIR.mkdir(parents=True, exist_ok=True)
    download_file(
        _PRIORS_FILE,
        'https://keeper.mpdl.mpg.de/d/0e9919daa8fa4e119ebe/files/?p=%2Fprior_summary.csv&dl=1'
    )


def apply_healpix() -> None:
    """Determines the HEALPix containing each XHIP star."""
    if (
        _HEALPIX_FILE.exists()
        and not confirm_action('XHIP-HEALPix mapping already generated, replace?')
    ):
        return

    print('Mapping XHIP star positions to HEALPix cells')
    with open_cds_tarfile(VIZIER_DIR/'xhip.tar.gz') as tf:
        data = tf.read_gzip(
            'main.dat',
            ['HIP', 'RAdeg', 'DEdeg', 'Plx', 'e_Plx']
        )

    hp = HEALPix(nside=1<<5, order='nested')
    data['healpix'] = np.apply_along_axis(
        lambda ra_dec: hp.lonlat_to_healpix(ra_dec[0]*u.deg, ra_dec[1]*u.deg),
        0,
        [data['RAdeg'], data['DEdeg']]
    )

    data.write(_HEALPIX_FILE, format='csv', overwrite=True)


def build_hip2_distances() -> None:
    """Estimates distances using the Gaia EDR3 geometric distance algorithm."""
    if (
        _DIST_FILE.exists()
        and not confirm_action('HIP2 distances file already generated, replace?')
    ):
        return

    estimate_distances(_PRIORS_FILE, _HEALPIX_FILE, _DIST_FILE)
