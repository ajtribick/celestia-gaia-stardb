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

"""Routines for parsing the HIP data."""

import warnings

import astropy.io.ascii as io_ascii
import astropy.units as u
import numpy as np
from astropy.coordinates import ICRS, SkyCoord
from astropy.table import MaskedColumn, Table, join, unique
from astropy.time import Time
from erfa import ErfaWarning

from .celestia_gaia import apply_distances
from .directories import GAIA_EDR3_DIR, HIP2_DIST_DIR, VIZIER_DIR, XMATCH_DIR
from .utils import open_cds_tarfile


def load_gaia_hip() -> Table:
    """Load the Gaia DR2 HIP sources."""
    print('Loading Gaia DR2 sources for HIP')

    gaia = Table.read(GAIA_EDR3_DIR/'xmatch-gaia-hip.vot.gz', format='votable')
    gaia['bp_rp'] = gaia['phot_bp_mean_mag'] - gaia['phot_rp_mean_mag']

    print('  Merging distance information')
    gaia['hip'] = gaia['hip'].astype('uint32')
    gaia['r_est'] = MaskedColumn(apply_distances(GAIA_EDR3_DIR, gaia['source_id']), unit='pc')
    gaia['r_est'].mask = np.isnan(gaia['r_est'])

    gaia.rename_column('hip', 'HIP')
    gaia.remove_columns([
        'hip_ra', 'hip_dec', 'hp_mag', 'dr2_radial_velocity', 'ref_epoch', 'pmra', 'pmra_error',
        'pmdec', 'pmdec_error', 'phot_bp_mean_mag', 'phot_rp_mean_mag',
        'astrometric_params_solved',
    ])
    return gaia


def load_xhip() -> Table:
    """Load the XHIP catalogue from the VizieR archive."""
    print('Loading XHIP')
    with open_cds_tarfile(VIZIER_DIR/'xhip.tar.gz') as tf:
        print('  Loading main catalog')
        hip_data = tf.read_gzip(
            'main.dat',
            [
                'HIP', 'Comp', 'RAdeg', 'DEdeg', 'Plx', 'pmRA', 'pmDE', 'e_Plx',
                'Dist', 'e_Dist', 'SpType', 'RV',
            ],
            fill_values=[('', '-1', 'Tc', 'Lc'), ('', 'NaN', 'phi')],
        )
        hip_data.add_index('HIP')

        print('  Loading photometric data')
        photo_data = tf.read_gzip(
            'photo.dat',
            [
                'HIP', 'Vmag', 'Jmag', 'Hmag', 'Kmag', 'e_Jmag', 'e_Hmag', 'e_Kmag',
                'B-V', 'V-I', 'e_B-V', 'e_V-I',
            ],
        )
        photo_data['HIP'].unit = None # for some reason it is set to parsecs in the ReadMe
        photo_data.add_index('HIP')
        hip_data = join(hip_data, photo_data, join_type='left', keys='HIP')

        print('  Loading bibliographic data')
        biblio_data = tf.read_gzip('biblio.dat', ['HIP', 'HD'])
        biblio_data.add_index('HIP')
        return join(hip_data, biblio_data, join_type='left', keys='HIP')


def load_tyc2specnew() -> Table:
    """Load revised spectral types."""
    print("Loading revised TYC2 spectral types")
    with open_cds_tarfile(VIZIER_DIR/'tyc2specnew.tar.gz') as tf:
        data = tf.read('table2.dat', ['HIP', 'SpType1'])
        return data[data['SpType1'] != '']


def load_sao() -> Table:
    """Load the SAO-HIP cross match."""
    print('Loading SAO-HIP cross match')
    data = io_ascii.read(
        XMATCH_DIR/'sao_hip_xmatch.csv',
        include_names=['HIP', 'SAO', 'angDist', 'delFlag'],
        format='csv',
    )

    data = data[data['delFlag'].mask]
    data.remove_column('delFlag')

    data = unique(data.group_by(['HIP', 'angDist']), keys=['HIP'])
    data.remove_column('angDist')

    data.add_index('HIP')
    return data


def load_distances() -> Table:
    """Loads the computed HIP2 distances."""

    print('Loading distances')
    dist = io_ascii.read(
        HIP2_DIST_DIR/'hip2dist.csv',
        include_names=['HIP', 'dist_med'],
        format='csv',
    )

    return dist


HIP_TIME = Time('J1991.25')
GAIA_TIME = Time('J2015.5')


def update_coordinates(hip_data: Table) -> None:
    """Update the coordinates from J1991.25 to J2015.5 to match Gaia."""
    print('Updating coordinates to J2015.5')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', ErfaWarning)
        coords = SkyCoord(
            frame=ICRS,
            ra=hip_data['RAdeg'],
            dec=hip_data['DEdeg'],
            pm_ra_cosdec=hip_data['pmRA'],
            pm_dec=hip_data['pmDE'],
            distance=hip_data['r_est'],
            radial_velocity=hip_data['RV'].filled(0),
            obstime=HIP_TIME,
        ).apply_space_motion(GAIA_TIME)

    hip_data['ra'] = coords.ra / u.deg
    hip_data['ra'].unit = u.deg
    hip_data['dec'] = coords.dec / u.deg
    hip_data['dec'].unit = u.deg


def process_xhip() -> Table:
    """Processes the XHIP data."""
    xhip = load_xhip()
    sptypes = load_tyc2specnew()
    xhip = join(xhip, sptypes, keys=['HIP'], join_type='left', metadata_conflicts='silent')
    xhip['SpType'] = xhip['SpType1'].filled(xhip['SpType'])
    xhip.remove_column('SpType1')

    xhip = join(
        xhip, load_distances(), keys=['HIP'], join_type='left', metadata_conflicts='silent'
    )

    # prefer cluster distances (e_Dist NULL), otherwise use computed distance
    is_cluster_distance = np.logical_and(
        np.logical_not(xhip['Dist'].mask),
        xhip['e_Dist'].mask,
    )
    xhip['r_est'] = np.where(is_cluster_distance, xhip['Dist'], xhip['dist_med'])
    xhip['r_est'].unit = u.pc
    xhip.remove_column('dist_med')

    update_coordinates(xhip)
    xhip.remove_columns(['RAdeg', 'DEdeg', 'pmRA', 'pmDE', 'RV', 'Dist', 'e_Dist'])
    return xhip


def process_hip() -> Table:
    """Process the Gaia and HIP data."""
    data = join(
        load_gaia_hip(),
        process_xhip(),
        keys=['HIP'],
        join_type='outer',
        table_names=['gaia', 'xhip'],
        metadata_conflicts='silent',
    )

    data = join(data, load_sao(), keys=['HIP'], join_type='left')

    data['r_gaia_score'] = np.where(
        data['r_est_gaia'].mask,
        -20000.0,
        np.where(data['parallax'] <= 0, -10000.0, data['parallax'] / data['parallax_error']),
    )
    data['r_xhip_score'] = np.where(data['Plx'] <= 0, -10000.0, data['Plx'] / data['e_Plx'])

    data['dist_use'] = np.where(
        data['r_gaia_score'] >= data['r_xhip_score'], data['r_est_gaia'], data['r_est_xhip']
    )
    data['dist_use'].unit = u.pc

    data['ra'] = data['ra_gaia'].filled(data['ra_xhip'])
    data['dec'] = data['dec_gaia'].filled(data['dec_xhip'])

    data.remove_columns([
        'ra_gaia', 'dec_gaia', 'r_est_gaia', 'ra_xhip', 'dec_xhip', 'r_est_xhip',
        'parallax', 'parallax_error', 'Plx', 'e_Plx',
        'r_gaia_score', 'r_xhip_score',
    ])

    return data
