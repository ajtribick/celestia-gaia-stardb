#!/usr/bin/env python3

"""Routines for parsing the HIP data."""

import gzip
import os
import tarfile
import warnings

import numpy as np
import astropy.units as u

from astropy import io
from astropy.coordinates import ICRS, SkyCoord
from astropy.table import join
from astropy.time import Time
from astropy.units import UnitsWarning

from spparse import CEL_UNKNOWN_STAR, parse_spectrum_vec

# temperatures from star.cpp, spectral types O3-M9

TEFF_SPEC = np.array([
    52500, 48000, 44500, 41000, 38000, 35800, 33000,
    30000, 25400, 22000, 18700, 17000, 15400, 14000, 13000, 11900, 10500,
    9520, 9230, 8970, 8720, 8460, 8200, 8020, 7850, 7580, 7390,
    7200, 7050, 6890, 6740, 6590, 6440, 6360, 6280, 6200, 6110,
    6030, 5940, 5860, 5830, 5800, 5770, 5700, 5630, 5570, 5410,
    5250, 5080, 4900, 4730, 4590, 4350, 4200, 4060, 3990, 3920,
    3850, 3720, 3580, 3470, 3370, 3240, 3050, 2940, 2640, 2000])

TEFF_BINS = (TEFF_SPEC[:-1] + TEFF_SPEC[1:]) // 2

CEL_SPECS = parse_spectrum_vec(['OBAFGKM'[i//10]+str(i%10) for i in range(3, 70)])

def load_gaia_hip():
    """Load the Gaia DR2 HIP sources."""
    print('Loading Gaia DR2 sources for HIP')
    col_names = ['source_id', 'hip_id', 'ra', 'dec', 'phot_g_mean_mag', 'bp_rp',
                 'teff_val', 'r_est']
    gaia = io.ascii.read(os.path.join('gaia', 'gaiadr2_hip-result.csv'), include_names=col_names)
    gaia.rename_column('hip_id', 'HIP')
    return gaia

def load_xhip():
    """Load the XHIP catalogue from the VizieR archive."""
    print('Loading XHIP')
    with tarfile.open(os.path.join('vizier', 'xhip.tar.gz'), 'r:gz') as tf:
        print('  Loading main catalog')
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['HIP', 'RAdeg', 'DEdeg', 'Plx', 'pmRA', 'pmDE',
                         'e_Plx', 'Dist', 'e_Dist', 'SpType', 'RV']
            fill_values = [
                ('', '-1', 'Tc', 'Lc'),
                ('', 'NaN', 'phi')
            ]
            reader = io.ascii.get_reader(io.ascii.Cds,
                                         readme=readme,
                                         include_names=col_names,
                                         # workaround for bug in astropy where nullability is
                                         # ignored for columns with limits, see
                                         # https://github.com/astropy/astropy/issues/9291
                                         fill_values=fill_values)
            reader.data.table_name = 'main.dat'
            with tf.extractfile('./main.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                # Suppress a warning generated because the reader does not handle logarithmic units
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UnitsWarning)
                    hip_data = reader.read(f)

        hip_data.add_index('HIP')

        print('  Loading photometric data')
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['HIP', 'Vmag', 'Jmag', 'Hmag', 'Kmag', 'e_Jmag',
                         'e_Hmag', 'e_Kmag', 'B-V', 'V-I', 'e_B-V', 'e_V-I']
            reader = io.ascii.get_reader(io.ascii.Cds,
                                         readme=readme,
                                         include_names=col_names)
            reader.data.table_name = 'photo.dat'
            with tf.extractfile('./photo.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                photo_data = reader.read(f)

        photo_data['HIP'].unit = None # for some reason it is set to parsecs in the ReadMe
        photo_data.add_index('HIP')
        hip_data = join(hip_data, photo_data, join_type='left', keys='HIP')

        print('  Loading bibliographic data')
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['HIP', 'HD']
            reader = io.ascii.get_reader(io.ascii.Cds,
                                         readme=readme,
                                         include_names=col_names)
            reader.data.table_name = 'biblio.dat'
            with tf.extractfile('./biblio.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                biblio_data = reader.read(f)

        biblio_data.add_index('HIP')

        return join(hip_data, biblio_data, join_type='left', keys='HIP')

def load_ubvri():
    """Load UBVRI Teff calibration from VizieR archive."""
    print('Loading UBVRI calibration')
    with tarfile.open(os.path.join('vizier', 'ubvriteff.tar.gz'), 'r:gz') as tf:
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['V-K', 'B-V', 'V-I', 'J-K', 'H-K', 'Teff']
            reader = io.ascii.get_reader(io.ascii.Cds,
                                         readme=readme,
                                         include_names=col_names)
            reader.data.table_name = 'table3.dat'
            with tf.extractfile('./table3.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                # Suppress a warning generated because the reader does not handle logarithmic units
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UnitsWarning)
                    return reader.read(f)

def compute_distances(hip_data, length_kpc=1.35):
    """Compute the distance using an exponentially-decreasing prior.

    The method is described in:

    Bailer-Jones (2015) "Estimating Distances from Parallaxes"
    https://ui.adsabs.harvard.edu/abs/2015PASP..127..994B/abstract

    Using a uniform length scale of 1.35 kpc as suggested for the TGAS
    catalogue of stars in Gaia DR1.

    Astraatmadja & Bailer-Jones "Estimating distances from parallaxes.
    III. Distances of two million stars in the Gaia DR1 catalogue"
    https://ui.adsabs.harvard.edu/abs/2016ApJ...833..119A/abstract
    """

    print('Computing distances')

    eplx2 = hip_data['e_Plx'] ** 2

    r3coeff = np.full_like(hip_data['Plx'], 1/length_kpc)
    r2coeff = np.full_like(hip_data['Plx'], -2)

    roots = np.apply_along_axis(np.roots,
                                0,
                                [r3coeff, r2coeff, hip_data['Plx'] / eplx2, -1 / eplx2])
    roots[np.logical_or(np.real(roots) < 0.0, abs(np.imag(roots)) > 1.0e-6)] = np.nan
    parallax_distance = np.nanmin(np.real(roots), 0) * 1000

    # prefer cluster distances (e_Dist NULL), otherwise use computed distance
    is_cluster_distance = np.logical_and(np.logical_not(hip_data['Dist'].mask),
                                         hip_data['e_Dist'].mask)

    hip_data['dist_use'] = np.where(is_cluster_distance, hip_data['Dist'], parallax_distance)

HIP_TIME = Time('J1991.25')
GAIA_TIME = Time('J2015.5')

def transform_coord(ra, dec, pmra, pmdec, dist, rv):
    """Update the coordinate of one star from J1991.25 to J2015.5 to match Gaia."""
    coord = SkyCoord(frame=ICRS,
                     ra=ra*u.deg,
                     dec=dec*u.deg,
                     pm_ra_cosdec=pmra*u.mas/u.year,
                     pm_dec=pmdec*u.mas/u.year,
                     distance=dist*u.pc,
                     radial_velocity=rv*u.km/u.s,
                     obstime=HIP_TIME).apply_space_motion(GAIA_TIME)
    return coord.ra / u.rad, coord.dec / u.rad

transform_coord_vec = np.vectorize(transform_coord) # pylint: disable=invalid-name

def update_coordinates(hip_data):
    """Update the coordinates from J1991.25 to J2015.5 to match Gaia."""
    print('Updating coordinates to J2015.5')

    result = transform_coord_vec(hip_data['RAdeg'],
                                 hip_data['DEdeg'],
                                 hip_data['pmRA'],
                                 hip_data['pmDE'],
                                 hip_data['dist_use'],
                                 hip_data['RV'].filled(0))

    hip_data['RArad'] = result[0]
    hip_data['DErad'] = result[1]

def parse_spectra(hip_data):
    """Parse the spectral types into the celestia.Sci format."""
    print('Parsing spectral types')
    hip_data['CelSpec'] = parse_spectrum_vec(hip_data['SpType'].filled(''))

def estimate_temperature(ubvri_data, b_v, e_bv, v_i, e_vi, v_k, e_vk, j_k, e_jk, h_k, e_hk):
    """Estimate the temperature of one star from color indices."""
    iweights = [((ubvri_data['B-V'] - b_v) / max(abs(e_bv), 0.0001)) ** 2,
                ((ubvri_data['V-I'] - v_i) / max(abs(e_vi), 0.0001)) ** 2,
                ((ubvri_data['V-K'] - v_k) / max(abs(e_vk), 0.0001)) ** 2,
                ((ubvri_data['J-K'] - j_k) / max(abs(e_jk), 0.0001)) ** 2,
                ((ubvri_data['H-K'] - h_k) / max(abs(e_hk), 0.0001)) ** 2]
    weights = 1.0 / np.maximum(np.nan_to_num(np.nansum(iweights, axis=0), nan=np.inf), 0.0001)
    return np.sum(ubvri_data['Teff'] * weights) / np.sum(weights)

estimate_temperature_vec = np.vectorize(estimate_temperature, excluded={0}) # pylint: disable=invalid-name

def estimate_temperatures(hip_data, ubvri_data):
    """Estimate the temperatures from color indices."""
    print('Estimating temperatures from color indices')
    vmag = hip_data['Vmag'].filled(np.nan)
    jmag = hip_data['Jmag'].filled(np.nan)
    e_jmag = hip_data['e_Jmag'].filled(np.nan)
    hmag = hip_data['Hmag'].filled(np.nan)
    e_hmag = hip_data['e_Hmag'].filled(np.nan)
    kmag = hip_data['Kmag'].filled(np.nan)
    e_kmag = hip_data['e_Kmag'].filled(np.nan)
    hip_data['Teff'] = estimate_temperature_vec(ubvri_data,
                                                hip_data['B-V'].filled(np.nan),
                                                hip_data['e_B-V'].filled(np.nan),
                                                hip_data['V-I'].filled(np.nan),
                                                hip_data['e_V-I'].filled(np.nan),
                                                vmag-kmag,
                                                np.sqrt(e_kmag**2 + 0.05**2),
                                                jmag-kmag,
                                                np.sqrt(e_jmag**2 + e_kmag**2),
                                                hmag-kmag,
                                                np.sqrt(e_hmag**2 + e_kmag**2))

    # use estimated temperatures to fill in missing spectral types
    est_types = CEL_SPECS[np.digitize(hip_data['Teff'], TEFF_BINS)]
    mask = hip_data['CelSpec'] == CEL_UNKNOWN_STAR
    hip_data['CelSpec'][mask] = est_types[mask]

def process_hip(hip_data, ubvri_data):
    """Post-process the HIP data."""
    compute_distances(hip_data)
    parse_spectra(hip_data)
    update_coordinates(hip_data)
    estimate_temperatures(hip_data, ubvri_data)

if __name__ == '__main__':
    #hip_data = load_xhip()
    #ubvri_data = load_ubvri()
    #process_hip(hip_data, ubvri_data)
    #print(hip_data)
    gaia_data = join(load_gaia_hip(),
                     load_xhip(),
                     keys=['HIP'],
                     join_type='outer',
                     table_names=['gaia', 'xhip'])
    print(gaia_data)
