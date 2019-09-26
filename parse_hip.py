#!/usr/bin/python3

import gzip, os, tarfile, warnings
import numpy as np
import astropy.units as u

from astropy.coordinates import ICRS, SkyCoord
from astropy.io import ascii
from astropy.table import join
from astropy.time import Time
from astropy.units import UnitsWarning
from spparse import parse_spectrum

def load_xhip():
    """Load the XHIP catalogue from the VizieR archive."""
    print('Loading XHIP catalogue')
    with tarfile.open(os.path.join('vizier', 'xhip.tar.gz'), 'r:gz') as tf:
        print("Loading main catalog")
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['HIP', 'RAdeg', 'DEdeg', 'Plx', 'pmRA', 'pmDE', 'e_Plx', 'Dist', 'e_Dist', 'SpType', 'RV']
            reader = ascii.get_reader(ascii.Cds,
                                      readme=readme,
                                      include_names=col_names,
                                      # workaround for bug in astropy where nullability is ignored for columns with limits
                                      # see https://github.com/astropy/astropy/issues/9291
                                      fill_values=[('', '-1', 'Tc', 'Lc'),
                                                   ('', 'NaN', 'phi')])
            reader.data.table_name = 'main.dat'
            with tf.extractfile('./main.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                # Suppress a warning generated because the reader does not handle logarithmic units
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UnitsWarning)
                    hip_data = reader.read(f)

        hip_data.add_index('HIP')

        print('Loading photometric data')
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['HIP', 'Vmag']
            reader = ascii.get_reader(ascii.Cds,
                                      readme=readme,
                                      include_names=col_names)
            reader.data.table_name = 'photo.dat'
            with tf.extractfile('./photo.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                photo_data = reader.read(f)

        photo_data['HIP'].unit = None # for some reason it is set to parsecs in the ReadMe
        photo_data.add_index('HIP')
        hip_data = join(hip_data, photo_data, join_type='left', keys='HIP')

        print('Loading bibliographic data')
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['HIP', 'HD']
            reader = ascii.get_reader(ascii.Cds,
                                      readme=readme,
                                      include_names=col_names)
            reader.data.table_name = 'biblio.dat'
            with tf.extractfile('./biblio.dat.gz') as gzf, gzip.open(gzf, 'rb') as f:
                biblio_data = reader.read(f)

        biblio_data.add_index('HIP')

        return join(hip_data, biblio_data, join_type='left', keys='HIP')


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

    roots = np.apply_along_axis(np.roots, 0, [r3coeff, r2coeff, hip_data['Plx'] / eplx2, -1 / eplx2])
    roots[np.logical_or(np.real(roots) < 0.0, abs(np.imag(roots)) > 1.0e-6)] = np.inf
    parallax_distance = np.amin(np.real(roots), 0) * 1000

    # prefer cluster distances (e_Dist NULL), otherwise use computed distance
    is_cluster_distance = np.logical_and(np.logical_not(hip_data['Dist'].mask),
                                         hip_data['e_Dist'].mask)

    hip_data['dist_use'] = np.where(is_cluster_distance, hip_data['Dist'], parallax_distance)

hip_time = Time('J1991.25')
gaia_time = Time('J2015.5')

def transform_coord(ra, dec, pmra, pmdec, dist, rv):
    coord = SkyCoord(frame=ICRS,
                     ra=ra*u.deg,
                     dec=dec*u.deg,
                     pm_ra_cosdec=pmra*u.mas/u.year,
                     pm_dec=pmdec*u.mas/u.year,
                     distance=dist*u.pc,
                     radial_velocity=rv*u.km/u.s,
                     obstime=hip_time).apply_space_motion(gaia_time)
    return coord.ra / u.deg, coord.dec / u.deg

transform_coord_vec = np.vectorize(transform_coord)

def update_coordinates(hip_data):
    """Update the coordinates from J1991.25 to J2015.5 to match Gaia."""
    print('Updating coordinates to J2015.5')

    result = transform_coord_vec(hip_data['RAdeg'],
                                 hip_data['DEdeg'],
                                 hip_data['pmRA'],
                                 hip_data['pmDE'],
                                 hip_data['dist_use'],
                                 hip_data['RV'].filled(0))

    hip_data['RAdeg2015'] = result[0]
    hip_data['DEdeg2015'] = result[1]

parse_spectrum_vec = np.vectorize(parse_spectrum)

def parse_spectra(hip_data):
    """Parse the spectral types into the celestia.Sci format."""
    print('Parsing spectral types')
    hip_data['CelSpec'] = parse_spectrum_vec(hip_data['SpType'].filled(''))

def process_hip(hip_data):
    compute_distances(hip_data)
    update_coordinates(hip_data)
    parse_spectra(hip_data)

if __name__ == '__main__':
    hip_data = load_xhip()
    process_hip(hip_data)
    print(hip_data)
