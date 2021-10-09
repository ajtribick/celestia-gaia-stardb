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

"""Makes the star database."""

import gzip
import struct
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

import astropy.units as u
import numpy as np
from astropy.table import MaskedColumn, Table, join, unique, vstack

from .celestia_gaia import apply_distances
from .directories import GAIA_EDR3_DIR, OUTPUT_DIR, VIZIER_DIR
from .parse_hip import process_hip
from .parse_tyc import make_tyc, process_tyc
from .spparse import CEL_UNKNOWN_STAR, parse_spectrum
from .utils import WorkaroundCDSReader, open_cds_tarfile

VERSION = "1.1.0-beta.1"

# remove the following objects from the output

EXCLUSIONS = [
    60936,  # quasar 3C 273
    114110, # non-existent star (see HIP1 errata)
    114176, # non-existent star (see HIP1 errata)
]

# temperatures from star.cpp, spectral types O3-M9

TEFF_SPEC = np.array([
    52500, 48000, 44500, 41000, 38000, 35800, 33000,
    30000, 25400, 22000, 18700, 17000, 15400, 14000, 13000, 11900, 10500,
    9520, 9230, 8970, 8720, 8460, 8200, 8020, 7850, 7580, 7390,
    7200, 7050, 6890, 6740, 6590, 6440, 6360, 6280, 6200, 6110,
    6030, 5940, 5860, 5830, 5800, 5770, 5700, 5630, 5570, 5410,
    5250, 5080, 4900, 4730, 4590, 4350, 4200, 4060, 3990, 3920,
    3850, 3720, 3580, 3470, 3370, 3240, 3050, 2940, 2640, 2000,
])

TEFF_BINS = (TEFF_SPEC[:-1] + TEFF_SPEC[1:]) // 2

parse_spectrum_vec = np.vectorize(parse_spectrum, otypes=[np.uint16]) # pylint: disable=invalid-name

CEL_SPECS = parse_spectrum_vec(['OBAFGKM'[i//10]+str(i%10) for i in range(3, 70)])


def load_ubvri() -> Table:
    """Load UBVRI Teff calibration from VizieR archive."""
    print('Loading UBVRI calibration')
    with open_cds_tarfile(VIZIER_DIR/'ubvriteff.tar.gz') as tf:
        return tf.read_gzip('table3.dat', ['V-K', 'B-V', 'V-I', 'J-K', 'H-K', 'Teff'])


def parse_spectra(data: Table) -> Table:
    """Parse the spectral types into the celestia.Sci format."""
    print('Parsing spectral types')
    data['SpType'] = data['SpType'].filled('')
    sptypes = unique(data['SpType',])
    sptypes['CelSpec'] = parse_spectrum_vec(sptypes['SpType'])
    return join(data, sptypes)


def estimate_magnitudes(data: Table) -> None:
    """Estimates magnitudes and color indices from G magnitude and BP-RP.

    Formula used is from Riello et al. (2021) "Gaia Early Data Release 3: Photometric content and
    validation".
    """
    print("Computing missing magnitudes and color indices")

    bp_rp = data['bp_rp'].filled(0)

    data['Vmag'] = MaskedColumn(
        data['Vmag'].filled(
            data['phot_g_mean_mag'].filled(np.nan)
            - (-0.02704 + bp_rp*(0.01424 + bp_rp*(-0.2156 + bp_rp*0.01426))),
        )
    )
    data['e_Vmag'] = MaskedColumn(data['e_Vmag'].filled(0.03017))
    data['Vmag'].mask = np.isnan(data['Vmag'])
    data['e_Vmag'].mask = data['Vmag'].mask

    bp_rp = data['bp_rp'].filled(np.nan)

    imag = (
        data['phot_g_mean_mag'].filled(np.nan)
        - (0.01753 + bp_rp*(0.76 - 0.0991*bp_rp))
    )
    e_imag = np.where(np.isnan(imag), np.nan, 0.03765)

    f_bmag = data['Bmag'].filled(np.nan)
    f_vmag = data['Vmag'].filled(np.nan)
    f_jmag = data['Jmag'].filled(np.nan)
    f_hmag = data['Hmag'].filled(np.nan)
    f_kmag = data['Kmag'].filled(np.nan)
    f_e_bmag = data['e_Bmag'].filled(np.nan)
    f_e_vmag = data['e_Vmag'].filled(np.nan)
    f_e_jmag = data['e_Jmag'].filled(np.nan)
    f_e_hmag = data['e_Hmag'].filled(np.nan)
    f_e_kmag = data['e_Kmag'].filled(np.nan)

    data['B-V'] = MaskedColumn(data['B-V'].filled(f_bmag - f_vmag))
    data['e_B-V'] = MaskedColumn(data['e_B-V'].filled(np.sqrt(f_e_bmag**2 + f_e_vmag**2)))
    data['V-I'] = MaskedColumn(data['V-I'].filled(f_vmag - imag))
    data['e_V-I'] = MaskedColumn(data['e_V-I'].filled(np.sqrt(f_e_vmag**2 + e_imag**2)))
    data['V-K'] = MaskedColumn(f_vmag - f_kmag)
    data['e_V-K'] = MaskedColumn(np.sqrt(f_e_vmag**2 + f_e_kmag**2))
    data['J-K'] = MaskedColumn(f_jmag - f_kmag)
    data['e_J-K'] = MaskedColumn(np.sqrt(f_e_jmag**2 + f_e_kmag**2))
    data['H-K'] = MaskedColumn(f_hmag - f_kmag)
    data['e_H-K'] = MaskedColumn(np.sqrt(f_e_hmag**2 + f_e_kmag**2))

    data['B-V'].mask = np.logical_or(data['B-V'].mask, np.isnan(data['B-V']))
    data['e_B-V'].mask = np.logical_or(data['e_B-V'].mask, np.isnan(data['e_B-V']))
    data['V-I'].mask = np.logical_or(data['V-I'].mask, np.isnan(data['V-I']))
    data['e_V-I'].mask = np.logical_or(data['e_V-I'].mask, np.isnan(data['e_V-I']))
    data['V-K'].mask = np.isnan(data['V-K'])
    data['e_V-K'].mask = np.isnan(data['e_V-K'])
    data['J-K'].mask = np.isnan(data['J-K'])
    data['e_J-K'].mask = np.isnan(data['e_J-K'])
    data['H-K'].mask = np.isnan(data['H-K'])
    data['e_H-K'].mask = np.isnan(data['e_H-K'])

    data.remove_columns([
        'Bmag', 'e_Bmag', 'e_Vmag', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag',
    ])


def estimate_temperatures(data: Table) -> None:
    """Estimate the temperature of stars."""
    ubvri_data = load_ubvri()
    print('Estimating temperatures from color indices')

    indices = Table(
        [
            data['B-V'].filled(np.nan),
            data['V-I'].filled(np.nan),
            data['V-K'].filled(np.nan),
            data['J-K'].filled(np.nan),
            data['H-K'].filled(np.nan),
            np.maximum(data['e_B-V'].filled(np.nan), 0.001),
            np.maximum(data['e_V-I'].filled(np.nan), 0.001),
            np.maximum(data['e_V-K'].filled(np.nan), 0.001),
            np.maximum(data['e_J-K'].filled(np.nan), 0.001),
            np.maximum(data['e_H-K'].filled(np.nan), 0.001),
        ],
        names=['B-V','V-I','V-K','J-K','H-K', 'e_B-V','e_V-I','e_V-K','e_J-K','e_H-K'],
    )

    weights = np.full_like(data['HIP'], 0, dtype=np.float64)
    teffs = np.full_like(data['HIP'], 0, dtype=np.float64)
    for row in ubvri_data:
        sumsq = np.maximum(
            np.nan_to_num(((indices['B-V']-row['B-V'])/indices['e_B-V'])**2)
            + np.nan_to_num(((indices['V-K']-row['V-K'])/indices['e_V-K'])**2)
            + np.nan_to_num(((indices['J-K']-row['J-K'])/indices['e_J-K'])**2)
            + np.nan_to_num(((indices['V-I']-row['V-I'])/indices['e_V-I'])**2)
            + np.nan_to_num(((indices['H-K']-row['H-K'])/indices['e_H-K'])**2),
            0.001,
        )
        teffs += row['Teff'] / sumsq
        weights += 1.0 / sumsq

    data['teff_est'] = teffs / weights
    data['teff_est'].unit = u.K


def estimate_spectra(data: Table) -> Table:
    """Estimate the spectral type of stars."""
    no_teff = data[data['teff_val'].mask]
    # temporarily disable no-member error in pylint, as it cannot see the reduce method
    # pylint: disable=no-member
    has_indices = np.logical_and.reduce((
        no_teff['B-V'].mask,
        no_teff['V-I'].mask,
        no_teff['V-K'].mask,
        no_teff['J-K'].mask,
        no_teff['H-K'].mask,
    ))
    # pylint: enable=no-member
    no_teff = no_teff[np.logical_not(has_indices)]
    estimate_temperatures(no_teff)
    data = join(
        data,
        no_teff['HIP', 'teff_est'],
        keys=['HIP'],
        join_type='left',
    )
    data['teff_val'] = data['teff_val'].filled(data['teff_est'].filled(np.nan))
    data = data[np.logical_not(np.isnan(data['teff_val']))]
    data['CelSpec'] = CEL_SPECS[np.digitize(data['teff_val'], TEFF_BINS)]
    return data


def load_sao() -> Table:
    """Loads the SAO catalog."""
    print("Loading SAO")

    with (VIZIER_DIR/'sao.readme').open('r') as readme:
        reader = WorkaroundCDSReader('sao.dat', ['SAO', 'HD'], [np.int64, np.int64], readme)

    with gzip.open(VIZIER_DIR/'sao.dat.gz', 'rt', encoding='ascii') as f:
        data = reader.read(f)

    data = unique(data.group_by('SAO'), keys=['HD'])
    data = unique(data.group_by('HD'), keys=['SAO'])
    return data


def merge_all() -> Table:
    """Merges the HIP and TYC data."""
    print("Loading Gaia crossmatch data")
    data = Table.read(GAIA_EDR3_DIR/'xmatch-gaia-hiptyc.vot.gz', format='votable')
    data.rename_column('id', 'HIP')
    data['HIP'] = np.where(
        data['HIP'] < 1000000,
        data['HIP'],
        make_tyc(data['HIP']//1000000, (data['HIP']//10) % 100000, data['HIP']%10)
    ).astype('uint32')

    data['bp_rp'] = data['phot_bp_mean_mag'] - data['phot_rp_mean_mag']

    data['r_est'] = MaskedColumn(apply_distances(GAIA_EDR3_DIR, data['source_id']), unit='pc')
    data['r_est'].mask = np.isnan(data['r_est'])

    data.remove_columns([
        'src_ra', 'src_dec', 'hp_mag', 'bt_mag', 'vt_mag', 'ep_ra1990', 'ep_de1990',
        'dr2_radial_velocity', 'ref_epoch', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error',
        'phot_bp_mean_mag', 'phot_rp_mean_mag', 'astrometric_params_solved'
    ])

    data = process_hip(data)
    data = process_tyc(data)

    # Merge SAO, preferring the values from the SAO catalogue
    sao = load_sao()
    sao = sao[np.isin(sao['HD'], data[np.logical_not(data['HD'].mask)]['HD'])]
    data['SAO'].mask = np.logical_or(
        data['SAO'].mask,
        np.isin(data['SAO'], sao['SAO']),
    )

    hd_sao = join(
        data[np.logical_not(data['HD'].mask)],
        sao,
        keys=['HD'],
        table_names=['xref', 'sao'],
        join_type='left',
    )
    hd_sao.rename_column('SAO_xref', 'SAO')
    hd_sao['SAO'] = MaskedColumn(
        hd_sao['SAO_sao'].filled(hd_sao['SAO']),
        mask=np.logical_and(hd_sao['SAO'].mask, hd_sao['SAO_sao'].mask),
    )
    hd_sao.remove_column('SAO_sao')

    return vstack([data[data['HD'].mask], hd_sao], join_type='exact')


OBLIQUITY = np.radians(23.4392911)
COS_OBLIQUITY = np.cos(OBLIQUITY)
SIN_OBLIQUITY = np.sin(OBLIQUITY)
ROT_MATRIX = np.array([
    [1, 0, 0],
    [0, COS_OBLIQUITY, SIN_OBLIQUITY],
    [0, -SIN_OBLIQUITY, COS_OBLIQUITY]
])


def process_data() -> Table:
    """Processes the missing data values."""
    data = merge_all()
    data = data[data['dist_use'] > 0]
    data = data[np.isin(data['HIP'], EXCLUSIONS, invert=True)]
    estimate_magnitudes(data)
    data = parse_spectra(data)
    unknown_spectra = data[data['CelSpec'] == CEL_UNKNOWN_STAR][
        'HIP', 'teff_val',
        'B-V', 'e_B-V', 'V-I', 'e_V-I', 'V-K', 'e_V-K', 'J-K', 'e_J-K', 'H-K', 'e_H-K',
    ]
    unknown_spectra = estimate_spectra(unknown_spectra)
    data = join(
        data,

        unknown_spectra['HIP', 'CelSpec'],
        keys=['HIP'],
        join_type='left',
        table_names=['data', 'est'],
    )
    data['CelSpec'] = np.where(
        data['CelSpec_data'] == CEL_UNKNOWN_STAR,
        data['CelSpec_est'].filled(CEL_UNKNOWN_STAR),
        data['CelSpec_data'],
    )
    data.remove_columns([
        'phot_g_mean_mag', 'bp_rp', 'teff_val', 'SpType', 'B-V', 'e_B-V', 'V-I', 'e_V-I', 'V-K',
        'e_V-K', 'J-K', 'e_J-K', 'H-K', 'e_H-K', 'CelSpec_est', 'CelSpec_data',
    ])

    data['Vmag_abs'] = data['Vmag'] - 5*(np.log10(data['dist_use'])-1)

    print('Converting coordinates to ecliptic frame')

    data['ra'].convert_unit_to(u.rad)
    data['dec'].convert_unit_to(u.rad)
    data['dist_use'].convert_unit_to(u.lyr)

    coords = np.matmul(
        ROT_MATRIX,
        np.array([
            data['dist_use']*np.cos(data['ra'])*np.cos(data['dec']),
            data['dist_use']*np.sin(data['dec']),
            -data['dist_use']*np.sin(data['ra'])*np.cos(data['dec']),
        ])
    )
    data['x'] = coords[0]
    data['y'] = coords[1]
    data['z'] = coords[2]

    data['x'].unit = u.lyr
    data['y'].unit = u.lyr
    data['z'].unit = u.lyr

    return data


def write_starsdat(data: Table, outfile: Path) -> None:
    """Write the stars.dat file."""
    print('Writing stars.dat')
    with outfile.open('wb') as f:
        f.write(struct.pack('<8sHL', b'CELSTARS', 0x0100, len(data)))
        print(f'  Writing {len(data)} records')
        fmt = struct.Struct('<L3fhH')
        for hip, x, y, z, vmag_abs, celspec in zip(
            data['HIP'], data['x'], data['y'], data['z'], data['Vmag_abs'], data['CelSpec'],
        ):
            f.write(fmt.pack(hip, x, y, z, int(round(vmag_abs*256)), celspec))


def write_xindex(data: Table, field: str, outfile: Path) -> None:
    """Write a cross-index file."""
    print('Writing '+field+' cross-index')

    print('  Extracting cross-index data')
    data = data[np.logical_not(data[field].mask)]['HIP', 'Comp', field]
    data['Comp'] = data['Comp'].filled('')
    data = unique(data.group_by([field, 'Comp', 'HIP']), keys=[field])

    print(f'  Writing {len(data)} records')
    with outfile.open('wb') as f:
        f.write(struct.pack('<8sH', b'CELINDEX', 0x0100))
        fmt = struct.Struct('<2L')
        for hip, cat in zip(data['HIP'], data[field]):
            f.write(fmt.pack(cat, hip))


def make_stardb() -> None:
    """Make the Celestia star database files."""
    data = process_data()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    write_starsdat(data, OUTPUT_DIR/'stars.dat')

    xindices = [
        ('HD', 'hdxindex.dat'),
        ('SAO', 'saoxindex.dat'),
    ]

    for fieldname, outfile in xindices:
        write_xindex(data, fieldname, OUTPUT_DIR/outfile)

    print("Creating archive")
    archivename = f'celestia-gaia-stardb-{VERSION}'
    with ZipFile(f'{archivename}.zip', 'w', compression=ZIP_DEFLATED, compresslevel=9) as zf:
        contents = ['stars.dat', 'hdxindex.dat', 'saoxindex.dat', 'LICENSE.txt', 'CREDITS.md']
        for f in contents:
            zf.write(OUTPUT_DIR/f, arcname=Path(archivename)/f)

    # archivename = f'celestia-gaia-auxiliary-{VERSION}'
    # with ZipFile(f'{archivename}.zip', 'w', compression=ZIP_DEFLATED, compresslevel=9) as zf:
    #     contents = [
    #         'hip2dist.csv', 'hip-gaia-xmatch.csv', 'tyc-gaia-xmatch.csv',
    #         'LICENSE.txt', 'CREDITS.md',
    #     ]
    #     for f in contents:
    #         zf.write(AUXFILES_DIR/f, arcname=Path(archivename)/f)
