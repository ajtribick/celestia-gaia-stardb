#!/usr/bin/env python3

import gzip, os, tarfile
import numpy as np
import astropy.units as u

from astropy.io import ascii
from astropy.table import join, unique, vstack

def load_gaia_tyc():
    """Load the Gaia DR2 TYC2 sources."""
    print('Loading Gaia DR2 sources for TYC2')
    col_names = ['source_id', 'tyc2_id', 'ra', 'dec', 'phot_g_mean_mag', 'bp_rp', 'teff_val', 'r_est']
    t = ascii.read(os.path.join('gaia', 'gaiadr2_tyc-result.csv'), include_names=col_names)
    tycs = np.array(np.char.split(t['tyc2_id'], '-').tolist()).astype(np.int64)
    t['TYC'] = tycs[:,0] + tycs[:,1]*10000 + tycs[:,2]*1000000000
    t.remove_column('tyc2_id')
    return t

def load_ascc():
    """Load ASCC from VizieR archive."""
    def load_section(tf, info):
        with tf.extractfile('./ReadMe') as readme:
            col_names = ['Bmag', 'Vmag', 'e_Bmag', 'e_Vmag', 'd3', 'TYC1',
                        'TYC2', 'TYC3', 'HD', 'Jmag', 'e_Jmag', 'Hmag',
                        'e_Hmag', 'Kmag', 'e_Kmag']
            reader = ascii.get_reader(ascii.Cds,
                                    readme=readme,
                                    include_names=col_names)
            reader.data.table_name = 'cc*.dat'
            print('  Loading ' + os.path.basename(info.name))
            with tf.extractfile(info) as gzf, gzip.open(gzf, 'rb') as f:
                section = reader.read(f)

        section = section[section['TYC1'] != 0]

        section['TYC'] = (section['TYC1'] +
                          section['TYC2']*10000 +
                          section['TYC3']*1000000000)

        section.remove_columns(['TYC1', 'TYC2', 'TYC3'])

        convert_cols = ['Bmag', 'Vmag', 'e_Bmag', 'e_Vmag', 'Jmag',
                        'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag']
        for col in convert_cols:
            section[col] = section[col].astype(np.float64)
            section[col].convert_unit_to(u.mag)

        return section

    def is_data(info):
        sections = os.path.split(info.name)
        return (len(sections) == 2 and
                sections[0] == '.' and
                sections[1].startswith('cc'))

    print('Loading ASCC')
    with tarfile.open(os.path.join('vizier', 'ascc.tar.gz'), 'r:gz') as tf:
        combined = vstack([load_section(tf, m) for m in tf if is_data(m)],
                          join_type='exact')

    return unique(combined.group_by(['TYC','d3']), keys=['TYC'])

