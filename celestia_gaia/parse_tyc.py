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

"""Routines for parsing the TYC2 data."""

from __future__ import annotations

import gzip
import tarfile
from pathlib import PurePath
from typing import IO

import astropy.io.ascii as io_ascii
import astropy.units as u
import numpy as np
from astropy.table import MaskedColumn, Table, join, unique, vstack

from .directories import VIZIER_DIR, XMATCH_DIR
from .utils import TarCds, WorkaroundCDSReader, open_cds_tarfile


def make_tyc(tyc1: int, tyc2: int, tyc3: int) -> int:
    """Build a synthetic HIP identifier from TYC parts."""
    return tyc1 + tyc2 * 10000 + tyc3 * 1000000000


TYC_HD_ERRATA = {
    "add": [
        # B. Skiff, 30-Jan-2007
        (make_tyc(8599, 1797, 1), 298954),
        # B. Skiff, 12-Jul-2007
        (make_tyc(6886, 1389, 1), 177868),
        # LMC inner region
        (make_tyc(9161, 685, 1), 269051),
        (make_tyc(9165, 548, 1), 269052),
        (make_tyc(9169, 1563, 1), 269207),
        (make_tyc(9166, 2, 1), 269367),
        (make_tyc(9166, 540, 1), 269382),
        (make_tyc(8891, 278, 1), 269537),
        (make_tyc(9162, 657, 1), 269599),
        (make_tyc(9167, 730, 1), 269858),
        (make_tyc(9163, 960, 2), 269928),
        (make_tyc(9163, 751, 1), 270005),
        (make_tyc(8904, 686, 1), 270078),
        (make_tyc(9163, 887, 1), 270128),
        (make_tyc(8904, 766, 1), 270342),
        (make_tyc(9168, 1217, 1), 270435),
        (make_tyc(8904, 911, 1), 270467),
        (make_tyc(8904, 5, 1), 270485),
        # LMC outer region
        (make_tyc(9157, 1, 1), 270502),
        (make_tyc(9160, 1142, 1), 270526),
        (make_tyc(8888, 928, 1), 270765),
        (make_tyc(8888, 910, 1), 270794),
        (make_tyc(9172, 791, 1), 272092),
        # VizieR annotations
        (make_tyc(7389, 1138, 1), 320669),
        # extras
        (make_tyc(1209, 1833, 1), 11502),
        (make_tyc(1209, 1835, 1), 11503),
    ],
    "delete": [
        # B. Skiff (13-Nov-2007)
        32228,
        # LMC inner region
        269686,
        269784,
        # LMC outer region
        270653,
        271058,
        271224,
        271264,
        271389,
        271600,
        271695,
        271727,
        271764,
        271802,
        271875,
        # VizieR annotations
        181060,
    ],
}


def parse_tyc_cols(
    data: Table,
    src_columns: tuple[str, str, str] = ("TYC1", "TYC2", "TYC3"),
    dest_column: str = "HIP",
) -> None:
    """Convert TYC identifier components into a synthetic HIP identifier."""
    data[dest_column] = make_tyc(
        data[src_columns[0]],
        data[src_columns[1]],
        data[src_columns[2]],
    )
    data.remove_columns(src_columns)


def load_tyc_spec() -> Table:
    """Load the TYC2 spectral type catalogue."""
    print("Loading TYC2 spectral types")
    with open_cds_tarfile(VIZIER_DIR / "tyc2spec.tar.gz") as tf:
        data = tf.read_gzip("catalog.dat", ["TYC1", "TYC2", "TYC3", "SpType"])

    parse_tyc_cols(data)
    data.add_index("HIP")
    return data


def _load_ascc_section(tf: TarCds, table: str) -> Table:
    print(f"  Loading {table}")
    section = tf.read_gzip(
        table,
        [
            "Bmag",
            "Vmag",
            "e_Bmag",
            "e_Vmag",
            "d3",
            "TYC1",
            "TYC2",
            "TYC3",
            "Jmag",
            "e_Jmag",
            "Hmag",
            "e_Hmag",
            "Kmag",
            "e_Kmag",
            "SpType",
        ],
    )

    section = section[section["TYC1"] != 0]
    parse_tyc_cols(section)

    convert_cols = [
        "Bmag",
        "Vmag",
        "e_Bmag",
        "e_Vmag",
        "Jmag",
        "e_Jmag",
        "Hmag",
        "e_Hmag",
        "Kmag",
        "e_Kmag",
    ]
    for col in convert_cols:
        section[col] = section[col].astype(np.float64)
        section[col].convert_unit_to(u.mag)
        section[col].format = ".3f"

    return section


def load_ascc() -> Table:
    """Load ASCC from VizieR archive."""

    print("Loading ASCC")
    with open_cds_tarfile(VIZIER_DIR / "ascc.tar.gz") as tf:
        data = None
        for data_file in tf.tf:
            path = PurePath(data_file.name)
            if path.parent != PurePath(".") or not path.stem.startswith("cc"):
                continue
            section_data = _load_ascc_section(tf, path.stem)
            if data is None:
                data = section_data
            else:
                data = vstack([data, section_data], join_type="exact")

    data = unique(data.group_by(["HIP", "d3"]), keys=["HIP"])
    data.rename_column("d3", "Comp")
    data.add_index("HIP")
    return data


def load_tyc2_suppl1() -> Table:
    """Loads the Tycho-2 supplement 1 data."""
    print("Loading TYC2 supplement 1")
    reader = io_ascii.get_reader(
        io_ascii.Cds,
        readme=str(VIZIER_DIR / "tyc2.readme"),
        include_names=["TYC1", "TYC2", "TYC3", "BTmag", "VTmag"],
    )
    reader.data.table_name = "suppl_1.dat"
    with gzip.open(VIZIER_DIR / "tyc2suppl_1.dat.gz") as gzf:
        data = reader.read(gzf)

    parse_tyc_cols(data)
    data = data[np.logical_not(np.logical_and(data["VTmag"].mask, data["BTmag"].mask))]
    # Magnitude transformation formulae from Section 2.2 "Contents of the Tycho Catalogue"
    data["BT-VT"] = data["BTmag"] - data["VTmag"]
    data["Vmag"] = data["VTmag"] - 0.090 * data["BT-VT"]
    data["Bmag"] = data["Vmag"] + 0.850 * data["BT-VT"]
    data.remove_columns(["BTmag", "VTmag", "BT-VT"])
    return data


def load_tyc_hd() -> Table:
    """Load the Tycho-HD cross index."""
    print("Loading TYC-HD cross index")
    with open_cds_tarfile(VIZIER_DIR / "tyc2hd.tar.gz") as tf:
        data = tf.read_gzip("tyc2_hd.dat", ["TYC1", "TYC2", "TYC3", "HD"])

    parse_tyc_cols(data)

    err_del = np.array(TYC_HD_ERRATA["delete"] + [a[1] for a in TYC_HD_ERRATA["add"]])
    data = data[np.logical_not(np.isin(data["HD"], err_del))]

    err_add = Table(
        np.array(TYC_HD_ERRATA["add"]), names=["HIP", "HD"], dtype=[np.int64, np.int64]
    )

    data = vstack([data, err_add], join_type="exact")

    data = unique(data.group_by("HD"), keys="HIP")
    data = unique(data.group_by("HIP"), keys="HD")

    return data


class TYCTeffReader(WorkaroundCDSReader):
    """Custom CDS loader for the TYC Teff table to reduce memory usage."""

    def __init__(self, readme: IO):
        super().__init__("tycall.dat", ["Tycho", "Teff"], [], readme)

    def create_table(self) -> Table:
        """Creates the table."""
        return Table(
            [
                np.empty(self.record_count, np.int64),
                np.empty(self.record_count, np.float64),
            ],
            names=["TYC", "teff_val"],
        )

    def process_line(self, table: Table, record: int, fields: dict[str, str]) -> bool:
        """Processes fields from a line of the input file."""
        try:
            tycsplit = fields["Tycho"].split("-")
            tyc = (
                int(tycsplit[0])
                + int(tycsplit[1]) * 10000
                + int(tycsplit[2]) * 1000000000
            )
            teff = float(fields["Teff"])
        except ValueError:
            tyc = 0
            teff = 99999

        if teff != 99999:
            table["TYC"][record] = tyc
            table["teff_val"][record] = teff
            return True

        return False


def load_tyc_teff() -> Table:
    """Load the Tycho-2 effective temperatures."""
    print("Loading TYC2 effective temperatures")
    with tarfile.open(VIZIER_DIR / "tyc2teff.tar.gz", "r:gz") as tf:
        with tf.extractfile("./ReadMe") as readme:
            reader = TYCTeffReader(readme)

        with tf.extractfile("./tycall.dat.gz") as gzf, gzip.open(
            gzf, "rt", encoding="ascii"
        ) as f:
            data = reader.read(f)

        data["teff_val"].unit = u.K
        data.rename_column("TYC", "HIP")
        return unique(data, keys=["HIP"])


def load_sao() -> Table:
    """Load the SAO-TYC2 cross match."""
    print("Loading SAO-TYC2 cross match")
    xmatch_files = [
        "sao_tyc2_xmatch.csv",
        "sao_tyc2_suppl1_xmatch.csv",
        "sao_tyc2_suppl2_xmatch.csv",
    ]
    data = vstack(
        [
            io_ascii.read(
                XMATCH_DIR / f,
                include_names=["SAO", "TYC1", "TYC2", "TYC3", "angDist", "delFlag"],
                format="csv",
                converters={"delFlag": [io_ascii.convert_numpy(np.str)]},
            )
            for f in xmatch_files
        ],
        join_type="exact",
    )

    data = data[data["delFlag"].mask]
    data.remove_column("delFlag")

    parse_tyc_cols(data)

    data = unique(data.group_by(["HIP", "angDist"]), keys=["HIP"])
    data.remove_column("angDist")

    return data


def process_tyc(data: Table) -> Table:
    """Processes the TYC data."""
    data = join(
        data,
        load_tyc_spec(),
        keys=["HIP"],
        join_type="left",
        table_names=("gaia", "tspec"),
        metadata_conflicts="silent",
    )

    data = join(
        data,
        load_ascc(),
        keys=["HIP"],
        join_type="left",
        table_names=("gaia", "ascc"),
        metadata_conflicts="silent",
    )

    data["SpType"] = MaskedColumn(
        data["SpType_gaia"].filled(
            data["SpType_tspec"].filled(data["SpType"].filled(""))
        )
    )
    data["SpType"].mask = data["SpType"] == ""
    data.remove_columns(["SpType_gaia", "SpType_tspec"])

    for base_col in (c[:-5] for c in data.colnames if c.endswith("_gaia")):
        gaia_col = base_col + "_gaia"
        ascc_col = base_col + "_ascc"
        data.rename_column(gaia_col, base_col)
        mask = np.logical_and(data[base_col].mask, data[ascc_col].mask)
        data[base_col] = MaskedColumn(data[base_col].filled(data[ascc_col]), mask=mask)
        data.remove_column(ascc_col)

    data = join(
        data,
        load_tyc2_suppl1(),
        keys=["HIP"],
        join_type="left",
        table_names=["gaia", "tyc2s1"],
        metadata_conflicts="silent",
    )

    data["Vmag"] = MaskedColumn(
        data["Vmag_gaia"].filled(data["Vmag_tyc2s1"]),
        mask=np.logical_and(data["Vmag_gaia"].mask, data["Vmag_tyc2s1"].mask),
    )

    data["Bmag"] = MaskedColumn(
        data["Bmag_gaia"].filled(data["Bmag_tyc2s1"]),
        mask=np.logical_and(data["Bmag_gaia"].mask, data["Bmag_tyc2s1"].mask),
    )

    data.remove_columns(["Vmag_gaia", "Vmag_tyc2s1", "Bmag_gaia", "Bmag_tyc2s1"])

    data = join(
        data,
        load_tyc_hd(),
        keys=["HIP"],
        join_type="left",
        table_names=["gaia", "hd"],
        metadata_conflicts="silent",
    )
    data["HD"] = MaskedColumn(data["HD_hd"].filled(data["HD_gaia"].filled(0)))
    data["HD"].mask = data["HD"] == 0
    data.remove_columns(["HD_gaia", "HD_hd"])

    data = join(data, load_tyc_teff(), keys=["HIP"], join_type="left")

    data = join(
        data, load_sao(), keys=["HIP"], table_names=["gaia", "sao"], join_type="left"
    )
    data["SAO"] = MaskedColumn(data["SAO_sao"].filled(data["SAO_gaia"].filled(0)))
    data["SAO"].mask = data["SAO"] == 0
    data.remove_columns(["SAO_gaia", "SAO_sao"])

    return data
