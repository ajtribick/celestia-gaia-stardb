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

"""Routines for parsing spectral types."""

import re
from enum import IntEnum

from arpeggio import (
    NoMatch,
    OneOrMore,
    Optional,
    ParserPython,
    PTNodeVisitor,
    RegExMatch,
    ZeroOrMore,
    visit_parse_tree,
)


class CelMkClass(IntEnum):
    """Celestia MK and WD classes."""

    O = 0x0000
    B = 0x0100
    A = 0x0200
    F = 0x0300
    G = 0x0400
    K = 0x0500
    M = 0x0600
    R = 0x0700
    S = 0x0800
    N = 0x0900
    WC = 0x0A00
    WN = 0x0B00
    UNKNOWN = 0x0C00
    L = 0x0D00
    T = 0x0E00
    C = 0x0F00
    DA = 0x1000
    DB = 0x1100
    DC = 0x1200
    DO = 0x1300
    DQ = 0x1400
    DZ = 0x1500
    D = 0x1600
    DX = 0x1700


CEL_UNKNOWN_SUBCLASS = 0x00A0


class CelLumClass(IntEnum):
    """Celestia luminosity classes."""

    IA0 = 0x0000
    IA = 0x0001
    IB = 0x0002
    II = 0x0003
    III = 0x0004
    IV = 0x0005
    V = 0x0006
    VI = 0x0007
    UNKNOWN = 0x0008


CEL_UNKNOWN_STAR = CelMkClass.UNKNOWN + CEL_UNKNOWN_SUBCLASS + CelLumClass.UNKNOWN

# pylint: disable=missing-function-docstring,multiple-statements

# Format specification


def spacer():
    return ZeroOrMore(" ")


def rangesym():
    return spacer, ["/", "-"], spacer


def uncertain():
    return Optional([":", "?"])


def numeric():
    return RegExMatch(r"[0-9]+(\.[0-9]*)?")


def roman():
    return RegExMatch(r"(I[VX]|VI{0,3}|I{1,3})([Aa]?[Bb]?)")


def prefix():
    return ["esd", "sd", "d", "g", "c"]


# MS stars, spectra written as, e.g. M3S


def msnorange():
    return "M", spacer, Optional(numeric, spacer), "S"


def msrange():
    return [
        ("M", spacer, numeric, rangesym, numeric, spacer, "S"),
        ("M", spacer, numeric, ["-", "+"], spacer, "S"),
    ]


# normal MK spectra, e.g. K4 plus Wolf-Rayet stars


def mkclass():
    return ["WN", "WC", "WO", "WR", "O", "B", "A", "F", "G", "K", "M", "L", "T", "Y"]


def mknorange():
    return [
        (mkclass, spacer, "(", numeric, ")"),
        (mkclass, Optional(spacer, numeric)),
    ]


def mkrange():
    return [
        (mkclass, spacer, numeric, rangesym, numeric),
        (mkclass, spacer, "(", numeric, rangesym, numeric, ")"),
        (mkclass, spacer, numeric, ["-", "+"]),
    ]


def mkmsnorange():
    return [msnorange, mknorange]


def mkmsrange():
    return [msrange, mkrange]


def mktype():
    return [
        (mkmsnorange, rangesym, mkmsnorange),
        mkmsrange,
        mkmsnorange,
    ]


def lumclass():
    return ["0", "Vz", roman]


def lumrange():
    return ["Ia0", (lumclass, Optional(rangesym, lumclass))]


def lumtype():
    return [
        (lumrange, Optional(uncertain)),
        ("(", lumrange, ")"),
    ]


def noprefixstar():
    return [
        (mktype, uncertain, spacer, Optional(lumtype)),
        ("(", mktype, ")", spacer, Optional(lumtype)),
    ]


def normalstar():
    return (Optional(prefix, uncertain), noprefixstar)


# metallic stars (kA5hF0mA4 etc.)


def metalprefix():
    return ["g", "h", "k", "m", "He"]


def metalsection():
    return metalprefix, noprefixstar


def metalstar():
    return Optional(noprefixstar), OneOrMore(metalsection)


# S stars and carbon stars - these can have an abundance index in addition to
# the subclass, furthermore carbon stars can be written with the subtype
# written as a suffix, e.g. C4,3J instead of C-J4,3


def cclass():
    return ["C-R", "C-N", "C-J", "C-Hd", "C-H"], uncertain


def scclass():
    return ["SC", "S", cclass, "C", "R", "N"]


def scsuffix():
    return ["R", "N", "J", "Hd", "H"]


def scrange():
    return [
        (numeric, "-", Optional(numeric)),
        (numeric, "+"),
        (numeric, Optional(OneOrMore(" "), "-", OneOrMore(" "), numeric)),
    ]


def scindices():
    return spacer, scrange, Optional(["/", ","], scrange)


def scsuffixed():
    return "C", Optional(scindices), spacer, scsuffix


def scnosuffix():
    return scclass, Optional(scindices)


def sctype():
    return [scsuffixed, scnosuffix]


def scstar():
    return Optional(prefix), sctype, spacer, Optional(lumtype)


# white dwarfs


def wdclass():
    return ["DA", "DB", "DC", "DO", "DZ", "DQ", "DX", "D"]


def wdstar():
    return [
        "WD",
        "wd",
        (wdclass, numeric, Optional(rangesym, numeric)),
        (wdclass, Optional(rangesym, wdclass)),
    ]


def starspectrum():
    return [metalstar, normalstar, scstar, wdstar]


def spectrum():
    return [
        starspectrum,
        ("(", starspectrum, ")"),
        ("[", starspectrum, "]"),
    ]


class SpecVisitor(PTNodeVisitor):
    """Parse tree visitor to compute Celestia spectral type."""

    # pylint: disable=unused-argument,redefined-outer-name,too-many-public-methods

    def visit_spacer(self, node, children):
        return None

    def visit_rangesym(self, node, children):
        return None

    def visit_uncertain(self, node, children):
        return None

    def visit_numeric(self, node, children):
        return int(float(node.value))

    def visit_prefix(self, node, children):
        if str(node) == "esd" or str(node) == "sd":
            lclass = CelLumClass.VI
        elif str(node) == "d":
            lclass = CelLumClass.V
        elif str(node) == "g":
            lclass = CelLumClass.III
        elif str(node) == "c":
            lclass = CelLumClass.IB
        else:
            raise ValueError
        return lclass

    def visit_msnorange(self, node, children):
        if len(children.numeric) > 0:
            return "M", children.numeric[0]
        return "M", None

    def visit_msrange(self, node, children):
        return "M", children.numeric[0]

    def visit_mknorange(self, node, children):
        if len(children.numeric) > 0:
            return children.mkclass[0], children.numeric[0]
        return children.mkclass[0], None

    def visit_mkrange(self, node, children):
        return children.mkclass[0], children.numeric[0]

    def visit_mktype(self, node, children):
        mkclass, subclass = children[0]
        if subclass is None:
            subclass = CEL_UNKNOWN_SUBCLASS
        elif subclass < 0:
            subclass = 0x00
        elif subclass > 9:
            subclass = 0x90
        else:
            subclass *= 0x10

        if mkclass == "Y":
            return CelMkClass.T, 0x90
        if mkclass in ("WR", "WO"):
            return CelMkClass.WC, subclass
        return CelMkClass[mkclass], subclass

    def visit_lumrange(self, node, children):
        if len(children) == 2 and (children[0] in ("Ia", "IA")) and children[1] == "0":
            lclass = CelLumClass.IA0
        elif children[0] in ("Ia0", "IA0", "0"):
            lclass = CelLumClass.IA0
        elif children[0].startswith("III"):
            lclass = CelLumClass.III
        elif children[0].startswith("II"):
            lclass = CelLumClass.II
        elif children[0].startswith("IV"):
            lclass = CelLumClass.IV
        elif children[0].startswith("IX"):
            lclass = CelLumClass.VI
        elif children[0] in ("Ia", "IA"):
            lclass = CelLumClass.IA
        elif children[0].startswith("I"):
            lclass = CelLumClass.IB
        elif children[0].startswith("VI"):  # VII, VIII as well
            lclass = CelLumClass.VI
        elif children[0].startswith("V"):
            lclass = CelLumClass.V
        else:
            raise ValueError
        return lclass

    def visit_lumtype(self, node, children):
        return children.lumrange[0]

    def visit_noprefixstar(self, node, children):
        mkclass, mksubclass = children.mktype[0]
        if len(children.lumtype) > 0:
            return mkclass, mksubclass, children.lumtype[0]
        return mkclass, mksubclass, CelLumClass.UNKNOWN

    def visit_normalstar(self, node, children):
        mkclass, mksubclass, lclass = children.noprefixstar[0]

        if len(children.prefix) > 0:
            lclass = children.prefix[0]

        return mkclass, mksubclass, lclass

    def visit_metalsection(self, node, children):
        return children[0], children[1]

    def visit_metalstar(self, node, children):
        sections = dict(children.metalsection)
        if len(children.noprefixstar) > 0:
            sections[" "] = children.noprefixstar[0]
            first_section = sections[" "]
        else:
            first_section = children.metalsection[0][1]

        overall_lclass = children.metalsection[-1][1][2]
        if len(sections) == 1:
            selected = sections[next(iter(sections))]
        elif " " in sections:
            selected = sections[" "]
        elif "h" in sections:
            selected = sections["h"]
        elif "k" in sections:
            selected = sections["k"]
        else:
            selected = first_section

        if selected[2] == CelLumClass.UNKNOWN:
            return selected[0], selected[1], overall_lclass
        return selected

    def visit_cclass(self, node, children):
        return children[0]

    def visit_scrange(self, node, children):
        return int(children[0])

    def visit_scindices(self, node, children):
        return children.scrange[0]

    def visit_scnosuffix(self, node, children):
        if len(children.scindices) > 0:
            return children.scclass[0], children.scindices[0]
        return children.scclass[0], None

    def visit_scsuffixed(self, node, children):
        ctype = "C-" + children.scsuffix[0]
        if len(children.scindices) > 0:
            return ctype, children.scindices[0]
        return ctype, None

    def visit_scstar(self, node, children):
        scclass, scsubclass = children.sctype[0]
        if len(children.prefix) > 0:
            lclass = children.prefix[0]
        elif len(children.lumtype) > 0:
            lclass = children.lumtype[0]
        else:
            lclass = CelLumClass.UNKNOWN

        if scsubclass is None:
            scsubclass = CEL_UNKNOWN_SUBCLASS
        elif scsubclass < 0:
            scsubclass = 0x00
        elif scsubclass > 9:
            scsubclass = 0x90
        else:
            scsubclass *= 0x10

        if scclass in ("C-R", "R"):
            return CelMkClass.R, scsubclass, lclass
        if scclass in ("C-N", "N"):
            return CelMkClass.N, scsubclass, lclass
        if scclass == "SC":
            return CelMkClass.S, scsubclass, lclass
        if scclass.startswith("C"):
            return CelMkClass.C, scsubclass, lclass
        return CelMkClass[scclass], scsubclass, lclass

    def visit_wdstar(self, node, children):
        if children[0] in ("WD", "wd"):
            return CelMkClass.D, CEL_UNKNOWN_SUBCLASS

        try:
            wdclass = CelMkClass[children.wdclass[0]]
        except KeyError:
            wdclass = CelMkClass.D

        if len(children.numeric) > 0:
            wdsubclass = children.numeric[0]
            if wdsubclass < 0:
                wdsubclass = 0
            elif wdsubclass > 9:
                wdsubclass = 0x90
            else:
                wdsubclass *= 0x10
        else:
            wdsubclass = CEL_UNKNOWN_SUBCLASS

        return wdclass, wdsubclass


# pylint: enable=missing-function-docstring

PARSER = ParserPython(spectrum, skipws=False)
VISITOR = SpecVisitor()
MULTISEPARATOR = re.compile(
    r"""
        \+\ *                            # plus sign followed by spaces
        (?:
            \.{2,}                       # ellipsis
            |
            (?:
                \(?
                (?:sd|d|g|c|k|h|m|g|He)? # luminosity or A-star prefix
                [OBAFGKM]|W[DNOCR]|wd    # normal spectrum
            )
            |
            (?:0(?!\.))                  # typo of 0 for O, exclude fractions
        )
    """,
    re.VERBOSE,
)


def parse_spectrum(sptype: str) -> int:
    """Parse a spectral type string into a Celestia spectral type."""

    # resolve ambiguity in grammar: B 0-Ia could be interpreted as (B 0-) Ia or B (0-Ia)
    # resolve in favour of latter
    processed_type = sptype.strip().replace("0-Ia", "Ia-0")

    if not processed_type:
        return CEL_UNKNOWN_STAR

    # remove outer brackets
    while (processed_type[0] == "(" and processed_type[-1] == ")") or (
        processed_type[0] == "[" and processed_type[-1] == "]"
    ):
        processed_type = processed_type[1:-1]

    # remove leading uncertainty indicator
    if processed_type[0] == ":":
        processed_type = processed_type[1:]

    # deal with cases where an O-type spectrum is represented using 0
    if processed_type[0] == "0":
        processed_type = "O" + processed_type[1:]

    # remove nebulae and novae (might otherwise be parsed as N-type)
    if processed_type.casefold().startswith(
        "neb".casefold()
    ) or processed_type.casefold().startswith("nova".casefold()):
        return CEL_UNKNOWN_STAR

    # resolve ambiguity about whether + is an open-ended range or identifies a component

    separator_match = MULTISEPARATOR.search(processed_type)
    if separator_match:
        processed_type = processed_type[: separator_match.span()[0]]

    try:
        parse_tree = PARSER.parse(processed_type)
    except NoMatch:
        return CEL_UNKNOWN_STAR

    return sum(visit_parse_tree(parse_tree, VISITOR))
