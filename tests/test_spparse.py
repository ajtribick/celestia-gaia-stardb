"""Tests for spectral type parsing"""

from celestia_gaia.spparse import parse_spectrum


def test_spparse_normal():
    """Tests for parsing normal star spectra"""
    assert parse_spectrum("K") == 0x05a8
    assert parse_spectrum("F5Ia") == 0x0351
    assert parse_spectrum("G2V") == 0x0426
    assert parse_spectrum("L8") == 0x0d88
    assert parse_spectrum("G(2) V") == 0x0426
    assert parse_spectrum("M3 (III)") == 0x0634
    assert parse_spectrum("[WN4]") == 0x0b48
    assert parse_spectrum("O3 0-Ia") == 0x0030
    assert parse_spectrum("O3- 0-Ia") == 0x0030
    assert parse_spectrum("O3+ 0-Ia") == 0x0030
    assert parse_spectrum("B 0.5IA") == 0x0101
    assert parse_spectrum("F5V Fe-0.5") == 0x0356
    assert parse_spectrum("F6/7IV/Vn") == 0x0365
    assert parse_spectrum("G0.  IA     E,V") == 0x0401
    assert parse_spectrum("dM1.5e") == 0x0616
    assert parse_spectrum("sdO2VIII He6") == 0x0027
    assert parse_spectrum("") == 0x0ca8


def test_spparse_peculiar():
    """Tests for parsing peculiar A-type star spectra"""
    assert parse_spectrum("kA2.5hF2mF2") == 0x0328
    assert parse_spectrum("kA5hF0VmF3") == 0x0306
    assert parse_spectrum("kA3mF3II") == 0x0233
    assert parse_spectrum("kB8HeA0hB8IImA0IbSi") == 0x0183


def test_spparse_white_dwarf():
    """Tests for parsing white dwarf spectra"""
    assert parse_spectrum("D") == 0x16a0
    assert parse_spectrum("wd") == 0x16a0
    assert parse_spectrum("WD") == 0x16a0
    assert parse_spectrum("DC6") == 0x1260
    assert parse_spectrum("DAe") == 0x10a0


def test_spparse_carbon():
    """Tests for parsing carbon star spectra"""
    assert parse_spectrum("C3R") == 0x0738
    assert parse_spectrum("C2,3R") == 0x0728
    assert parse_spectrum("C2,4J") == 0x0f28
    assert parse_spectrum("C3,4Hd III") == 0x0f34
    assert parse_spectrum("S4/3 II") == 0x0843
    assert parse_spectrum("S3+/1- III:") == 0x0834


def test_spparse_composite():
    """Tests for parsing composite spectra"""
    assert parse_spectrum("WN3+O") == 0x0b38
    assert parse_spectrum("B+...") == 0x01a8


def test_spparse_zero_for_o():
    """Tests for parsing 0 as O"""
    assert parse_spectrum("02V") == 0x0026
    assert parse_spectrum("07+08") == 0x0078
    assert parse_spectrum("WN7+0") == 0x0b78
