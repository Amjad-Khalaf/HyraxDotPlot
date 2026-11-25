import pytest

from hyraxdotplot.models import Contig


def test_contig_creation():
    contig = Contig("contig_1", 1000, 500)
    assert contig.name == "contig_1"
    assert contig.length == 1000
    assert contig.offset == 500


def test_contig_default_offset():
    contig = Contig("contig_1", 1000)
    assert contig.offset == 0


def test_contig_end():
    contig = Contig("contig_1", 1000, 500)
    assert contig.end == 1500


def test_invalid_contig_length_1():
    with pytest.raises(ValueError):
        Contig("contig_1", 0)


def test_invalid_contig_length_2():
    with pytest.raises(ValueError):
        Contig("contig_1", -100)


def test_invalid_contig_name():
    with pytest.raises(ValueError):
        Contig(" ", 1000)


def test_invalid_contig_offset():
    with pytest.raises(ValueError):
        Contig("contig_1", 1000, -500)
