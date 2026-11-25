import pytest

from hyraxdotplot.assembly import Assembly
from hyraxdotplot.models import Contig


@pytest.fixture()
def example_assembly():
    return Assembly(
        contigs={
            "contig_1": Contig("contig_1", 1000, 0),
            "contig_2": Contig("contig_2", 600, 1000),
            "contig_3": Contig("contig_3", 400, 1600),
        }
    )


def test_total_length(example_assembly):
    assert example_assembly.total_length == 2000


def test_contig_count(example_assembly):
    assert example_assembly.contig_count == 3


def test_global_position(example_assembly):
    with pytest.raises(ValueError):
        example_assembly.to_global("contig_2", 700)


def test_unknown_contig(example_assembly):
    with pytest.raises(KeyError):
        example_assembly.to_global("unknown", 100)
