from pathlib import Path

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


FIXTURE_DIR = Path(__file__).parent / "fixtures"


def test_from_fai():
    assembly = Assembly.from_fai(FIXTURE_DIR / "example.fai")

    assert assembly.contig_count == 3
    assert assembly.total_length == 2000


def test_from_fai_contig_order():
    assembly = Assembly.from_fai(FIXTURE_DIR / "example.fai")

    assert list(assembly.contigs) == [
        "contig_1",
        "contig_2",
        "contig_3",
    ]


def test_from_fai_offsets():
    assembly = Assembly.from_fai(FIXTURE_DIR / "example.fai")

    assert assembly.contigs["contig_1"].offset == 0
    assert assembly.contigs["contig_2"].offset == 1000
    assert assembly.contigs["contig_3"].offset == 1600


def test_from_fai_global_position():
    assembly = Assembly.from_fai(FIXTURE_DIR / "example.fai")

    assert assembly.to_global("contig_2", 200) == 1200


def test_from_fai_empty_file(tmp_path):
    path = tmp_path / "empty.fai"
    path.write_text("")

    with pytest.raises(ValueError, match="empty"):
        Assembly.from_fai(path)


def test_from_fai_duplicate_contigs(tmp_path):
    path = tmp_path / "duplicate.fai"
    path.write_text("contig_1\t1000\t10\t60\t61\ncontig_1\t500\t1040\t60\t61\n")

    with pytest.raises(ValueError, match="Duplicate"):
        Assembly.from_fai(path)


def test_from_fai_invalid_length(tmp_path):
    path = tmp_path / "invalid.fai"
    path.write_text("contig_1\t-100\t10\t60\t61\n")

    with pytest.raises(ValueError):
        Assembly.from_fai(path)


def test_from_fai_malformed_record(tmp_path):
    path = tmp_path / "malformed.fai"
    path.write_text("contig_1\n")

    with pytest.raises(ValueError, match="Invalid FAI"):
        Assembly.from_fai(path)
