from pathlib import Path

import pytest

from hyraxdotplot.models import Alignment
from hyraxdotplot.parsers.nucmer import parse_nucmer
from hyraxdotplot.parsers.paf import parse_paf

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def example_alignments():
    return parse_nucmer(FIXTURE_DIR / "example.coords")


def test_parse_nucmer_count(example_alignments):
    assert len(example_alignments) == 3


def test_parse_nucmer_object_type(example_alignments):
    assert all(isinstance(alignment, Alignment) for alignment in example_alignments)


def test_parse_nucmer_contig_names(example_alignments):
    alignment = example_alignments[0]

    assert alignment.query == "contig_1"
    assert alignment.target == "contig_5"


def test_parse_nucmer_coordinates(example_alignments):
    alignment = example_alignments[0]

    assert alignment.query_start == 200
    assert alignment.query_end == 600
    assert alignment.target_start == 400
    assert alignment.target_end == 800


def test_parse_nucmer_forward_strand(example_alignments):
    assert example_alignments[0].strand == "+"


def test_parse_nucmer_reverse_strand(example_alignments):
    alignment = example_alignments[1]

    assert alignment.strand == "-"
    assert alignment.query_start == 100
    assert alignment.query_end == 400
    assert alignment.target_start == 500
    assert alignment.target_end == 800


def test_parse_nucmer_identity(example_alignments):
    assert example_alignments[0].identity == pytest.approx(99.0)
    assert example_alignments[1].identity == pytest.approx(95.0)
    assert example_alignments[2].identity == pytest.approx(90.0)


def test_parse_nucmer_span(example_alignments):
    assert example_alignments[0].query_span == 400
    assert example_alignments[0].target_span == 400


def test_parse_nucmer_first_base(example_alignments):
    assert example_alignments[2].query_start == 0


@pytest.mark.parametrize(
    "raw_start,raw_end,expected",
    [
        (1, 100, (0, 100)),
        (100, 1, (0, 100)),
        (201, 600, (200, 600)),
        (600, 201, (200, 600)),
        (1, 1, (0, 1)),
        (100, 100, (99, 100)),
    ],
)
def test_normalise_coordinates(raw_start, raw_end, expected):
    from hyraxdotplot.parsers.nucmer import _normalise_coordinates

    assert _normalise_coordinates(raw_start, raw_end) == expected


@pytest.fixture
def valid_nucmer_fields():
    return [
        "401",
        "800",
        "201",
        "600",
        "400",
        "400",
        "99.00",
        "1500",
        "1000",
        "contig_5",
        "contig_1",
    ]


@pytest.mark.parametrize(
    "column,value",
    [
        (0, "0"),
        (0, "-1"),
        (1, "2000"),
        (2, "0"),
        (3, "2000"),
        (3, "abc"),
        (6, "-1"),
        (6, "101"),
        (6, "nan"),
        (7, "0"),
        (8, "-1"),
        (9, ""),
        (10, ""),
    ],
)
def test_parse_nucmer_invalid_record(tmp_path, valid_nucmer_fields, column, value):
    valid_nucmer_fields[column] = value

    path = tmp_path / "invalid.coords"
    path.write_text("\t".join(valid_nucmer_fields) + "\n")

    with pytest.raises(ValueError, match="line 1"):
        parse_nucmer(path)


@pytest.mark.parametrize(
    "coordinates,expected_strand",
    [
        (["401", "800", "201", "600"], "+"),
        (["401", "800", "600", "201"], "-"),
        (["800", "401", "201", "600"], "-"),
        (["800", "401", "600", "201"], "+"),
    ],
)
def test_parse_nucmer_orientation(
    tmp_path, valid_nucmer_fields, coordinates, expected_strand
):
    valid_nucmer_fields[:4] = coordinates

    path = tmp_path / "orientation.coords"
    path.write_text("\t".join(valid_nucmer_fields) + "\n")

    alignment = parse_nucmer(path)[0]

    assert alignment.strand == expected_strand
    assert alignment.query_start == 200
    assert alignment.query_end == 600
    assert alignment.target_start == 400
    assert alignment.target_end == 800


def test_parse_nucmer_empty_file(tmp_path):
    path = tmp_path / "empty.coords"
    path.write_text("")

    assert parse_nucmer(path) == []


def test_parse_nucmer_missing_columns(tmp_path):
    path = tmp_path / "incomplete.coords"
    path.write_text("401\t800\t201\n")

    with pytest.raises(ValueError, match="11 columns"):
        parse_nucmer(path)


def test_parse_nucmer_header(tmp_path, valid_nucmer_fields):
    path = tmp_path / "header.coords"

    record = "\t".join(valid_nucmer_fields)

    header = (
        "reference.fasta query.fasta\n"
        "NUCMER\n"
        "\n"
        "[S1]\t[E1]\t[S2]\t[E2]\t[LEN 1]\t[LEN 2]\t"
        "[% IDY]\t[LEN R]\t[LEN Q]\t[TAGS]\n"
    )

    path.write_text(header + record + "\n")

    assert len(parse_nucmer(path)) == 1


def test_paf_nucmer_equivalence():
    paf_alignments = parse_paf(FIXTURE_DIR / "example.paf")

    nucmer_alignments = parse_nucmer(FIXTURE_DIR / "example.coords")

    assert paf_alignments == nucmer_alignments
