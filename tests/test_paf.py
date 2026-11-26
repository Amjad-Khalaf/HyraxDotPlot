from pathlib import Path

import pytest

from hyraxdotplot.models import Alignment
from hyraxdotplot.parsers.paf import parse_paf

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def example_alignments():
    return parse_paf(FIXTURE_DIR / "example.paf")


def test_parse_paf_count(example_alignments):
    assert len(example_alignments) == 3


def test_parse_paf_object_type(example_alignments):
    assert all(isinstance(alignment, Alignment) for alignment in example_alignments)


def test_parse_paf_contig_names(example_alignments):
    alignment = example_alignments[0]

    assert alignment.query == "contig_1"
    assert alignment.target == "contig_5"


def test_parse_paf_coordinates(example_alignments):
    alignment = example_alignments[0]

    assert alignment.query_start == 200
    assert alignment.query_end == 600
    assert alignment.target_start == 400
    assert alignment.target_end == 800


def test_parse_paf_forward_strand(example_alignments):
    assert example_alignments[0].strand == "+"


def test_parse_paf_reverse_strand(example_alignments):
    alignment = example_alignments[1]

    assert alignment.strand == "-"
    assert alignment.is_reverse is True
    assert alignment.target_start == 500
    assert alignment.target_end == 800


def test_parse_paf_identity(example_alignments):
    assert example_alignments[0].identity == pytest.approx(99.0)
    assert example_alignments[1].identity == pytest.approx(95.0)
    assert example_alignments[2].identity == pytest.approx(90.0)


def test_parse_paf_empty_file(tmp_path):
    path = tmp_path / "empty.paf"
    path.write_text("")

    assert parse_paf(path) == []


def test_parse_paf_blank_lines(tmp_path):
    path = tmp_path / "blank_lines.paf"
    record = (
        "contig_1\t1000\t200\t600\t+\tcontig_5\t1500\t400\t800\t396\t400\t60\tde:f:0.01"
    )
    path.write_text(f"\n{record}\n\n")

    assert len(parse_paf(path)) == 1


@pytest.fixture
def valid_paf_fields():
    return [
        "contig_1",
        "1000",
        "200",
        "600",
        "+",
        "contig_5",
        "1500",
        "400",
        "800",
        "396",
        "400",
        "60",
        "de:f:0.01",
    ]


@pytest.mark.parametrize(
    "column,value",
    [
        (1, "abc"),
        (1, "0"),
        (2, "-1"),
        (3, "1500"),
        (4, "?"),
        (6, "-1"),
        (8, "2000"),
        (12, "de:f:-0.1"),
        (12, "de:f:nan"),
        (12, "de:f:1.5"),
        (12, "de:f:abc"),
    ],
)
def test_parse_paf_invalid_record(tmp_path, valid_paf_fields, column, value):
    valid_paf_fields[column] = value
    path = tmp_path / "invalid.paf"
    path.write_text("\t".join(valid_paf_fields) + "\n")
    with pytest.raises(ValueError, match="line 1"):
        parse_paf(path)


def test_parse_paf_mising_columns(tmp_path):
    path = tmp_path / "incomplete.paf"
    path.write_text("contig_1\t1000\t200\n")

    with pytest.raises(ValueError, match="12 columns"):
        parse_paf(path)


def test_parse_paf_missing_de_tag(tmp_path, valid_paf_fields):
    path = tmp_path / "missing_de.paf"
    path.write_text("\t".join(valid_paf_fields[:12]) + "\n")

    with pytest.raises(ValueError, match="missing the de"):
        parse_paf(path)


def test_parse_paf_missing_de_after_valid_record(tmp_path, valid_paf_fields):
    path = tmp_path / "missing_de_second_record.paf"

    valid_record = "\t".join(valid_paf_fields)
    invalid_record = "\t".join(valid_paf_fields[:12])

    path.write_text(f"{valid_record}\n{invalid_record}\n")

    with pytest.raises(ValueError, match="line 2"):
        parse_paf(path)


def test_parse_paf_optional_tag_order(tmp_path, valid_paf_fields):
    path = tmp_path / "tags.paf"

    fields = valid_paf_fields[:12] + ["tp:A:P", "NM:i:4", "de:f:0.01"]

    path.write_text("\t".join(fields) + "\n")

    alignments = parse_paf(path)
    assert alignments[0].identity == pytest.approx(99.0)
