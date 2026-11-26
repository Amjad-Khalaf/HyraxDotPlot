from dataclasses import replace
from pathlib import Path

import pytest

from hyraxdotplot.assembly import Assembly
from hyraxdotplot.models import Alignment
from hyraxdotplot.parsers.nucmer import parse_nucmer
from hyraxdotplot.parsers.paf import parse_paf
from hyraxdotplot.validation import validate_alignment, validate_alignments

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def query_assembly():
    return Assembly.from_fai(FIXTURE_DIR / "example.fai")


@pytest.fixture
def target_assembly():
    return Assembly.from_fai(FIXTURE_DIR / "target.fai")


@pytest.fixture
def example_alignment():
    return Alignment(
        query="contig_1",
        target="contig_5",
        query_start=200,
        query_end=600,
        target_start=400,
        target_end=800,
        identity=99.0,
        strand="+",
    )


def test_valid_alignment(example_alignment, query_assembly, target_assembly):
    assert (
        validate_alignment(example_alignment, query_assembly, target_assembly) is None
    )


def test_alignment_at_contig_boundary(
    example_alignment, query_assembly, target_assembly
):
    alignment = replace(
        example_alignment,
        query_start=0,
        query_end=1000,
        target_start=0,
        target_end=1500,
    )

    validate_alignment(alignment, query_assembly, target_assembly)


def test_valid_reverse_alignment(example_alignment, query_assembly, target_assembly):
    alignment = replace(example_alignment, strand="-")

    validate_alignment(alignment, query_assembly, target_assembly)


@pytest.mark.parametrize(
    "changes,error_message",
    [
        (
            {"query": "unknown_contig"},
            "Query contig",
        ),
        (
            {"target": "unknown_contig"},
            "Target contig",
        ),
        (
            {"query_end": 1001},
            "Query alignment end",
        ),
        (
            {"target_end": 1501},
            "Target alignment end",
        ),
    ],
)
def test_invalid_alignment(
    example_alignment, query_assembly, target_assembly, changes, error_message
):
    alignment = replace(example_alignment, **changes)

    with pytest.raises(ValueError, match=error_message):
        validate_alignment(
            alignment,
            query_assembly,
            target_assembly,
        )


def test_validate_multiple_alignments(
    example_alignment, query_assembly, target_assembly
):
    alignments = [example_alignment, replace(example_alignment, strand="-")]

    validate_alignments(alignments, query_assembly, target_assembly)


def test_validate_empty_alignments(query_assembly, target_assembly):
    validate_alignments([], query_assembly, target_assembly)


def test_validate_invalid_second_alignment(
    example_alignment, query_assembly, target_assembly
):
    alignments = [example_alignment, replace(example_alignment, query_end=1001)]

    with pytest.raises(ValueError, match="Alignment 2:"):
        validate_alignments(alignments, query_assembly, target_assembly)


def test_validate_paf_alignments(query_assembly, target_assembly):
    alignments = parse_paf(FIXTURE_DIR / "example.paf")

    validate_alignments(alignments, query_assembly, target_assembly)


def test_validate_nucmer_alignments(query_assembly, target_assembly):
    alignments = parse_nucmer(FIXTURE_DIR / "example.coords")

    validate_alignments(alignments, query_assembly, target_assembly)


def test_swapped_assemblies(query_assembly, target_assembly):
    alignments = parse_paf(FIXTURE_DIR / "example.paf")

    with pytest.raises(ValueError, match="Query contig"):
        validate_alignments(alignments, target_assembly, query_assembly)
