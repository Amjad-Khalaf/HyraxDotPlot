from dataclasses import replace

import pytest

from hyraxdotplot.models import Alignment


@pytest.fixture
def example_alignment():
    return Alignment(
        query="contig_1",
        target="contig_2",
        query_start=200,
        query_end=600,
        target_start=400,
        target_end=800,
        identity=99.5,
        strand="+",
    )


def test_alignment_creation(example_alignment):
    assert example_alignment.query == "contig_1"
    assert example_alignment.target == "contig_2"
    assert example_alignment.identity == 99.5
    assert example_alignment.strand == "+"


def test_query_span(example_alignment):
    assert example_alignment.query_span == 400


def test_target_span(example_alignment):
    assert example_alignment.target_span == 400


def test_forward_alignment(example_alignment):
    assert example_alignment.is_reverse is False


def test_reverse_alignment():
    alignment = Alignment(
        query="contig_1",
        target="contig_2",
        query_start=200,
        query_end=600,
        target_start=400,
        target_end=800,
        identity=99.5,
        strand="-",
    )

    assert alignment.query_span == 400
    assert alignment.target_span == 400
    assert alignment.is_reverse is True


def test_unequal_alignment_spans():
    alignment = Alignment(
        query="contig_1",
        target="contig_2",
        query_start=200,
        query_end=600,
        target_start=400,
        target_end=900,
        identity=99.5,
        strand="-",
    )

    assert alignment.query_span == 400
    assert alignment.target_span == 500


@pytest.mark.parametrize(
    "changes",
    [
        {"query": ""},
        {"target": ""},
        {"query_start": -1},
        {"query_end": 100},
        {"query_end": 200},
        {"target_start": -1},
        {"target_end": 300},
        {"target_end": 400},
        {"identity": -1},
        {"identity": 101.9},
        {"identity": float("nan")},
        {"strand": "?"},
    ],
)
def test_invalid_alignment(example_alignment, changes):
    with pytest.raises(ValueError):
        replace(example_alignment, **changes)
