from dataclasses import replace
from pathlib import Path

import pytest

from hyraxdotplot.filtering import AlignmentFilter, filter_alignments
from hyraxdotplot.parsers.paf import parse_paf

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def example_alignments():
    return parse_paf(FIXTURE_DIR / "example.paf")


def test_default_filter():
    config = AlignmentFilter()

    assert config.min_identity == 0
    assert config.max_identity == 100
    assert config.min_query_span == 0
    assert config.max_query_span is None
    assert config.strand is None


@pytest.mark.parametrize(
    "changes",
    [
        {"min_identity": -1},
        {"max_identity": 101},
        {"min_identity": 95, "max_identity": 90},
        {"min_identity": float("nan")},
        {"max_identity": float("nan")},
        {"min_query_span": -1},
        {"max_query_span": -1},
        {"min_query_span": 500, "max_query_span": 100},
        {"min_target_span": -1},
        {"max_target_span": -1},
        {"min_target_span": 500, "max_target_span": 100},
        {"strand": "?"},
    ],
)
def test_invalid_filter_configuration(changes):
    with pytest.raises(ValueError):
        AlignmentFilter(**changes)


def test_contig_filter_normalisation():
    config = AlignmentFilter(query_contigs={"contig_1"}, target_contigs={"contig_5"})

    assert isinstance(config.query_contigs, frozenset)
    assert isinstance(config.target_contigs, frozenset)


@pytest.mark.parametrize(
    "config,expected_queries",
    [
        (AlignmentFilter(), ["contig_1", "contig_2", "contig_3"]),
        (AlignmentFilter(min_identity=95), ["contig_1", "contig_2"]),
        (AlignmentFilter(max_identity=95), ["contig_2", "contig_3"]),
        (AlignmentFilter(min_query_span=300), ["contig_1", "contig_2"]),
        (AlignmentFilter(max_query_span=300), ["contig_2", "contig_3"]),
        (AlignmentFilter(min_target_span=300), ["contig_1", "contig_2"]),
        (AlignmentFilter(max_target_span=300), ["contig_2", "contig_3"]),
        (AlignmentFilter(strand="+"), ["contig_1", "contig_3"]),
        (AlignmentFilter(strand="-"), ["contig_2"]),
        (AlignmentFilter(query_contigs=frozenset({"contig_1"})), ["contig_1"]),
        (
            AlignmentFilter(target_contigs=frozenset({"contig_5"})),
            ["contig_1", "contig_3"],
        ),
        (AlignmentFilter(query_contigs=frozenset()), []),
        (AlignmentFilter(target_contigs=frozenset()), []),
    ],
)
def test_filter_alignments(example_alignments, config, expected_queries):
    filtered = filter_alignments(example_alignments, config)

    assert [alignment.query for alignment in filtered] == expected_queries


def test_combined_filters(example_alignments):
    config = AlignmentFilter(
        min_identity=95,
        min_query_span=350,
        strand="+",
        target_contigs=frozenset({"contig_5"}),
    )

    filtered = filter_alignments(example_alignments, config)

    assert len(filtered) == 1
    assert filtered[0].query == "contig_1"


def test_no_matching_alignments(example_alignments):
    config = AlignmentFilter(min_identity=100)
    assert filter_alignments(example_alignments, config) == []


def test_empty_input():
    assert filter_alignments([], AlignmentFilter()) == []


def test_identity_boundary(example_alignments):
    config = AlignmentFilter(min_identity=95, max_identity=95)

    filtered = filter_alignments(example_alignments, config)

    assert len(filtered) == 1
    assert filtered[0].identity == pytest.approx(95)


def test_span_boundary(example_alignments):
    config = AlignmentFilter(min_query_span=300, max_query_span=300)
    filtered = filter_alignments(example_alignments, config)

    assert len(filtered) == 1
    assert filtered[0].query_span == 300


def test_unequal_query_target_spans(example_alignments):
    alignment = replace(example_alignments[0], target_end=900)

    config = AlignmentFilter(min_target_span=450)
    assert config.matches(alignment)

    config = AlignmentFilter(min_query_span=450)
    assert not config.matches(alignment)


def test_filter_preserves_original_data(example_alignments):
    original = example_alignments.copy()
    filtered = filter_alignments(example_alignments, AlignmentFilter(min_identity=95))

    assert example_alignments == original
    assert len(example_alignments) == 3
    assert len(filtered) == 2
    assert filtered is not example_alignments


def test_default_filter_returns_new_list(example_alignments):
    filtered = filter_alignments(example_alignments, AlignmentFilter())

    assert filtered == example_alignments
    assert filtered is not example_alignments
