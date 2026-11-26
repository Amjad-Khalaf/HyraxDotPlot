from dataclasses import replace
from pathlib import Path

import pytest

from hyraxdotplot.assembly import Assembly
from hyraxdotplot.filtering import AlignmentFilter, filter_alignments
from hyraxdotplot.parsers.nucmer import parse_nucmer
from hyraxdotplot.parsers.paf import parse_paf
from hyraxdotplot.projection import (
    ProjectedAlignment,
    project_alignment,
    project_alignments,
)
from hyraxdotplot.validation import validate_alignments

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def query_assembly():
    return Assembly.from_fai(FIXTURE_DIR / "example.fai")


@pytest.fixture
def target_assembly():
    return Assembly.from_fai(FIXTURE_DIR / "target.fai")


@pytest.fixture
def example_alignments():
    return parse_paf(FIXTURE_DIR / "example.paf")


def test_forward_projection(example_alignments, query_assembly, target_assembly):
    result = project_alignment(example_alignments[0], query_assembly, target_assembly)

    assert isinstance(result, ProjectedAlignment)
    assert result.x_coordinates == (200, 600)
    assert result.y_coordinates == (400, 800)


def test_reverse_projection(example_alignments, query_assembly, target_assembly):
    result = project_alignment(example_alignments[1], query_assembly, target_assembly)

    assert result.x_coordinates == (1100, 1400)
    assert result.y_coordinates == (2300, 2000)


def test_projection_with_query_offset(
    example_alignments, query_assembly, target_assembly
):
    result = project_alignment(example_alignments[2], query_assembly, target_assembly)

    assert result.x_coordinates == (1600, 1800)
    assert result.y_coordinates == (100, 300)


@pytest.mark.parametrize(
    "strand,expected_y",
    [
        ("+", (400, 900)),
        ("-", (900, 400)),
    ],
)
def test_unequal_alignment_spans(
    example_alignments, query_assembly, target_assembly, strand, expected_y
):
    alignment = replace(example_alignments[0], target_end=900, strand=strand)

    result = project_alignment(alignment, query_assembly, target_assembly)

    assert result.x_coordinates == (200, 600)
    assert result.y_coordinates == expected_y

    assert abs(result.x_end - result.x_start) == 400
    assert abs(result.y_end - result.y_start) == 500


def test_projection_preserves_alignment(
    example_alignments, query_assembly, target_assembly
):
    alignment = example_alignments[1]

    result = project_alignment(alignment, query_assembly, target_assembly)

    assert result.alignment is alignment
    assert result.alignment.query_start == 100
    assert result.alignment.query_end == 400
    assert result.alignment.target_start == 500
    assert result.alignment.target_end == 800
    assert result.alignment.identity == pytest.approx(95)
    assert result.alignment.strand == "-"


def test_projection_at_assembly_boundaries(
    example_alignments, query_assembly, target_assembly
):
    alignment = replace(
        example_alignments[0],
        query="contig_3",
        target="contig_6",
        query_start=0,
        query_end=400,
        target_start=0,
        target_end=1000,
    )

    result = project_alignment(alignment, query_assembly, target_assembly)

    assert result.x_coordinates == (1600, 2000)
    assert result.y_coordinates == (1500, 2500)


def test_single_base_projections(example_alignments, query_assembly, target_assembly):
    alignment = replace(
        example_alignments[0], query_start=0, query_end=1, target_start=0, target_end=1
    )

    result = project_alignment(alignment, query_assembly, target_assembly)

    assert result.x_coordinates == (0, 1)
    assert result.y_coordinates == (0, 1)


@pytest.mark.parametrize(
    "changes,error_message",
    [
        ({"query": "unknown"}, "Query contig"),
        ({"target": "unknown"}, "Target contig"),
        ({"query_end": 1001}, "Query alignment end"),
        ({"target_end": 1501}, "Target alignment end"),
    ],
)
def test_invalid_projection(
    example_alignments, query_assembly, target_assembly, changes, error_message
):
    alignment = replace(example_alignments[0], **changes)

    with pytest.raises(ValueError, match=error_message):
        project_alignment(alignment, query_assembly, target_assembly)


def test_project_multiple_alignments(
    example_alignments, query_assembly, target_assembly
):
    results = project_alignments(example_alignments, query_assembly, target_assembly)

    assert len(results) == 3

    assert [item.alignment for item in results] == (example_alignments)

    assert [item.x_coordinates for item in results] == [
        (200, 600),
        (1100, 1400),
        (1600, 1800),
    ]

    assert [item.y_coordinates for item in results] == [
        (400, 800),
        (2300, 2000),
        (100, 300),
    ]


def test_project_empty_alignments(query_assembly, target_assembly):
    assert project_alignments([], query_assembly, target_assembly) == []


@pytest.mark.parametrize(
    "parser,filename",
    [
        (parse_paf, "example.paf"),
        (parse_nucmer, "example.coords"),
    ],
)
def test_complete_projection_pipeline(
    query_assembly, target_assembly, parser, filename
):
    alignments = parser(FIXTURE_DIR / filename)

    validate_alignments(alignments, query_assembly, target_assembly)

    filtered = filter_alignments(alignments, AlignmentFilter(min_identity=95))

    results = project_alignments(filtered, query_assembly, target_assembly)

    assert len(results) == 2
    assert results[0].x_coordinates == (200, 600)
    assert results[0].y_coordinates == (400, 800)

    assert results[1].x_coordinates == (1100, 1400)
    assert results[1].y_coordinates == (2300, 2000)
