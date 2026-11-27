from pathlib import Path

import pytest
from bokeh.models import (
    ColorBar,
    HoverTool,
    LinearColorMapper,
    Span,
)
from bokeh.plotting import figure

from hyraxdotplot.assembly import Assembly
from hyraxdotplot.parsers.paf import parse_paf
from hyraxdotplot.plotting.config import DotPlotConfig
from hyraxdotplot.plotting.dotplot import (
    create_alignment_source,
    create_dotplot,
    create_identity_mapper,
)
from hyraxdotplot.projection import project_alignments

FIXTURE_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def query_assembly():
    return Assembly.from_fai(FIXTURE_DIR / "example.fai")


@pytest.fixture
def target_assembly():
    return Assembly.from_fai(FIXTURE_DIR / "target.fai")


@pytest.fixture
def projected_alignments(query_assembly, target_assembly):
    alignments = parse_paf(FIXTURE_DIR / "example.paf")
    return project_alignments(alignments, query_assembly, target_assembly)


@pytest.fixture
def example_plot(projected_alignments, query_assembly, target_assembly):
    return create_dotplot(projected_alignments, query_assembly, target_assembly)


def test_default_config():
    config = DotPlotConfig()

    assert config.width == 800
    assert config.height == 600
    assert config.show_contig_boundaries
    assert config.show_contig_labels


@pytest.mark.parametrize(
    "changes",
    [
        {"width": 0},
        {"height": -1},
        {"line_width": 0},
        {"boundary_width": -1},
        {"identity_min": -1},
        {"identity_max": 101},
        {"identity_min": 95, "identity_max": 90},
        {"identity_min": 95, "identity_max": 95},
    ],
)
def test_invalid_config(changes):
    with pytest.raises(ValueError):
        DotPlotConfig(**changes)


def test_alignment_source_column_lengths(projected_alignments):
    source = create_alignment_source(projected_alignments)

    assert all(len(column) == 3 for column in source.data.values())


def test_empty_alignment_source():
    source = create_alignment_source([])
    assert all(len(column) == 0 for column in source.data.values())


def test_identity_mapper():
    config = DotPlotConfig(identity_min=90, identity_max=100)
    mapper = create_identity_mapper(config)

    assert isinstance(mapper, LinearColorMapper)
    assert mapper.low == 90
    assert mapper.high == 100
    assert len(mapper.palette) == 256

    assert mapper.palette[0] == "#ffa500"
    assert mapper.palette[-1] == "#000000"


def test_create_dotplot(example_plot):
    assert isinstance(example_plot, figure)


def test_dotplot_axis_ranges(example_plot):
    assert example_plot.x_range.start == 0
    assert example_plot.x_range.end == 2000

    assert example_plot.y_range.start == 0
    assert example_plot.y_range.end == 2500


def test_dotplot_alignment_renderer(example_plot):
    renderer = example_plot.select_one({"name": "alignments"})

    assert renderer is not None
    assert len(renderer.data_source.data["x"]) == 3


def test_dotplot_hover_tool(example_plot):
    hover = example_plot.select_one(HoverTool)

    assert hover is not None
    assert ("Query (X)", "@query") in hover.tooltips
    assert ("Target (Y)", "@target") in hover.tooltips

    renderer = example_plot.select_one({"name": "alignments"})
    assert renderer in hover.renderers


def test_dotplot_color_bar(example_plot):
    color_bar = example_plot.select_one(ColorBar)

    assert color_bar is not None
    assert color_bar.title == "Identity (%)"


def test_contig_boundaries(example_plot):
    boundaries = list(example_plot.select({"type": Span}))

    vertical = sorted(
        boundary.location for boundary in boundaries if boundary.dimension == "height"
    )
    horizontal = sorted(
        boundary.location for boundary in boundaries if boundary.dimension == "width"
    )
    print(vertical, horizontal)
    assert vertical == [1000, 1600]
    assert horizontal == [1500]


def test_contig_labels(example_plot):
    x_axis = example_plot.above[0]

    assert x_axis.major_label_overrides == {
        500: "contig_1",
        1300: "contig_2",
        1800: "contig_3",
    }

    y_axis = next(
        axis for axis in example_plot.right if hasattr(axis, "major_label_overrides")
    )

    assert y_axis.major_label_overrides == {750: "contig_5", 2000: "contig_6"}


def test_disable_contig_boundaries(projected_alignments, query_assembly, target_assembly):
    plot = create_dotplot(
        projected_alignments,
        query_assembly,
        target_assembly,
        DotPlotConfig(show_contig_boundaries=False),
    )

    assert not list(plot.select({"type": Span}))


def test_disable_contig_labels(projected_alignments, query_assembly, target_assembly):
    plot = create_dotplot(
        projected_alignments,
        query_assembly,
        target_assembly,
        DotPlotConfig(show_contig_labels=False),
    )

    assert len(plot.xaxis) == 1
    assert len(plot.yaxis) == 1


def test_custom_plot_configuration(
    projected_alignments, query_assembly, target_assembly
):
    config = DotPlotConfig(title="Custom", width=1200, height=900, line_width=4)
    plot = create_dotplot(projected_alignments, query_assembly, target_assembly, config)

    assert plot.title.text == "Custom"
    assert plot.width == 1200
    assert plot.height == 900

    renderer = plot.select_one({"name": "alignments"})
    assert renderer.glyph.line_width == 4


def test_empty_dotplot(query_assembly, target_assembly):
    plot = create_dotplot([], query_assembly, target_assembly)
    renderer = plot.select_one({"name": "alignments"})

    assert renderer.data_source.data["x"] == []
    assert plot.x_range.end == 2000
    assert plot.y_range.end == 2500
