from bokeh.models import (
    ColorBar,
    ColumnDataSource,
    FixedTicker,
    HoverTool,
    LinearAxis,
    LinearColorMapper,
    Range1d,
    Span,
)
from bokeh.plotting import figure
from bokeh.transform import transform
from matplotlib.colors import LinearSegmentedColormap, to_hex

from hyraxdotplot.assembly import Assembly
from hyraxdotplot.plotting.config import DotPlotConfig
from hyraxdotplot.projection import ProjectedAlignment


def create_alignment_source(alignments: list[ProjectedAlignment]) -> ColumnDataSource:
    """Convert projected alignments into Bokeh data source"""

    data = {
        "x": [],
        "y": [],
        "query": [],
        "target": [],
        "query_start": [],
        "query_end": [],
        "target_start": [],
        "target_end": [],
        "query_span": [],
        "target_span": [],
        "identity": [],
        "strand": [],
    }

    for projected in alignments:
        alignment = projected.alignment

        data["x"].append(list(projected.x_coordinates))
        data["y"].append(list(projected.y_coordinates))

        data["query"].append(alignment.query)
        data["target"].append(alignment.target)

        data["query_start"].append(alignment.query_start)
        data["query_end"].append(alignment.query_end)

        data["target_start"].append(alignment.target_start)
        data["target_end"].append(alignment.target_end)

        data["query_span"].append(alignment.query_span)
        data["target_span"].append(alignment.target_span)

        data["identity"].append(alignment.identity)
        data["strand"].append(alignment.strand)

    return ColumnDataSource(data=data)


def create_identity_mapper(config: DotPlotConfig) -> LinearColorMapper:
    """Create a continuous colour mapper for alignment identities"""

    colormap = LinearSegmentedColormap.from_list(
        "hyrax_identity", ["orange", "red", "black"], N=256
    )

    palette = [to_hex(colormap(index / 255)) for index in range(256)]

    return LinearColorMapper(
        palette=palette, low=config.identity_min, high=config.identity_max
    )


def add_contig_boundaries(
    plot: figure,
    query_assembly: Assembly,
    target_assembly: Assembly,
    config: DotPlotConfig,
) -> None:
    """Add vertical and horizontal contig boundaries"""

    if not config.show_contig_boundaries:
        return

    for contig in query_assembly.contigs.values():
        if contig.offset == 0:
            continue

        boundary = Span(
            location=contig.offset,
            dimension="height",
            line_color=config.boundary_colour,
            line_width=config.boundary_width,
        )

        plot.add_layout(boundary)

    for contig in target_assembly.contigs.values():
        if contig.offset == 0:
            continue
    
        boundary = Span(
            location=contig.offset,
            dimension="width",
            line_color=config.boundary_colour,
            line_width=config.boundary_width,
        )
    
        plot.add_layout(boundary)

def add_contig_labels(
    plot: figure,
    query_assembly: Assembly,
    target_assembly: Assembly,
    config: DotPlotConfig,
) -> None:
    """Add contig names at their cumulative midpoints"""
    if not config.show_contig_labels:
        return

    query_labels = {
        contig.offset + contig.length / 2: contig.name
        for contig in query_assembly.contigs.values()
    }

    target_labels = {
        contig.offset + contig.length / 2: contig.name
        for contig in target_assembly.contigs.values()
    }

    x_axis = LinearAxis(
        ticker=FixedTicker(ticks=list(query_labels.keys())),
        major_label_overrides=query_labels,
        major_label_orientation=0.8,
        major_label_text_font_size="9pt",
        major_tick_line_color=None,
        minor_tick_line_color=None,
    )

    y_axis = LinearAxis(
        ticker=FixedTicker(ticks=list(query_labels.keys())),
        major_label_overrides=target_labels,
        major_label_orientation=0.8,
        major_label_text_font_size="9pt",
        major_tick_line_color=None,
        minor_tick_line_color=None,
    )

    plot.add_layout(x_axis, "above")
    plot.add_layout(y_axis, "right")


def create_dotplot(
    alignments: list[ProjectedAlignment],
    query_assembly: Assembly,
    target_assembly: Assembly,
    config: DotPlotConfig | None = None,
) -> figure:
    """Code to actually create the interactive dotplot!"""

    if config is None:
        config = DotPlotConfig()

    if query_assembly.total_length <= 0 or target_assembly.total_length <= 0:
        raise ValueError("Assemblies must not be empty, something is wrong...")

    source = create_alignment_source(alignments)

    mapper = create_identity_mapper(config)

    plot = figure(
        title=config.title,
        width=config.width,
        height=config.height,
        x_axis_label=config.x_axis_label,
        y_axis_label=config.y_axis_label,
        x_range=Range1d(start=0, end=query_assembly.total_length),
        y_range=Range1d(start=0, end=target_assembly.total_length),
        tools="pan,wheel_zoom,box_zoom,reset,save",
        active_scroll="wheel_zoom",
    )

    plot.xgrid.visible = False
    plot.ygrid.visible = False

    add_contig_boundaries(plot, query_assembly, target_assembly, config)

    add_contig_labels(plot, query_assembly, target_assembly, config)

    renderer = plot.multi_line(
        xs="x",
        ys="y",
        source=source,
        line_color=transform("identity", mapper),
        line_width=config.line_width,
        name="alignments",
    )

    hover = HoverTool(
        renderers=[renderer],
        tooltips=[
            ("Query (X)", "@query"),
            ("Target (Y)", "@target"),
            ("Query interval", "@query_start{0,0}-@query_end{0,0}"),
            ("Target interval", "@target_start{0,0}-@target_end{0,0}"),
            ("Query span", "@query_span{0,0} bp"),
            ("Target span", "@target_span{0,0} bp"),
            ("Identity", "@identity{0.00}%"),
            ("Strand", "@strand"),
        ],
    )

    plot.add_tools(hover)

    color_bar = ColorBar(color_mapper=mapper, title="Identity (%)")

    plot.add_layout(color_bar, "right")

    return plot
