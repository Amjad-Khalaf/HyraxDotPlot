from dataclasses import dataclass


@dataclass(frozen=True)
class DotPlotConfig:
    """Configuration for the appearance of the dotplot"""

    title: str = "HyraxDotPlot"
    width: int = 800
    height: int = 600

    x_axis_label: str = "Query assembly (bp)"
    y_axis_label: str = "Target assembly (bp)"

    line_width: float = 1.5

    show_contig_boundaries: bool = True
    show_contig_labels: bool = True

    boundary_colour: str = "grey"
    boundary_width: float = 0.5

    identity_min: float = 50.0
    identity_max: float = 100.0

    def __post_init__(self) -> None:
        if self.width <= 0 or self.height <= 0:
            raise ValueError(
                "Please provide valid plot dimensions, they must be positive integers..."
            )

        if not self.line_width > 0:
            raise ValueError("Line width must be a positive integer...")

        if not self.boundary_width > 0:
            raise ValueError("Boundary width must be a positive integer...")

        if not 0 <= self.identity_min < self.identity_max <= 100:
            raise ValueError("Invalid identity colour range...")
