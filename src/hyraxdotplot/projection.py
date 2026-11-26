from dataclasses import dataclass

from hyraxdotplot.assembly import Assembly
from hyraxdotplot.models import Alignment
from hyraxdotplot.validation import validate_alignment


@dataclass(frozen=True)
class ProjectedAlignment:
    """
    Represents an alignment with the global coordinate system.

    The x-axis represents the query assembly, and the y-axis represents the target assembly.

    Coordinates describe the endpoints of a plotted segment.
    """

    alignment: Alignment

    x_start: int
    x_end: int

    y_start: int
    y_end: int

    @property
    def x_coordinates(self) -> tuple[int, int]:
        """Return the x-coordinates of the segment"""
        return self.x_start, self.x_end

    @property
    def y_coordinates(self) -> tuple[int, int]:
        """Return the y-coordinates of the segment"""
        return self.y_start, self.y_end


def project_alignment(
    alignment: Alignment, query_assembly: Assembly, target_assembly: Assembly
):
    """Project local alignment coordinates onto the two assemblies"""

    validate_alignment(alignment, query_assembly, target_assembly)

    x_start = query_assembly.to_global(alignment.query, alignment.query_start)
    x_end = query_assembly.to_global(alignment.query, alignment.query_end)

    y_start = target_assembly.to_global(alignment.target, alignment.target_start)
    y_end = target_assembly.to_global(alignment.target, alignment.target_end)

    if alignment.is_reverse:
        y_start, y_end = y_end, y_start

    return ProjectedAlignment(
        alignment=alignment, x_start=x_start, x_end=x_end, y_start=y_start, y_end=y_end
    )


def project_alignments(
    alignments: list[Alignment], query_assembly: Assembly, target_assembly: Assembly
) -> list[ProjectedAlignment]:
    """Project a collection of alignments onto their two assemblies"""
    return [
        project_alignment(alignment, query_assembly, target_assembly)
        for alignment in alignments
    ]
