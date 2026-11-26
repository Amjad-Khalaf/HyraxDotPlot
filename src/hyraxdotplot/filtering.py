from dataclasses import dataclass
from typing import Literal

from hyraxdotplot.models import Alignment


@dataclass(frozen=True)
class AlignmentFilter:
    """Configuration for filtering pairwise alignments"""

    min_identity: float = 0.0
    max_identity: float = 100.0

    min_query_span: int = 0
    max_query_span: int | None = None

    min_target_span: int = 0
    max_target_span: int | None = None

    strand: Literal["+", "-"] | None = None

    query_contigs: frozenset[str] | None = None
    target_contigs: frozenset[str] | None = None

    def __post_init__(self) -> None:
        if not 0 <= self.min_identity <= self.max_identity <= 100:
            raise ValueError(
                "Identity thresholds must satisfy 0 <= self.min_identity <= self.max_identity <= 100"
            )

        if self.min_query_span < 0:
            raise ValueError("Minimum query span cannot be negative...")

        if (
            self.max_query_span is not None
            and self.max_query_span < self.min_query_span
        ):
            raise ValueError(
                "Maximum query span cannot be smaller than minimum query span..."
            )

        if self.min_target_span < 0:
            raise ValueError("Minimum target span cannot be negative...")

        if (
            self.max_target_span is not None
            and self.max_target_span < self.min_target_span
        ):
            raise ValueError(
                "Maximum target span cannot be smaller than minimum target span..."
            )

        if self.strand not in (None, "+", "-"):
            raise ValueError("Invalid strand filter...")

        if self.query_contigs is not None:
            object.__setattr__(self, "query_contigs", frozenset(self.query_contigs))

        if self.target_contigs is not None:
            object.__setattr__(self, "target_contigs", frozenset(self.target_contigs))

    def matches(self, alignment: Alignment) -> bool:
        """Return whether an alignment satisfies all filters requested"""

        if not (self.min_identity <= alignment.identity <= self.max_identity):
            return False

        if alignment.query_span < self.min_query_span:
            return False

        if (
            self.max_query_span is not None
            and alignment.query_span > self.max_query_span
        ):
            return False

        if alignment.target_span < self.min_target_span:
            return False

        if (
            self.max_target_span is not None
            and alignment.target_span > self.max_target_span
        ):
            return False

        if self.strand is not None and alignment.strand != self.strand:
            return False

        if self.query_contigs is not None and alignment.query not in self.query_contigs:
            return False

        if (
            self.target_contigs is not None
            and alignment.target not in self.target_contigs
        ):
            return False

        return True


def filter_alignments(
    alignments: list[Alignment], config: AlignmentFilter
) -> list[Alignment]:
    """Return alignments satisfying filtering conditions requested"""

    return [alignment for alignment in alignments if config.matches(alignment)]
