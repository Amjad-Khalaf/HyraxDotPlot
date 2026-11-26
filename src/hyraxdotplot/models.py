from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class Contig:
    """Represents a contig sequence in an assembly"""

    name: str
    length: int
    offset: int = 0

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("Contig name cannot be empty")
        if self.length <= 0:
            raise ValueError("Contig length must be a positive integer")
        if self.offset < 0:
            raise ValueError("Contig offset cannot be negative")

    @property
    def end(self) -> int:
        """Return the cumulative end coordinate"""
        return self.offset + self.length


@dataclass(frozen=True)
class Alignment:
    """
    Represents a pairwise alignment between two contigs

    Coordinates are zero-based and half-open,
    strand represents the orientation of the query sequence relative to the target sequence
    """

    query: str
    target: str

    query_start: int
    query_end: int

    target_start: int
    target_end: int

    identity: float
    strand: Literal["+", "-"]

    def __post_init__(self) -> None:
        if not self.query.strip():
            raise ValueError("Query name cannot be blank!")

        if not self.target.strip():
            raise ValueError("Target name cannot be blank!")

        if not 0 <= self.query_start < self.query_end:
            raise ValueError(
                "Query coordinates are weird, start is negative or larger than end coordinate!"
            )

        if not 0 <= self.target_start < self.target_end:
            raise ValueError(
                "Target coordinates are weird, start is negative or larger than end coordinate!"
            )

        if not 0 <= self.identity <= 100:
            raise ValueError("Identity must be between 0 and 100!")

        if self.strand not in ("-", "+"):
            raise ValueError("Strand must be '+' or '-'!")

    @property
    def query_span(self) -> int:
        "Return the length of the query alignment"
        return self.query_end - self.query_start

    @property
    def target_span(self) -> int:
        "Return the length of the target alignment"
        return self.target_end - self.target_start

    @property
    def is_reverse(self) -> bool:
        "Return True if alignment is reverse-oriented"
        return self.strand == "-"
