from dataclasses import dataclass


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
