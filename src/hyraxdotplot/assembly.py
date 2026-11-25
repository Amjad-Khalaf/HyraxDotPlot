from dataclasses import dataclass

from hyraxdotplot.models import Contig


@dataclass
class Assembly:
    """Represent a genome assembly containing ordered contigs"""

    contigs: dict[str, Contig]

    @property
    def total_length(self) -> int:
        """Return the combined length of all contigs"""
        return sum(contig.length for contig in self.contigs.values())

    @property
    def contig_count(self) -> int:
        """Return the number of contigs"""
        return len(self.contigs)

    def to_global(self, contig_name: str, position: int) -> int:
        """Convert a local boundary coordinate to an assembly coordinate"""
        contig = self.contigs[contig_name]
        if not 0 <= position <= contig.length:
            raise ValueError("Position is outside contig boundaries")
        return contig.offset + position
