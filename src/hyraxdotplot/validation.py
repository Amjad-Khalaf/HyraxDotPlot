from hyraxdotplot.assembly import Assembly
from hyraxdotplot.models import Alignment

"""Up to this point, all internal checks and unit tests confirm the integrity of the generated objects. The point of this part is to check whether the Alignment and Assembly objects are coherent"""


def validate_alignment(
    alignment: Alignment, query_assembly: Assembly, target_assembly: Assembly
) -> None:
    "Validate an alignment against its corresponding assemblies"

    if alignment.query not in query_assembly.contigs:
        raise ValueError(
            f"Query contig '{alignment.query}' not found in query assembly. Something has gone wrong..."
        )

    if alignment.target not in target_assembly.contigs:
        raise ValueError(
            f"Target contig '{alignment.target}' not found in target assembly. Something has gone wrong..."
        )

    query_contig = query_assembly.contigs[alignment.query]
    target_contig = target_assembly.contigs[alignment.target]

    if alignment.query_end > query_contig.length:
        raise ValueError(
            f"Query alignment end ({alignment.query_end}) exceeds contig '{alignment.query}' length ({query_contig.length}). Something has gone wrong..."
        )

    if alignment.target_end > target_contig.length:
        raise ValueError(
            f"Target alignment end ({alignment.target_end}) exceeds contig '{alignment.target}' length ({target_contig.length}). Something has gone wrong..."
        )


def validate_alignments(
    alignments: list[Alignment], query_assembly: Assembly, target_assembly: Assembly
) -> None:
    "Validate a collection of alignments against their two corresponding assemblies"

    for index, alignment in enumerate(alignments, start=1):
        try:
            validate_alignment(alignment, query_assembly, target_assembly)
        except ValueError as exc:
            raise ValueError(f"Alignment {index}: {exc}") from exc
