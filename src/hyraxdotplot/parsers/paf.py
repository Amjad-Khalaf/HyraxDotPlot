from pathlib import Path

from hyraxdotplot.models import Alignment


def _extract_identity(tags: list[str]) -> float:
    """Calculate gap-compressed identity from the PAF de tag"""

    for tag in tags:
        if tag.startswith("de:f:"):
            divergence = float(tag.split(":", 2)[2])

            if not 0 <= divergence <= 1:
                raise ValueError(
                    "Divergence between two contigs needs to be between 0 and 1, check your PAF!"
                )

            return 100 * (1 - divergence)

    raise ValueError(
        "PAF record is missing the de:f: tag, please check how it was generated..."
    )


def _parse_paf_line(line: str) -> Alignment:
    """Convert a single PAF line into an Alignment object"""

    fields = line.rstrip("\n").split("\t")

    if len(fields) < 12:
        raise ValueError(
            "PAF must be tab-formatted, and should contain at least 12 columns. Please check how it was generated..."
        )

    query_length = int(fields[1])
    query_start = int(fields[2])
    query_end = int(fields[3])

    target_length = int(fields[6])
    target_start = int(fields[7])
    target_end = int(fields[8])

    if query_length <= 0 or target_length <= 0:
        raise ValueError("Contig lengths must be positive integers...")

    if query_end > query_length:
        raise ValueError(
            "Your query alignment end is greater than your query length. Something is wrong...."
        )

    if target_end > target_length:
        raise ValueError(
            "Your target alignment end is greater than your target length. Something is wrong...."
        )

    strand = fields[4]

    if strand not in ("+", "-"):
        raise ValueError(
            "Invalid strand value. Your strand values should be either '+' or '-'."
        )

    identity = _extract_identity(fields[12:])

    return Alignment(
        query=fields[0],
        target=fields[5],
        query_start=query_start,
        query_end=query_end,
        target_start=target_start,
        target_end=target_end,
        identity=identity,
        strand=strand,
    )


def parse_paf(path: str | Path) -> list[Alignment]:
    """Read a PAF file and return a list of Alignment objects"""

    alignments: list[Alignment] = []

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue

            try:
                alignment = _parse_paf_line(line)

            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: {exc}") from exc

            alignments.append(alignment)

    return alignments
