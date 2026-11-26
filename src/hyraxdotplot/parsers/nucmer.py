from pathlib import Path

from hyraxdotplot.models import Alignment


def _normalise_coordinates(start: int, end: int) -> tuple[int, int]:
    """Convert 1-based inclusive coordinates to zero-based half-open, this matches the PAF style"""

    if start < 1 or end < 1:
        raise ValueError("NUCMER coordinates must be positive integers!")

    return min(start, end) - 1, max(start, end)


def _parse_nucmer_line(line: str) -> Alignment:
    """Convert a NUCMER coordinate record into an Alignment object"""

    fields = line.rstrip("\r\n").split("\t")

    if len(fields) != 11:
        raise ValueError(
            "NUCMER record must contain 11 columns. Please generate this file with show-coords -T -l"
        )

    raw_target_start = int(fields[0])
    raw_target_end = int(fields[1])

    raw_query_start = int(fields[2])
    raw_query_end = int(fields[3])

    identity = float(fields[6])

    target_length = int(fields[7])
    query_length = int(fields[8])

    if target_length <= 0 or query_length <= 0:
        raise ValueError("Sequence lengths must be positive integers!")

    target_forward = raw_target_start <= raw_target_end
    query_forward = raw_query_start <= raw_query_end

    strand = "+" if target_forward == query_forward else "-"

    target_start, target_end = _normalise_coordinates(raw_target_start, raw_target_end)
    query_start, query_end = _normalise_coordinates(raw_query_start, raw_query_end)

    if target_end > target_length:
        raise ValueError(
            "Target end coordinates exceed target sequence length. Something is wrong..."
        )

    if query_end > query_length:
        raise ValueError(
            "Query end coordinates exceed query sequence length. Something is wrong..."
        )

    return Alignment(
        query=fields[10],
        target=fields[9],
        query_start=query_start,
        query_end=query_end,
        target_start=target_start,
        target_end=target_end,
        identity=identity,
        strand=strand,
    )


def parse_nucmer(path: str | Path) -> list[Alignment]:
    """Read a NUCMER coordinates file, and return a list of Alignment objects"""

    alignments: list[Alignment] = []

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()

            if not stripped:
                continue

            if stripped == "NUCMER":
                continue

            if stripped.startswith("[S1]"):
                continue

            if line_number == 1 and "\t" not in stripped:
                parts = stripped.split()
                if len(parts) == 2:
                    continue

            try:
                alignment = _parse_nucmer_line(line)
            except ValueError as exc:
                raise ValueError(f"{path}: line {line_number}: {exc}") from exc

            alignments.append(alignment)

    return alignments
