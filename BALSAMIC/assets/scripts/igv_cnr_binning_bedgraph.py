#!/usr/bin/env python3
from __future__ import annotations

from typing import TextIO, List, Tuple

import click
import pandas as pd


def bin_denoised_segments(
    infile: TextIO,
    outfile: TextIO,
    bin_size: int,
) -> None:
    """
    Read a denoisedCR-style TSV (contig, start, end, value) and
    average every N bins per chromosome into a bedGraph-like output.

    Parameters
    ----------
    infile : TextIO
        Input handle (open text file).
    outfile : TextIO
        Output handle (open text file).
    bin_size : int
        Number of consecutive bins to aggregate.
    """
    chunks: List[Tuple[str, int, int, float]] = []
    current_chr: str | None = None
    buffer: list[tuple[int, int, float]] = []

    for line in infile:
        # Skip GATK-style header / metadata lines
        if line.startswith("@") or line.startswith("CONTIG"):
            continue

        parts = line.strip().split()
        if len(parts) < 4:
            continue

        chrom, start_str, end_str, value_str = parts[:4]

        try:
            start = int(start_str)
            end = int(end_str)
            log2 = float(value_str)
        except ValueError:
            # Malformed numeric fields; skip
            continue

        if chrom != current_chr:
            # Flush previous chromosome
            if buffer:
                for i in range(0, len(buffer), bin_size):
                    chunk = buffer[i : i + bin_size]
                    starts = [c[0] for c in chunk]
                    ends = [c[1] for c in chunk]
                    vals = [c[2] for c in chunk]
                    chunks.append(
                        (current_chr, starts[0], ends[-1], sum(vals) / len(vals))
                    )
            buffer = []
            current_chr = chrom

        buffer.append((start, end, log2))

    # Flush last chromosome
    if buffer and current_chr is not None:
        for i in range(0, len(buffer), bin_size):
            chunk = buffer[i : i + bin_size]
            starts = [c[0] for c in chunk]
            ends = [c[1] for c in chunk]
            vals = [c[2] for c in chunk]
            chunks.append((current_chr, starts[0], ends[-1], sum(vals) / len(vals)))

    df = pd.DataFrame(chunks, columns=["chrom", "start", "end", "value"])
    df.to_csv(outfile, sep="\t", header=False, index=False)


@click.command(context_settings={"show_default": True})
@click.argument(
    "infile",
    type=click.File("r"),
)
@click.argument(
    "outfile",
    type=click.File("w"),
)
@click.option(
    "--bins-per-window",
    "-b",
    type=int,
    default=5,
    show_default=True,
    help="Number of consecutive bins to average together.",
)
def cli(infile: TextIO, outfile: TextIO, bin_size: int) -> None:
    """
    Bin a denoisedCR-style TSV into larger bedGraph intervals.

    INFILE and OUTFILE can be "-" for stdin/stdout.
    """
    if bin_size <= 0:
        raise click.ClickException("--bin-size must be a positive integer")

    bin_denoised_segments(infile=infile, outfile=outfile, bin_size=bin_size)


if __name__ == "__main__":
    cli()
