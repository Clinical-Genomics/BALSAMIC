#!/usr/bin/env python3
from __future__ import annotations

from typing import TextIO, List, Tuple

import click
import pandas as pd


def _parse_denoisedcr_row(line: str) -> tuple[str, int, int, float] | None:
    """Parse a denoisedCR-style row into (chrom, start, end, value), or None if invalid."""
    if line.startswith("@") or line.startswith("CONTIG"):
        return None

    parts = line.strip().split()
    if len(parts) < 4:
        return None

    chrom, start_str, end_str, value_str = parts[:4]
    try:
        return chrom, int(start_str), int(end_str), float(value_str)
    except ValueError:
        return None


def _flush_buffer(
    out: list[tuple[str, int, int, float]],
    chrom: str | None,
    buffer: list[tuple[int, int, float]],
    bin_size: int,
) -> None:
    """Aggregate buffered bins for a chromosome into bin_size chunks and append to out."""
    if not buffer or chrom is None:
        return

    for i in range(0, len(buffer), bin_size):
        chunk = buffer[i : i + bin_size]
        start = chunk[0][0]
        end = chunk[-1][1]
        mean_val = sum(v for _, _, v in chunk) / len(chunk)
        out.append((chrom, start, end, mean_val))


def bin_denoised_segments(
    infile: TextIO,
    outfile: TextIO,
    bin_size: int,
    track_name: str | None = None,
) -> None:
    """
    Read a denoisedCR-style TSV (contig, start, end, value) and
    average every N bins per chromosome into a bedGraph-like output.

    If track_name is provided, write a UCSC bedGraph track header line:
        track type=bedGraph name="..."
    """
    chunks: List[Tuple[str, int, int, float]] = []
    current_chr: str | None = None
    buffer: list[tuple[int, int, float]] = []

    for line in infile:
        parsed = _parse_denoisedcr_row(line)
        if parsed is None:
            continue

        chrom, start, end, log2 = parsed

        if chrom != current_chr:
            _flush_buffer(chunks, current_chr, buffer, bin_size)
            buffer = []
            current_chr = chrom

        buffer.append((start, end, log2))

    _flush_buffer(chunks, current_chr, buffer, bin_size)

    # Optional UCSC track header
    if track_name:
        outfile.write(f'track type=bedGraph name="{track_name}"\n')

    pd.DataFrame(chunks, columns=["chrom", "start", "end", "value"]).to_csv(
        outfile, sep="\t", header=False, index=False
    )


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
@click.option(
    "--track-name",
    type=str,
    required=False,
    help='Optional bedGraph "track" name header.',
)
def cli(
    infile: TextIO, outfile: TextIO, bins_per_window: int, track_name: str | None
) -> None:
    """
    Bin a denoisedCR-style TSV into larger bedGraph intervals.

    INFILE and OUTFILE can be "-" for stdin/stdout.
    """
    if bins_per_window <= 0:
        raise click.ClickException("--bins-per-window must be a positive integer")

    bin_denoised_segments(
        infile=infile,
        outfile=outfile,
        bin_size=bins_per_window,
        track_name=track_name,
    )


if __name__ == "__main__":
    cli()
