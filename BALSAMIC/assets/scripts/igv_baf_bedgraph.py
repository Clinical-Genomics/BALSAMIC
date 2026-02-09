#!/usr/bin/env python3
from __future__ import annotations

import sys
from contextlib import closing, nullcontext
from pathlib import Path
from typing import ContextManager, Optional, TextIO

import click
import numpy as np
from cyvcf2 import VCF


def open_output(path: str | Path) -> ContextManager[TextIO]:
    """Open output file; '-' means stdout (not closed)."""
    p = str(path)
    if p == "-":
        return nullcontext(sys.stdout)
    return open(p, "w", encoding="utf-8")


def variant_to_record(variant) -> str:
    """
    Convert a cyvcf2 Variant to a BedGraphRecord.

    Assumptions:
      - exactly one sample
      - FORMAT/AD exists and has at least ref+alt counts
      - FORMAT/DP exists and is > 0
    """

    chrom = str(variant.CHROM)

    ad = variant.format("AD")
    dp = variant.format("DP")

    counts = ad[0]
    alt_sum = int(np.sum(counts[1:]))

    dp_val = int(np.asarray(dp).reshape(-1)[0])
    af = round(float(alt_sum / dp_val), 5)

    start = int(variant.POS)
    end = start + 1

    return f"{chrom}\t{start}\t{end}\t{af}\n"


def convert_vcf_to_bedgraph(
    vcf_path: str | Path,
    bedgraph_path: str | Path,
    track_name: Optional[str] = None,
) -> int:
    """Convert a VCF into a bedGraph of AF computed from AD/DP."""
    n_written = 0
    with closing(VCF(str(vcf_path))) as vcf, open_output(bedgraph_path) as fout:
        if track_name:
            fout.write(f'track type=bedGraph name="{track_name}"\n')

        for variant in vcf:
            fout.write(variant_to_record(variant))
            n_written += 1

    return n_written


@click.command(context_settings={"show_default": True})
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False, allow_dash=False))
@click.argument("bedgraph", type=click.Path(dir_okay=False, allow_dash=True))
@click.option("--track-name", help='Optional bedGraph "track" name header.')
def cli(vcf: str, bedgraph: str, track_name: Optional[str]) -> None:
    """Convert VCF (1 sample) to bedGraph of AF from AD/DP."""
    n = convert_vcf_to_bedgraph(vcf, bedgraph, track_name=track_name)
    click.echo(f"Done. Wrote {n} rows → {bedgraph}", err=True)


if __name__ == "__main__":
    cli()
