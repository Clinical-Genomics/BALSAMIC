#!/usr/bin/env python3
from __future__ import annotations

from contextlib import closing
from typing import Optional, Tuple

import click
import numpy as np
from cyvcf2 import VCF


def variant_to_record(variant) -> Optional[str]:
    """
    Convert a cyvcf2 Variant to a bedGraph line.

    Assumptions:
      - exactly one sample
      - FORMAT/AD exists and has at least ref+alt counts
      - FORMAT/DP exists
    """
    ad = variant.format("AD")
    if ad is None:
        return None

    counts = ad[0]
    alt_sum = int(np.sum(counts[1:]))

    dp = variant.format("DP")
    dp_val = int(np.asarray(dp).reshape(-1)[0])

    # Compute AF safely
    if dp_val <= 0:
        af = 1.0
    elif alt_sum == 0:
        af = 0.0
    else:
        af = round(alt_sum / dp_val, 5)

    start = int(variant.POS)
    end = start + 1

    return f"{variant.CHROM}\t{start}\t{end}\t{af}\n"


def convert_vcf_to_bedgraph(
    vcf_path: str,
    bedgraph_path: str,
    track_name: Optional[str] = None,
) -> Tuple[int, int]:
    """
    Convert a VCF into a bedGraph of AF computed from AD/DP.

    Returns:
        (n_written, n_skipped)
    """
    n_written = 0
    n_skipped = 0

    with closing(VCF(str(vcf_path))) as vcf, open(
        bedgraph_path, "w", encoding="utf-8"
    ) as fout:
        if track_name:
            fout.write(f'track type=bedGraph name="{track_name}"\n')

        for variant in vcf:
            record = variant_to_record(variant)
            if record is None:
                n_skipped += 1
                continue
            fout.write(record)
            n_written += 1

    return n_written, n_skipped


@click.command(context_settings={"show_default": True})
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False))
@click.argument("bedgraph", type=click.Path(dir_okay=False))
@click.option("--track-name", help='Optional bedGraph "track" name header.')
def cli(vcf: str, bedgraph: str, track_name: Optional[str]) -> None:
    """Convert VCF (1 sample) to bedGraph of AF from AD/DP."""
    n_written, n_skipped = convert_vcf_to_bedgraph(vcf, bedgraph, track_name=track_name)
    click.echo(
        f"Done. Wrote {n_written} rows → {bedgraph} " f"({n_skipped} variants skipped)",
        err=True,
    )


if __name__ == "__main__":
    cli()
