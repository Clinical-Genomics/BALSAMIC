#!/usr/bin/env python3
from __future__ import annotations

import gzip
import sys
from pathlib import Path
from typing import Optional, TextIO

import click


def open_maybe_gzip(path: str | Path) -> TextIO:
    """Open a file that may be gzipped; '-' means stdin."""
    p = str(path)
    if p == "-":
        return sys.stdin
    if p.endswith(".gz"):
        return gzip.open(p, "rt")
    return open(p, "r", encoding="utf-8")


def open_output(path: str | Path) -> TextIO:
    """Open output file; '-' means stdout."""
    p = str(path)
    if p == "-":
        return sys.stdout
    return open(p, "w", encoding="utf-8")


def convert_vcf_to_bedgraph(
    vcf_path: str | Path,
    bedgraph_path: str | Path,
    track_name: Optional[str] = None,
) -> int:
    """
    Convert a VCF with a single sample column into a bedGraph of AF.
    Assumes FORMAT contains AD and optionally DP.
    """
    n_written = 0

    with open_maybe_gzip(vcf_path) as fin, open_output(bedgraph_path) as fout:
        if track_name:
            fout.write(f'track type=bedGraph name="{track_name}"\n')

        sample_seen = False

        for line in fin:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                if len(fields) < 10:
                    raise click.ClickException(
                        "Expected exactly one sample column (10+ VCF fields)."
                    )
                sample_seen = True
                continue

            if not sample_seen:
                raise click.ClickException("Malformed VCF: missing #CHROM header.")

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom, pos_str, _id, ref, alt, *_rest = parts[:9]
            sample_field = parts[9]

            # Skip symbolic/no-call
            if alt in (".", "*"):
                continue

            fmt = parts[8]
            fmt_keys = fmt.split(":")
            fmt_vals = sample_field.split(":")
            kv = dict(zip(fmt_keys, fmt_vals))

            ad = kv.get("AD")
            dp = kv.get("DP")

            if not ad or ad == ".":
                continue

            try:
                ad_counts = [int(x) for x in ad.split(",") if x != "."]
            except ValueError:
                continue

            if len(ad_counts) < 2:
                continue

            ref_ad = ad_counts[0]
            alt_sum = sum(ad_counts[1:])

            # DP fallback to sum(AD)
            if dp and dp != ".":
                try:
                    dp_val = int(dp)
                except ValueError:
                    dp_val = ref_ad + alt_sum
            else:
                dp_val = ref_ad + alt_sum

            if dp_val <= 0:
                continue

            af = alt_sum / dp_val

            # Skip trivial 0 or 1 AF
            if af == 0 or af == 1:
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                continue

            span = max(1, len(ref))
            start0 = pos - 1
            end0 = start0 + span

            fout.write(f"{chrom}\t{start0}\t{end0}\t{af:.6f}\n")
            n_written += 1

    return n_written


@click.command(context_settings={"show_default": True})
@click.argument(
    "vcf",
    type=click.Path(exists=True, dir_okay=False, allow_dash=True),
)
@click.argument(
    "bedgraph",
    type=click.Path(dir_okay=False, allow_dash=True),
)
@click.option(
    "--track-name",
    help='Optional bedGraph "track" name header.',
)
def cli(vcf: str, bedgraph: str, track_name: Optional[str]) -> None:
    """
    Convert VCF (with 1 sample) to bedGraph of AF from AD/DP.

    VCF and BEDGRAPH can be "-" for stdin/stdout.
    """
    n = convert_vcf_to_bedgraph(
        vcf_path=vcf,
        bedgraph_path=bedgraph,
        track_name=track_name,
    )
    click.echo(f"Done. Wrote {n} rows → {bedgraph}", err=True)


if __name__ == "__main__":
    cli()
