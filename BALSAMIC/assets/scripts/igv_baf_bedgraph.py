#!/usr/bin/env python3
from __future__ import annotations

import gzip
import sys
from dataclasses import dataclass
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


@dataclass
class BedGraphRecord:
    chrom: str
    start: int  # 0-based
    end: int  # half-open
    value: float  # AF


def _parse_header_line(line: str) -> None:
    """Validate that the #CHROM header has a sample column."""
    fields = line.strip().split("\t")
    if len(fields) < 10:
        raise click.ClickException(
            "Expected exactly one sample column (10+ VCF fields)."
        )


def _parse_ad_counts(ad: str) -> Optional[tuple[int, int]]:
    """
    Parse AD string into (ref_ad, alt_sum).

    Returns None if AD is missing or malformed.
    """
    if not ad or ad == ".":
        return None
    try:
        counts = [int(x) for x in ad.split(",") if x != "."]
    except ValueError:
        return None
    if len(counts) < 2:
        return None
    ref_ad = counts[0]
    alt_sum = sum(counts[1:])
    return ref_ad, alt_sum


def _compute_dp(ref_ad: int, alt_sum: int, dp_raw: Optional[str]) -> int:
    """
    Compute DP given raw DP string and AD-derived counts.

    Falls back to ref_ad + alt_sum when DP is missing or malformed.
    """
    if dp_raw and dp_raw != ".":
        try:
            return int(dp_raw)
        except ValueError:
            pass
    return ref_ad + alt_sum


def _extract_bedgraph_record(parts: list[str]) -> Optional[BedGraphRecord]:
    """
    Convert a single VCF record (split into parts) into a BedGraphRecord.

    Returns None if the record should be skipped (e.g., missing AD/DP,
    symbolic ALT, AF=0/1, malformed numeric fields).
    """
    if len(parts) < 10:
        return None

    chrom, pos_str, _id, ref, alt, _, _, _, fmt = parts[:9]
    sample_field = parts[9]

    # Skip symbolic / no-call ALT
    if alt in (".", "*"):
        return None

    # Parse FORMAT and sample fields
    fmt_keys = fmt.split(":")
    fmt_vals = sample_field.split(":")
    kv = dict(zip(fmt_keys, fmt_vals))

    ad_dp = _parse_ad_counts(kv.get("AD", ""))
    if ad_dp is None:
        return None
    ref_ad, alt_sum = ad_dp

    dp_val = _compute_dp(ref_ad, alt_sum, kv.get("DP"))
    if dp_val <= 0:
        return None

    af = alt_sum / dp_val
    # Skip trivial AFs
    if af == 0 or af == 1:
        return None

    try:
        pos = int(pos_str)
    except ValueError:
        return None

    span = max(1, len(ref))
    start0 = pos - 1
    end0 = start0 + span

    return BedGraphRecord(chrom=chrom, start=start0, end=end0, value=af)


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

        header_seen = False

        for line in fin:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                _parse_header_line(line)
                header_seen = True
                continue

            if not header_seen:
                raise click.ClickException("Malformed VCF: missing #CHROM header.")

            record = _extract_bedgraph_record(line.rstrip("\n").split("\t"))
            if record is None:
                continue

            fout.write(
                f"{record.chrom}\t{record.start}\t{record.end}\t{record.value:.6f}\n"
            )
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
