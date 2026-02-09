#!/usr/bin/env python3
from __future__ import annotations

import sys
import tempfile
from contextlib import closing, nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, TextIO, ContextManager

import click
import numpy as np
from cyvcf2 import VCF


def open_output(path: str | Path) -> ContextManager[TextIO]:
    """Return a context manager for output; '-' means stdout (not closed)."""
    p = str(path)
    if p == "-":
        return nullcontext(sys.stdout)
    return open(p, "w", encoding="utf-8")


@dataclass(frozen=True)
class BedGraphRecord:
    chrom: str
    start: int
    end: int
    value: float


def _is_symbolic_or_missing_alt(alt: Optional[str]) -> bool:
    """Return True for missing ALT, '*' no-call, or symbolic alleles like <DEL>."""
    if alt is None:
        return True
    if alt in (".", "*"):
        return True
    return alt.startswith("<") and alt.endswith(">")


def _sample_ad_dp(variant) -> Optional[tuple[int, int]]:
    """Return (ref_ad, alt_sum) from FORMAT/AD for the first sample, else None."""
    try:
        ad = variant.format("AD")
    except Exception:
        return None

    if ad is None:
        return None

    if not isinstance(ad, np.ndarray) or ad.ndim != 2 or ad.shape[0] < 1:
        return None

    counts = ad[0]
    if counts.size < 2 or np.any(counts < 0):
        return None

    ref_ad = int(counts[0])
    alt_sum = int(np.sum(counts[1:]))
    return ref_ad, alt_sum


def _sample_dp(variant) -> Optional[int]:
    """Return DP for the first sample if present and sane, else None."""
    try:
        dp = variant.format("DP")
    except Exception:
        dp = None

    if isinstance(dp, np.ndarray) and dp.size >= 1:
        v = int(dp.reshape(-1)[0])
        if v > 0:
            return v

    try:
        depths = variant.gt_depths
    except Exception:
        depths = None

    if depths is not None and len(depths) >= 1:
        v = int(depths[0])
        if v > 0:
            return v

    return None


def _variant_to_bedgraph_record(variant) -> Optional[BedGraphRecord]:
    """Convert a cyvcf2 Variant into a BedGraphRecord, or None to skip."""
    alts = variant.ALT or []
    if not alts or any(_is_symbolic_or_missing_alt(a) for a in alts):
        return None

    ad_pair = _sample_ad_dp(variant)
    if ad_pair is None:
        return None
    ref_ad, alt_sum = ad_pair

    dp_val = _sample_dp(variant)
    if dp_val is None or dp_val <= 0:
        dp_val = ref_ad + alt_sum
    if dp_val <= 0:
        return None

    af = alt_sum / dp_val
    if af == 0.0 or af == 1.0:
        return None

    start0 = int(variant.POS) - 1
    span = max(1, len(variant.REF or "N"))
    end0 = start0 + span

    return BedGraphRecord(
        chrom=str(variant.CHROM),
        start=start0,
        end=end0,
        value=float(af),
    )


def convert_vcf_to_bedgraph(
    vcf_path: str | Path,
    bedgraph_path: str | Path,
    track_name: Optional[str] = None,
) -> int:
    """
    Convert a VCF into a bedGraph of AF from AD/DP.
    Assumes exactly one sample is present.
    """
    p = str(vcf_path)

    if p == "-":
        # cyvcf2 expects a filename; buffer stdin to a temporary VCF file
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".vcf") as tmp:
            tmp.write(sys.stdin.buffer.read())
            tmp.flush()

            with closing(VCF(tmp.name)) as vcf, open_output(bedgraph_path) as fout:
                return _stream_write_bedgraph(vcf, fout, track_name)

    with closing(VCF(p)) as vcf, open_output(bedgraph_path) as fout:
        return _stream_write_bedgraph(vcf, fout, track_name)


def _stream_write_bedgraph(vcf: VCF, fout: TextIO, track_name: Optional[str]) -> int:
    """Write bedGraph records from an open cyvcf2 VCF reader to fout."""
    n_written = 0
    if track_name:
        fout.write(f'track type=bedGraph name="{track_name}"\n')

    for variant in vcf:
        rec = _variant_to_bedgraph_record(variant)
        if rec is None:
            continue
        fout.write(f"{rec.chrom}\t{rec.start}\t{rec.end}\t{rec.value:.6f}\n")
        n_written += 1

    return n_written


@click.command(context_settings={"show_default": True})
@click.argument("vcf", type=click.Path(exists=True, dir_okay=False, allow_dash=True))
@click.argument("bedgraph", type=click.Path(dir_okay=False, allow_dash=True))
@click.option("--track-name", help='Optional bedGraph "track" name header.')
def cli(vcf: str, bedgraph: str, track_name: Optional[str]) -> None:
    """Convert VCF (with 1 sample) to bedGraph of AF from AD/DP."""
    n = convert_vcf_to_bedgraph(
        vcf_path=vcf, bedgraph_path=bedgraph, track_name=track_name
    )
    click.echo(f"Done. Wrote {n} rows → {bedgraph}", err=True)


if __name__ == "__main__":
    cli()
