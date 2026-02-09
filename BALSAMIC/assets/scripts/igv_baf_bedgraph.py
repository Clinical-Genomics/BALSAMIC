#!/usr/bin/env python3
from __future__ import annotations

import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, TextIO

import click
import numpy as np
from cyvcf2 import VCF


def open_output(path: str | Path) -> TextIO:
    """Open output file; '-' means stdout."""
    p = str(path)
    if p == "-":
        return sys.stdout
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
    # symbolic alleles like <DEL>, <NON_REF>, etc.
    if alt.startswith("<") and alt.endswith(">"):
        return True
    return False


def _sample_ad_dp(variant) -> Optional[tuple[int, int]]:
    """
    Return (ref_ad, alt_sum) for the first (and only) sample using FORMAT/AD.

    Returns None if AD is missing/malformed.
    """
    try:
        ad = variant.format("AD")
    except Exception:
        return None

    if ad is None:
        return None

    # Expected shape: (n_samples, n_alleles)
    if not isinstance(ad, np.ndarray) or ad.ndim != 2 or ad.shape[0] < 1:
        return None

    counts = ad[0]

    # cyvcf2/htslib uses sentinel for missing; treat any negative as missing
    # (most commonly -2147483648)
    if counts.size < 2 or np.any(counts < 0):
        return None

    ref_ad = int(counts[0])
    alt_sum = int(np.sum(counts[1:]))
    return ref_ad, alt_sum


def _sample_dp(variant) -> Optional[int]:
    """
    Return DP for the first sample if present and sane, else None.
    """
    # Prefer FORMAT/DP if present
    try:
        dp = variant.format("DP")
    except Exception:
        dp = None

    if isinstance(dp, np.ndarray) and dp.size >= 1:
        v = int(dp.reshape(-1)[0])
        if v > 0:
            return v

    # Fallback: cyvcf2 convenience (may be -1 if missing)
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
    """
    Convert a cyvcf2 Variant into a BedGraphRecord, or None to skip.
    """
    # Skip any record with symbolic/missing ALT(s)
    # Keep behavior similar to original script: skip if ALT is "." or "*" (and also <...>)
    alts = variant.ALT or []
    if not alts:
        return None
    if any(_is_symbolic_or_missing_alt(a) for a in alts):
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


def _open_vcf_for_streaming(
    vcf_path: str | Path,
) -> tuple[VCF, Optional[tempfile.NamedTemporaryFile]]:
    """
    Open a VCF via cyvcf2. If vcf_path == '-', buffer stdin to a temp file.
    Returns (vcf_reader, temp_file_or_None).
    """
    p = str(vcf_path)
    tmp = None

    if p == "-":
        # cyvcf2 expects a filename; buffer stdin to a temporary VCF file
        tmp = tempfile.NamedTemporaryFile(mode="wb", suffix=".vcf", delete=True)
        data = sys.stdin.buffer.read()
        tmp.write(data)
        tmp.flush()
        vcf = VCF(tmp.name)
        return vcf, tmp

    vcf = VCF(p)
    return vcf, None


def convert_vcf_to_bedgraph(
    vcf_path: str | Path,
    bedgraph_path: str | Path,
    track_name: Optional[str] = None,
) -> int:
    """
    Convert a VCF into a bedGraph of AF from AD/DP.
    Assumes exactly one sample is present.
    """
    vcf, tmp = _open_vcf_for_streaming(vcf_path)
    try:
        n_written = 0
        with open_output(bedgraph_path) as fout:
            if track_name:
                fout.write(f'track type=bedGraph name="{track_name}"\n')

            for variant in vcf:
                rec = _variant_to_bedgraph_record(variant)
                if rec is None:
                    continue
                fout.write(f"{rec.chrom}\t{rec.start}\t{rec.end}\t{rec.value:.6f}\n")
                n_written += 1

        return n_written
    finally:
        try:
            vcf.close()
        except Exception:
            pass
        if tmp is not None:
            try:
                tmp.close()
            except Exception:
                pass


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
        vcf_path=vcf, bedgraph_path=bedgraph, track_name=track_name
    )
    click.echo(f"Done. Wrote {n} rows → {bedgraph}", err=True)


if __name__ == "__main__":
    cli()
