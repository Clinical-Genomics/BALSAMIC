#!/usr/bin/env python3
from __future__ import annotations

import gzip
import io
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import click

# --- Constants & headers ------------------------------------------------------

INFO_ALLOWLISTED_FILTERS_HDR = '##INFO=<ID=AllowlistedFilters,Number=1,Type=String,Description="Original FILTER value moved here because this variant was allow-listed">'
INFO_ALLOWLIST_STATUS_HDR = '##INFO=<ID=AllowlistStatus,Number=1,Type=String,Description="Reason(s) for allow-listing; pipe-separated (e.g., ManuallyCuratedList|ClinvarOnc|ClinvarPathogenic)">'

CLNSIG_PATHOGENIC = "Pathogenic"
CLNSIG_LIKELY_PATHOGENIC = "Likely_pathogenic"
CLINVAR_REASON_ONC = "ClinvarOnc"
CLINVAR_REASON_PATH = "ClinvarPathogenic"
CLINVAR_REASON_LIKELY_PATH = "ClinvarLikelyPathogenic"
MANUAL_REASON = "ManuallyCuratedList"


# --- IO helpers ---------------------------------------------------------------


def open_maybe_gzip(path: str | Path) -> io.TextIOBase:
    p = str(path)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", newline="")
    return open(p, "r", encoding="utf-8", newline="")


def open_out_text(path: str | Path | None) -> io.TextIOBase:
    if path is None or str(path) == "-":
        return sys.stdout
    return open(path, "w", encoding="utf-8", newline="")


# --- VCF parsing utilities ----------------------------------------------------


def parse_info(info_field: str) -> Dict[str, Optional[str]]:
    if info_field == "." or not info_field:
        return {}
    out: Dict[str, Optional[str]] = {}
    for entry in info_field.split(";"):
        if not entry:
            continue
        if "=" in entry:
            k, v = entry.split("=", 1)
            out[k] = v
        else:
            out[entry] = None
    return out


def format_info(info: Dict[str, Optional[str]]) -> str:
    if not info:
        return "."
    parts: List[str] = []
    for k in sorted(info.keys()):
        v = info[k]
        parts.append(f"{k}={v}" if v is not None else k)
    return ";".join(parts)


# --- Allow-list building & reasons -------------------------------------------


def build_allowlist_keyset(
    allow_vcf_path: str | Path,
) -> Set[Tuple[str, int, str, str]]:
    keyset: Set[Tuple[str, int, str, str]] = set()
    with open_maybe_gzip(allow_vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, pos_str, _id, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            for a in alt.split(","):
                keyset.add((chrom, pos, ref, a))
    return keyset


def determine_clinvar_reasons(info: Dict[str, Optional[str]]) -> List[str]:
    """
    Return ClinVar allow-list reasons as separate labels (no merging):
      - ONC present -> ClinvarOnc
      - CLNSIG contains 'Pathogenic' (case-sensitive) -> ClinvarPathogenic
      - CLNSIG contains 'Likely_pathogenic' (case-sensitive) -> ClinvarLikelyPathogenic
    """
    reasons: List[str] = []
    if "ONC" in info:
        reasons.append(CLINVAR_REASON_ONC)
    clnsig = info.get("CLNSIG") or ""
    if CLNSIG_PATHOGENIC in clnsig:
        reasons.append(CLINVAR_REASON_PATH)
    if CLNSIG_LIKELY_PATHOGENIC in clnsig:
        reasons.append(CLINVAR_REASON_LIKELY_PATH)
    return reasons


def ensure_info_headers(headers: List[str]) -> List[str]:
    has_filters = any(h.startswith("##INFO=<ID=AllowlistedFilters,") for h in headers)
    has_status = any(h.startswith("##INFO=<ID=AllowlistStatus,") for h in headers)
    if has_filters and has_status:
        return headers
    try:
        chrom_idx = next(i for i, h in enumerate(headers) if h.startswith("#CHROM"))
    except StopIteration:
        chrom_idx = len(headers)
    insertion: List[str] = []
    if not has_filters:
        insertion.append(INFO_ALLOWLISTED_FILTERS_HDR)
    if not has_status:
        insertion.append(INFO_ALLOWLIST_STATUS_HDR)
    return headers[:chrom_idx] + insertion + headers[chrom_idx:]


# --- Core processing ----------------------------------------------------------


def process_vcf(
    allow_keys: Set[Tuple[str, int, str, str]],
    in_vcf_path: str | Path,
    out_fh: io.TextIOBase,
) -> None:
    headers: List[str] = []
    header_done = False

    with open_maybe_gzip(in_vcf_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if not header_done:
                if line.startswith("#"):
                    headers.append(line)
                    continue
                headers = ensure_info_headers(headers)
                for h in headers:
                    out_fh.write(h + "\n")
                header_done = True  # fall through

            cols = line.split("\t")
            if len(cols) < 8:
                out_fh.write(line + "\n")
                continue

            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError:
                out_fh.write(line + "\n")
                continue
            ref = cols[3]
            alts = cols[4].split(",") if cols[4] else ["."]
            filt = cols[6]
            info = parse_info(cols[7])

            # Reasons
            reasons: List[str] = []
            if any((chrom, pos, ref, a) in allow_keys for a in alts):
                reasons.append(MANUAL_REASON)
            reasons.extend(determine_clinvar_reasons(info))

            if reasons:
                if filt not in ("PASS", ".", ""):
                    info["AllowlistedFilters"] = filt
                    cols[6] = "PASS"
                info["AllowlistStatus"] = "|".join(reasons)
                cols[7] = format_info(info)

            out_fh.write("\t".join(cols) + "\n")

    if not header_done and headers:
        headers = ensure_info_headers(headers)
        for h in headers:
            out_fh.write(h + "\n")


# --- CLI ----------------------------------------------------------------------


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "--allow-list",
    "allow_list",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    required=True,
    help="VCF (.vcf or .vcf.gz) of allow-listed variants (matched on CHROM, POS, REF, ALT).",
)
@click.option(
    "--vcf",
    "vcf_path",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    required=True,
    help="Input VCF (.vcf or .vcf.gz) to process.",
)
@click.option(
    "-o",
    "--out",
    "out_path",
    type=click.Path(dir_okay=False, writable=True, allow_dash=True, path_type=Path),
    default="-",
    show_default=True,
    help="Output VCF path (use '-' for stdout).",
)
def cli(allow_list: Path, vcf_path: Path, out_path: Path | None) -> None:
    """
    Keep ONC as its own allow-list reason (ClinvarOnc) alongside Pathogenic/Likely_pathogenic:
      AllowlistStatus example: ManuallyCuratedList|ClinvarOnc|ClinvarPathogenic
    """
    try:
        allow_keys = build_allowlist_keyset(allow_list)
        with open_out_text(out_path) as out_fh:
            process_vcf(allow_keys, vcf_path, out_fh)
    except BrokenPipeError:
        try:
            sys.stdout.close()
        finally:
            sys.exit(0)
    except Exception as e:
        click.echo(f"[error] {e}", err=True)
        sys.exit(1)


if __name__ == "__main__":
    cli()
