#!/usr/bin/env python3

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import OrderedDict

import click
import vcfpy

INFO_ALLOWLISTED_FILTERS_ID = "AllowlistedFilters"
INFO_ALLOWLIST_STATUS_ID = "AllowlistStatus"

INFO_ALLOWLISTED_FILTERS_DESC = (
    "Original FILTER value moved here because this variant was allow-listed"
)
INFO_ALLOWLIST_STATUS_DESC = (
    "Reason(s) for allow-listing; pipe-separated (e.g., ClinvarOnc|ClinvarPathogenic)"
)

CLNSIG_PATHOGENIC = "Pathogenic"
CLNSIG_LIKELY_PATHOGENIC = "Likely_pathogenic"
CLINVAR_REASON_ONC = "ClinvarOnc"
CLINVAR_REASON_PATH = "ClinvarPathogenic"
CLINVAR_REASON_LIKELY_PATH = "ClinvarLikelyPathogenic"
MANUAL_REASON = "ManuallyCuratedList"


def _alt_to_str(alt) -> str:
    """
    Convert a vcfpy ALT object to its string representation.
    Handles Substitution, Breakend, Symbolic, and SV.
    """
    # Substitution: e.g., 'T'
    if hasattr(alt, "value"):
        return alt.value
    # Symbolic/SV/Breakend: e.g., '<DEL>' or 'N[chr1:123[' etc.
    if hasattr(alt, "to_string"):
        return alt.to_string()
    return str(alt)


def _normalize_info_value_to_str(value) -> str:
    """
    vcfpy will parse INFO values as scalars or lists depending on Number.
    This normalizes to a comma-joined string for searching substrings.
    """
    if value is None:
        return ""
    if isinstance(value, list):
        return ",".join("" if v is None else str(v) for v in value)
    return str(value)


def build_allowlist_keyset(allow_vcf: Path) -> Set[Tuple[str, int, str, str]]:
    """
    Read a VCF (gz ok) via vcfpy and return a set of (CHROM, POS, REF, ALT).
    One entry per ALT allele in multi-allelic records.
    """
    keys: Set[Tuple[str, int, str, str]] = set()
    with vcfpy.Reader.from_path(str(allow_vcf)) as rdr:
        for rec in rdr:
            chrom = rec.CHROM
            pos = rec.POS
            ref = rec.REF
            for alt in rec.ALT:
                keys.add((chrom, pos, ref, _alt_to_str(alt)))
    return keys


def ensure_info_headers(header: vcfpy.Header) -> None:
    """
    Ensure the custom INFO headers exist; add if missing.
    """
    info_ids = {line.id for line in header.get_lines("INFO")}
    if INFO_ALLOWLISTED_FILTERS_ID not in info_ids:
        header.add_info_line(
            OrderedDict(
                [
                    ("ID", INFO_ALLOWLISTED_FILTERS_ID),
                    ("Number", "1"),
                    ("Type", "String"),
                    ("Description", INFO_ALLOWLISTED_FILTERS_DESC),
                ]
            )
        )
    if INFO_ALLOWLIST_STATUS_ID not in info_ids:
        header.add_info_line(
            OrderedDict(
                [
                    ("ID", INFO_ALLOWLIST_STATUS_ID),
                    ("Number", "1"),
                    ("Type", "String"),
                    ("Description", INFO_ALLOWLIST_STATUS_DESC),
                ]
            )
        )


def determine_clinvar_reasons(info: Dict[str, object]) -> List[str]:
    """
    Determine allow-list reasons based on ClinVar annotations.
    - If INFO has ONC flag/present → ClinvarOnc.
    - If CLNSIG contains 'Pathogenic' → ClinvarPathogenic.
    - If CLNSIG contains 'Likely_pathogenic' → ClinvarLikelyPathogenic.
    """
    reasons: List[str] = []

    # ONC can be a flag (present with True/None) or a value; presence is enough
    if "ONC" in info:
        val = info.get("ONC")
        if val is True or val is None or val == "" or val:  # any presence
            reasons.append(CLINVAR_REASON_ONC)

    clnsig = _normalize_info_value_to_str(info.get("CLNSIG"))
    if CLNSIG_PATHOGENIC in clnsig:
        reasons.append(CLINVAR_REASON_PATH)
    if CLNSIG_LIKELY_PATHOGENIC in clnsig:
        reasons.append(CLINVAR_REASON_LIKELY_PATH)

    return reasons


def any_alt_in_allowlist(
    allow_keys: Set[Tuple[str, int, str, str]],
    rec: vcfpy.Record,
) -> bool:
    for alt in rec.ALT:
        if (rec.CHROM, rec.POS, rec.REF, _alt_to_str(alt)) in allow_keys:
            return True
    return False


def record_filter_string(rec: vcfpy.Record) -> str:
    """
    Turn the record's FILTER field into a VCF string representation.
    - PASS → "PASS"
    - no filters ('.') → "."
    - otherwise semicolon-joined filters.
    """
    filt = rec.FILTER or []
    if not filt:
        # vcfpy uses [] to mean PASS in writing, but for our logic we
        # need to distinguish '.' (no filters applied) from 'PASS'.
        # We check if original raw was '.' is not directly available here.
        # Practical approach: treat empty list as PASS for writing,
        # and rely on rec.INFO.get('AllowlistedFilters') existence when we set it.
        return "PASS"
    if filt == ["PASS"]:
        return "PASS"
    return ";".join(filt)


def set_record_pass(rec: vcfpy.Record) -> None:
    """
    Force FILTER to PASS in a way vcfpy will serialize as 'PASS'.
    """
    rec.FILTER = ["PASS"]


def process_vcf(
    allow_keys: Optional[Set[Tuple[str, int, str, str]]],
    in_vcf: Path,
    out_path: Optional[Path],
) -> None:
    """
    Read, annotate, and write VCF using vcfpy.
    - Adds headers if missing.
    - Adds INFO/AllowlistStatus when reasons exist.
    - If original FILTER is not PASS/'.', moves it to INFO/AllowlistedFilters and sets FILTER to PASS.
    """
    with vcfpy.Reader.from_path(str(in_vcf)) as reader:
        header = reader.header.copy()
        ensure_info_headers(header)

        # Writer: stdout if '-' else to path (gz supported by vcfpy by extension)
        if out_path is None or str(out_path) == "-":
            writer = vcfpy.Writer.from_stream(sys.stdout, header)
        else:
            writer = vcfpy.Writer.from_path(str(out_path), header)

        with writer:
            for rec in reader:
                reasons: List[str] = []

                # Manual allow-list
                if allow_keys and any_alt_in_allowlist(allow_keys, rec):
                    reasons.append(MANUAL_REASON)

                # ClinVar-based
                reasons.extend(determine_clinvar_reasons(rec.INFO))

                if reasons:
                    # Move FILTER -> INFO/AllowlistedFilters if filtered
                    filt_str = record_filter_string(rec)
                    if filt_str not in ("PASS", ".", ""):
                        rec.INFO[INFO_ALLOWLISTED_FILTERS_ID] = filt_str
                        set_record_pass(rec)

                    # Set AllowlistStatus
                    rec.INFO[INFO_ALLOWLIST_STATUS_ID] = "|".join(reasons)

                writer.write_record(rec)



@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--allow-list",
    "allow_list",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    required=False,
    help="Optional: VCF (.vcf or .vcf.gz) of allow-listed variants (matched on CHROM, POS, REF, ALT).",
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
def cli(allow_list: Path | None, vcf_path: Path, out_path: Path | None) -> None:
    """Annotate variants with INFO/AllowlistStatus and optionally INFO/AllowlistedFilters (vcfpy-based)."""
    try:
        allow_keys: Optional[Set[Tuple[str, int, str, str]]] = None
        if allow_list:
            allow_keys = build_allowlist_keyset(allow_list)
        process_vcf(allow_keys, vcf_path, out_path)
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