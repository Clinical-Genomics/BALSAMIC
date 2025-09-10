#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import OrderedDict

import click
import vcfpy


INFO_RESCUE_FILTERS_ID = "RescueFilters"
INFO_RESCUE_STATUS_ID = "RescueStatus"

INFO_RESCUE_FILTERS_DESC = (
    "Original FILTER value moved here because this variant was rescue-listed"
)
INFO_RESCUE_STATUS_DESC = (
    "Reason(s) for rescue-listing; pipe-separated (e.g., ClinvarOnc|ClinvarPathogenic)"
)

CLNSIG_PATHOGENIC = "Pathogenic"
CLNSIG_LIKELY_PATHOGENIC = "Likely_pathogenic"
CLINVAR_REASON_ONC = "ClinvarOnc"
CLINVAR_REASON_PATH = "ClinvarPathogenic"
CLINVAR_REASON_LIKELY_PATH = "ClinvarLikelyPathogenic"
MANUAL_REASON = "ManuallyCuratedList"


def _alt_to_str(alt) -> str:
    """Convert a vcfpy ALT object to a string (handles Substitution, Symbolic/SV, Breakend)."""
    if hasattr(alt, "value"):  # Substitution
        return alt.value
    if hasattr(alt, "to_string"):  # Symbolic/SV/Breakend
        return alt.to_string()
    return str(alt)


def _info_value_as_str(value) -> str:
    """Normalize scalar/list INFO values to a comma-joined string for substring checks."""
    if value is None:
        return ""
    if isinstance(value, list):
        return ",".join("" if v is None else str(v) for v in value)
    return str(value)


def _filter_as_string(rec: vcfpy.Record) -> str:
    """
    Render FILTER as it *originally means* in VCF:
    - PASS           -> "PASS"      (vcfpy: usually ['PASS'])
    - no filters '.' -> "."         (vcfpy: usually [])
    - named filters  -> "f1;f2"     (vcfpy: ['f1','f2'])
    """
    flt = rec.FILTER or []
    if flt == ["PASS"]:
        return "PASS"
    if len(flt) == 0:
        return "."
    return ";".join(flt)


def _set_pass(rec: vcfpy.Record) -> None:
    """Force FILTER to PASS."""
    rec.FILTER = ["PASS"]


def build_rescue_keyset(rescue_vcf: Path) -> Set[Tuple[str, int, str, str]]:
    """Read a VCF and return a set of (CHROM, POS, REF, ALT) tuples (per ALT)."""
    keys: Set[Tuple[str, int, str, str]] = set()
    with vcfpy.Reader.from_path(str(rescue_vcf)) as rdr:
        for rec in rdr:
            for alt in rec.ALT:
                keys.add((rec.CHROM, rec.POS, rec.REF, _alt_to_str(alt)))
    return keys


def any_alt_in_rescue(
    rescue_keys: Set[Tuple[str, int, str, str]],
    rec: vcfpy.Record,
) -> bool:
    return any(
        (rec.CHROM, rec.POS, rec.REF, _alt_to_str(alt)) in rescue_keys for alt in rec.ALT
    )


def ensure_info_headers(header: vcfpy.Header) -> None:
    """Ensure our two INFO headers exist (add if missing)."""
    info_ids = {line.id for line in header.get_lines("INFO")}
    if INFO_RESCUE_FILTERS_ID not in info_ids:
        header.add_info_line(
            OrderedDict(
                ID=INFO_RESCUE_FILTERS_ID,
                Number="1",
                Type="String",
                Description=INFO_RESCUE_FILTERS_DESC,
            )
        )
    if INFO_RESCUE_STATUS_ID not in info_ids:
        header.add_info_line(
            OrderedDict(
                ID=INFO_RESCUE_STATUS_ID,
                Number="1",
                Type="String",
                Description=INFO_RESCUE_STATUS_DESC,
            )
        )


def determine_clinvar_reasons(info: Dict[str, object]) -> List[str]:
    """
    Return ClinVar-based reasons:
      - ONC present → ClinvarOnc
      - CLNSIG contains 'Pathogenic' → ClinvarPathogenic
      - CLNSIG contains 'Likely_pathogenic' → ClinvarLikelyPathogenic
    """
    reasons: List[str] = []

    if "ONC" in info:
        # treat presence of ONC (flag or any value) as true
        reasons.append(CLINVAR_REASON_ONC)

    clnsig = _info_value_as_str(info.get("CLNSIG"))
    if CLNSIG_PATHOGENIC in clnsig:
        reasons.append(CLINVAR_REASON_PATH)
    if CLNSIG_LIKELY_PATHOGENIC in clnsig:
        reasons.append(CLINVAR_REASON_LIKELY_PATH)

    return reasons


def process_vcf(
    rescue_keys: Optional[Set[Tuple[str, int, str, str]]],
    in_vcf: Path,
    out_path: Optional[Path],
) -> None:
    """
    Read, annotate, and write VCF using vcfpy.

    Policy change implemented here:
      - We consider ONLY literal 'PASS' as pass.
      - '.' (no filters applied) is treated as *not PASS*.
      - When a record is rescue-listed and its FILTER != 'PASS',
        we move the original FILTER string ('.' or 'f1;f2;...') to INFO/RescuedFilters
        and set FILTER to 'PASS'.
    """
    with vcfpy.Reader.from_path(str(in_vcf)) as reader:
        header = reader.header.copy()
        ensure_info_headers(header)

        writer = (
            vcfpy.Writer.from_stream(sys.stdout, header)
            if out_path is None or str(out_path) == "-"
            else vcfpy.Writer.from_path(str(out_path), header)
        )

        with writer:
            for rec in reader:
                reasons: List[str] = []

                if rescue_keys and any_alt_in_rescue(rescue_keys, rec):
                    reasons.append(MANUAL_REASON)

                reasons.extend(determine_clinvar_reasons(rec.INFO))

                if reasons:
                    # Only 'PASS' is pass; '.' and any named filters are moved.
                    original_filter = _filter_as_string(rec)
                    if original_filter != "PASS":
                        rec.INFO[INFO_RESCUE_FILTERS_ID] = original_filter
                        _set_pass(rec)

                    rec.INFO[INFO_RESCUE_STATUS_ID] = "|".join(reasons)

                writer.write_record(rec)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--rescue-list",
    "rescue_list",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    required=False,
    help="Optional: VCF (.vcf or .vcf.gz) of rescue-listed variants (match on CHROM, POS, REF, ALT).",
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
def cli(rescue_list: Path | None, vcf_path: Path, out_path: Path | None) -> None:
    """Annotate variants with INFO/RescueStatus and move any non-PASS FILTER into INFO when rescue-listed."""
    try:
        rescue_keys: Optional[Set[Tuple[str, int, str, str]]] = None
        if rescue_list:
            rescue_keys = build_rescue_keyset(rescue_list)
        process_vcf(rescue_keys, vcf_path, out_path)
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
