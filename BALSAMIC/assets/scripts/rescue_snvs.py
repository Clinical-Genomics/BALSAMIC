#!/usr/bin/env python3
from __future__ import annotations

import gzip
import io
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import click

# --- Constants & headers ------------------------------------------------------

INFO_RESCUE_FILTERS_HDR = '##INFO=<ID=RescueFilters,Number=1,Type=String,Description="Original FILTER value moved here because this variant was rescued">'
INFO_RESCUE_STATUS_HDR = '##INFO=<ID=RescueStatus,Number=1,Type=String,Description="Reason(s) for rescuing variant; pipe-separated (e.g., ClinvarOnc|ClinvarPathogenic)">'

CLNSIG_PATHOGENIC = "Pathogenic"
CLNSIG_LIKELY_PATHOGENIC = "Likely_pathogenic"
CLINVAR_REASON_ONC = "ClinvarOnc"
CLINVAR_REASON_PATH = "ClinvarPathogenic"
CLINVAR_REASON_LIKELY_PATH = "ClinvarLikelyPathogenic"
MANUAL_REASON = "ManuallyCuratedList"


def open_maybe_gzip(path: str | Path) -> io.TextIOBase:
    """
    Open a text file that may optionally be gzip-compressed.

    - If the file ends with ".gz", it is opened with gzip in binary mode
      and wrapped with a TextIOWrapper for UTF-8 text decoding.
    - Otherwise, the file is opened normally in text mode.
    """
    p = str(path)
    if p.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", newline="")
    return open(p, "r", encoding="utf-8", newline="")


def open_out_text(path: str | Path | None) -> io.TextIOBase:
    """
    Open a writable text stream for output.

    - If path is None or "-", return sys.stdout (for writing to console).
    - Otherwise, open the file at the given path in write mode ("w"),
      with UTF-8 encoding.
    """
    if path is None or str(path) == "-":
        return sys.stdout
    return open(path, "w", encoding="utf-8", newline="")


def parse_info(info_field: str) -> Dict[str, Optional[str]]:
    """
    Parse a VCF INFO field into a dictionary.

    - Splits on semicolons (";") to get entries.
    - Each entry is either:
      * "KEY=VALUE" → stored as {KEY: VALUE}
      * "FLAG" (no "=") → stored as {FLAG: None}
    - Special case: if the field is "." or empty, returns an empty dict.
    """
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
    """
    Format a dictionary back into a VCF INFO field string.

    - Keys are sorted alphabetically.
    - Each item is rendered as "KEY=VALUE" if a value is present,
      or just "KEY" if the value is None.
    - If the dictionary is empty, returns ".".
    """
    if not info:
        return "."
    parts: List[str] = []
    for k in sorted(info.keys()):
        v = info[k]
        parts.append(f"{k}={v}" if v is not None else k)
    return ";".join(parts)


def build_rescue_keyset(
    rescue_vcf_path: str | Path,
) -> Set[Tuple[str, int, str, str]]:
    """
    Build a set of (CHROM, POS, REF, ALT) keys from a VCF file.

    - Reads through the VCF, skipping headers and malformed lines.
    - For each valid record, extracts chrom, pos, ref, and each alt allele.
    - Adds one tuple per allele to the returned set.
    """
    keyset: Set[Tuple[str, int, str, str]] = set()
    with open_maybe_gzip(rescue_vcf_path) as fh:
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
    Determine rescue reasons based on ClinVar annotations.

    - If INFO contains "ONC" (oncogenic flag), adds `CLINVAR_REASON_ONC`.
    - If INFO["CLNSIG"] contains "pathogenic", adds `CLINVAR_REASON_PATH`.
    - If INFO["CLNSIG"] contains "likely_pathogenic", adds
      `CLINVAR_REASON_LIKELY_PATH`.
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
    """
    Ensure VCF headers include Rescue INFO field definitions.

    - Looks for existing INFO headers for `RescueFilters`
      and `RescueStatus`.
    - If missing, inserts them before the `#CHROM` line (or at the end if
      no `#CHROM` is found).
    - Returns a new list of headers (does not mutate in place).

    Args:
        headers: List of VCF header lines (including "##" and "#CHROM").

    Returns:
        Updated header list with required INFO definitions present.
    """
    has_filters = any(h.startswith("##INFO=<ID=RescueFilters,") for h in headers)
    has_status = any(h.startswith("##INFO=<ID=RescueStatus,") for h in headers)
    if has_filters and has_status:
        return headers
    try:
        chrom_idx = next(i for i, h in enumerate(headers) if h.startswith("#CHROM"))
    except StopIteration:
        chrom_idx = len(headers)
    insertion: List[str] = []
    if not has_filters:
        insertion.append(INFO_RESCUE_FILTERS_HDR)
    if not has_status:
        insertion.append(INFO_RESCUE_STATUS_HDR)
    return headers[:chrom_idx] + insertion + headers[chrom_idx:]


def process_vcf(
    rescue_keys: Optional[Set[Tuple[str, int, str, str]]],
    in_vcf_path: str | Path,
    out_fh: io.TextIOBase,
) -> None:
    """
    Stream a VCF, rescue variants and normalizing FILTER/INFO.

    Writes headers once (after reading them) and then writes each processed record.
    If the file contains only headers, they’re written at the end (same as before).
    """
    headers: Optional[List[str]] = []

    with open_maybe_gzip(in_vcf_path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if line.startswith("#"):
                # Still in the header section.
                assert headers is not None
                headers.append(line)
                continue

            # First non-header line: flush headers once.
            if headers is not None:
                for h in ensure_info_headers(headers):
                    out_fh.write(h + "\n")
                headers = None  # mark as flushed

            out_fh.write(_process_record_line(line, rescue_keys))

    # File had only headers (no records): write them now.
    if headers is not None:
        for h in ensure_info_headers(headers):
            out_fh.write(h + "\n")


def _process_record_line(
    line: str,
    rescue_keys: Optional[Set[Tuple[str, int, str, str]]],
) -> str:
    """Return the (possibly modified) VCF record line with trailing newline."""
    cols = line.split("\t")
    if len(cols) < 8:
        return line + "\n"

    chrom = cols[0]
    try:
        pos = int(cols[1])
    except ValueError:
        return line + "\n"

    ref = cols[3]
    alts = cols[4].split(",") if cols[4] else ["."]
    filt = cols[6]
    info = parse_info(cols[7])

    reasons: List[str] = []

    # Manual rescue (only if provided).
    if rescue_keys and _any_alt_in_rescue(rescue_keys, chrom, pos, ref, alts):
        reasons.append(MANUAL_REASON)

    # ClinVar rescue.
    reasons.extend(determine_clinvar_reasons(info))

    if not reasons:
        return line + "\n"

    # Update FILTER/INFO when rescued.
    if filt not in ("PASS", ".", ""):
        info["RescueFilters"] = filt
        cols[6] = "PASS"
    info["RescueStatus"] = "|".join(reasons)
    cols[7] = format_info(info)

    return "\t".join(cols) + "\n"


def _any_alt_in_rescuelist(
    rescue_keys: Set[Tuple[str, int, str, str]],
    chrom: str,
    pos: int,
    ref: str,
    alts: List[str],
) -> bool:
    return any((chrom, pos, ref, alt) in rescue_keys for alt in alts)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--rescue-list",
    "rescue_list",
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=Path),
    required=False,
    help="Optional: VCF (.vcf or .vcf.gz) of variants to rescue (matched on CHROM, POS, REF, ALT).",
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
    """
    Annotate variants with INFO/RescueStatus and optionally INFO/RescueFilters.

    Rescue criteria:
      - If --rescue-list is provided: manual rescue variant match (CHROM, POS, REF, ALT).
      - ClinVar: ONC present -> ClinvarOnc; CLNSIG contains 'Pathogenic' or 'Likely_pathogenic'.
    """
    try:
        rescue_keys: Set[Tuple[str, int, str, str]] | None = None
        if rescue_list:
            rescue_keys = build_rescuelist_keyset(rescue_list)
        with open_out_text(out_path) as out_fh:
            process_vcf(rescue_keys, vcf_path, out_fh)
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
