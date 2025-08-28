#!/usr/bin/env python
from __future__ import annotations

import sys
from typing import Dict, Iterable, List, Optional, Set, Tuple

import click

# --- Constants ---
INFO_HEADERS = {
    "AllowlistStatus": "Variant allowlisted based on CLNSIG/ONC or manually curated clinical list",
    "AllowlistedFilters": "Original filters for allowlisted variants that were overridden",
}
VCF_FIELDS = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"]  # (format unused here)

AllowKey = Tuple[str, int, str, str]  # (CHROM, POS, REF, ALT)


def parse_info(info_str: str) -> Dict[str, Optional[str]]:
    """Parse a VCF INFO string into a dict. Flags get value None. '.' -> {}."""
    if not info_str or info_str == ".":
        return {}
    out: Dict[str, Optional[str]] = {}
    for entry in info_str.split(";"):
        if not entry:
            continue
        if "=" in entry:
            k, v = entry.split("=", 1)
            out[k] = v
        else:
            out[entry] = None
    return out


def format_info(info: Dict[str, Optional[str]]) -> str:
    """Format INFO dict back into a VCF INFO string (sorted by key for stability)."""
    if not info:
        return "."
    parts: List[str] = []
    for k in sorted(info.keys()):
        v = info[k]
        parts.append(k if v is None else f"{k}={v}")
    return ";".join(parts)


def load_allowlist(allowlist_path: Optional[str]) -> Set[AllowKey]:
    """Read an allowlist VCF and return set of (chrom, pos, ref, alt) tuples."""
    allow: Set[AllowKey] = set()
    if not allowlist_path:
        return allow
    with open(allowlist_path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            chrom, pos_str, _id, ref, alt_field = cols[:5]
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            # ALT can have multiple alleles comma-separated
            for alt in alt_field.split(","):
                allow.add((chrom, pos, ref, alt))
    return allow


def needs_headers_injection(existing_header_lines: Iterable[str]) -> Dict[str, bool]:
    existing = "\n".join(existing_header_lines)
    return {k: (f"##INFO=<ID={k}," not in existing) for k in INFO_HEADERS.keys()}


def clinical_whitelist(info: Dict[str, Optional[str]]) -> Tuple[bool, List[str]]:
    reasons: List[str] = []
    onc_val = info.get("ONC")
    cln_val = info.get("CLNSIG")
    if onc_val and "oncogenic" in onc_val.lower():
        reasons.append("ONC=oncogenic")
    if cln_val and "pathogenic" in cln_val.lower():
        reasons.append("CLNSIG=pathogenic")
    return (len(reasons) > 0, reasons)


@click.command()
@click.argument("vcf_path", type=click.Path(exists=True))
@click.option(
    "--allowlist-path",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to allowlist VCF.",
)
@click.argument("output_file", type=click.Path())
def main(vcf_path: str, allowlist_path: Optional[str], output_file: str) -> None:
    """VCF processor with allowlist annotation support."""
    allow = load_allowlist(allowlist_path)

    header_lines: List[str] = []
    chrom_header: Optional[str] = None

    with open(vcf_path, "r", encoding="utf-8") as inp, open(output_file, "w", encoding="utf-8") as out:
        # First, buffer headers to decide what to inject
        data_lines: List[str] = []
        for raw in inp:
            if raw.startswith("##"):
                header_lines.append(raw.rstrip("\n"))
            elif raw.startswith("#CHROM"):
                chrom_header = raw.rstrip("\n")
                break
            else:
                # Malformed VCF with no #CHROM yet, treat as data
                data_lines.append(raw.rstrip("\n"))
                break

        if chrom_header is None:
            # Continue reading to find #CHROM if not yet found
            for raw in inp:
                if raw.startswith("##"):
                    header_lines.append(raw.rstrip("\n"))
                elif raw.startswith("#CHROM"):
                    chrom_header = raw.rstrip("\n")
                    break
                else:
                    data_lines.append(raw.rstrip("\n"))
                    break

        if chrom_header is None:
            click.echo("ERROR: No #CHROM header line found in VCF.", err=True)
            sys.exit(1)

        # Inject INFO headers if missing
        inject_flags = needs_headers_injection(header_lines)
        for h in header_lines:
            out.write(h + "\n")
        for key, missing in inject_flags.items():
            if missing:
                desc = INFO_HEADERS[key].replace('"', '\\"')
                out.write(f'##INFO=<ID={key},Number=1,Type=String,Description="{desc}">\n')

        # Write the #CHROM header
        out.write(chrom_header + "\n")

        # If we buffered any early data lines (edge-case), process them plus the rest
        def yield_data_lines() -> Iterable[str]:
            if data_lines:
                for dl in data_lines:
                    yield dl + "\n"
            for rest in inp:
                yield rest

        # Process variants
        for raw in yield_data_lines():
            if not raw or raw.startswith("#"):
                out.write(raw)
                continue

            line = raw.rstrip("\n")
            cols = line.split("\t")
            # Ensure we have at least 8 columns for VCF
            if len(cols) < 8:
                out.write(line + "\n")
                continue

            chrom, pos_str, vid, ref, alt_field, qual, filt, info_str = cols[:8]

            try:
                pos = int(pos_str)
            except ValueError:
                out.write(line + "\n")
                continue

            alts = alt_field.split(",")
            # Allowlist match?
            is_allow = any((chrom, pos, ref, alt) in allow for alt in alts)

            info = parse_info(info_str)
            has_clinical, clinical_reasons = clinical_whitelist(info)

            whitelisted = is_allow or has_clinical

            if whitelisted:
                reasons = []
                if is_allow:
                    reasons.append("allowlist")
                if has_clinical:
                    reasons.append("clinical(" + ",".join(clinical_reasons) + ")")
                # Move overridden filters (if any) to INFO
                if filt and filt != "PASS" and filt != ".":
                    # If AllowlistedFilters already exists, append unique values
                    prev = info.get("AllowlistedFilters")
                    merged = filt if not prev else ",".join(sorted(set((prev + "," + filt).split(","))))
                    info["AllowlistedFilters"] = merged
                # Add/overwrite status
                info["AllowlistStatus"] = "|".join(reasons) if reasons else "true"
                # Force PASS
                filt = "PASS"

            # Rebuild INFO and write line
            new_info = format_info(info)
            new_cols = cols[:]  # copy
            new_cols[6] = filt
            new_cols[7] = new_info
            out.write("\t".join(new_cols) + "\n")


if __name__ == "__main__":
    main()