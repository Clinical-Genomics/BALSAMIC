#!/usr/bin/env python
from typing import Dict, List, Tuple
import click
import gzip
import csv

# --- Constants ---
INFO_HEADERS = {
    "AllowlistStatus": "Variant allowlisted based on CLNSIG/ONC or manually curated clinical list",
    "AllowlistedFilters": "Original filters for allowlisted variants that were overridden",
}
VCF_FIELDS = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"]


# --- Utility functions ---
def parse_info(info_str: str) -> Dict[str, str]:
    """Parse a semicolon-separated VCF INFO string into a dictionary."""
    info_dict = {}
    for field in info_str.split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            info_dict[key] = value
    return info_dict


def partial_match(query: str, target: str) -> bool:
    """Check for exact or partial match, supporting '...' wildcard."""
    return query.strip("...") in target if "..." in query else query == target


# --- Core functionality ---
def onc_clnsig(info: str) -> bool:
    """Check if INFO string contains CLNSIG or ONC evidence."""
    info_fields = parse_info(info)
    return (
        "CLNSIG" in info_fields and "pathogenic" in info_fields["CLNSIG"].lower()
    ) or ("ONC" in info_fields and "oncogenic" in info_fields["ONC"].lower())


def matches_allowlist(variant, allowlist):
    """Check if a variant is in the allowlist (with '...' support)."""
    key = (variant["chrom"], variant["pos"], variant["ref"], variant["alt"])
    if key in allowlist:
        return True

    # Fallback for partial matches with '...'
    for allow in allowlist.values():
        if partial_match(allow["ref"], variant["ref"]) and partial_match(
            allow["alt"], variant["alt"]
        ):
            if variant["chrom"] == allow["chrom"] and variant["pos"] == allow["pos"]:
                return True
    return False


def allowlist_variants(variants, allowlist: dict = {}):
    """Annotate and optionally rescue variants using allowlist info."""
    for key, v in variants.items():
        allowlist_reasons = []

        if onc_clnsig(v["info"]):
            allowlist_reasons.append("ClinvarPathogenicOncogenic")

        if allowlist and matches_allowlist(v, allowlist):
            allowlist_reasons.append("ClinicalList")

        if allowlist_reasons:
            v["info"] += f";AllowlistStatus={','.join(allowlist_reasons)}"
            if v["filter"] != "PASS":
                v["info"] += f";AllowlistedFilters={v['filter']}"
                v["filter"] = "PASS"

    return variants


def write_info_headers(out_vcf):
    """Write INFO field definitions to the VCF header."""
    for info_id, desc in INFO_HEADERS.items():
        out_vcf.write(
            f'##INFO=<ID={info_id},Number=1,Type=String,Description="{desc}">\n'
        )


def write_vcf(output_path, input_path, variants):
    """Write updated VCF with allowlist annotations."""
    with gzip.open(input_path, "rt") as input_vcf, open(output_path, "w") as out_vcf:
        for line in input_vcf:
            if line.startswith("#CHROM"):
                write_info_headers(out_vcf)
                out_vcf.write(line)
                break
            else:
                out_vcf.write(line)

        for v in variants.values():
            line_fields = [str(v[field]) for field in VCF_FIELDS] + [
                str(s) for s in v["samples"]
            ]
            out_vcf.write("\t".join(line_fields) + "\n")


def read_variants(vcf: str) -> Dict[Tuple[str, str, str, str], Dict[str, str]]:
    """Parse VCF into a dictionary of variant entries."""
    variants = {}
    with gzip.open(vcf, "rt") as gz_vcf:
        for line in gz_vcf:
            if line.startswith("#"):
                continue
            fields: List[str] = line.strip().split("\t")
            variant_data = dict(zip(VCF_FIELDS, fields[:9]))
            variant_data["samples"] = fields[9:]
            key = (
                variant_data["chrom"],
                variant_data["pos"],
                variant_data["ref"],
                variant_data["alt"],
            )
            variants[key] = variant_data
    return variants


def read_allowlist(csv_path):
    """Load allowlist CSV into dictionary keyed by (chrom, pos, ref, alt)."""
    allowlist = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row["chrom"], row["pos"], row["ref"], row["alt"])
            allowlist[key] = row
    return allowlist


# --- CLI Entry Point ---
@click.command()
@click.argument("vcf_path", type=click.Path(exists=True))
@click.option(
    "--allowlist-path",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to allowlist file.",
)
@click.argument("output_file", type=click.Path())
def main(vcf_path: str, allowlist_path: str, output_file: str) -> None:
    """VCF processor with allowlist annotation support."""
    variants = read_variants(vcf_path)

    if allowlist_path:
        allowlist = read_allowlist(allowlist_path)
        variants = allowlist_variants(variants, allowlist)
    else:
        variants = allowlist_variants(variants)

    write_vcf(output_file, vcf_path, variants)


if __name__ == "__main__":
    main()
