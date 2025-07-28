#!/usr/bin/env python
from typing import Dict, List, Tuple
import click
import gzip
import csv

# --- Constants for INFO headers ---
INFO_HEADERS = {
    "WhitelistStatus": "Variant whitelisted based on CLNSIG/ONC or manually curated clinical list",
    "WhitelistedFilters": "Original filters for whitelisted variants that were overridden",
}


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


def matches_whitelist(variant, white):
    """Check if a variant matches a whitelist entry (with '...' support)."""
    if variant["chrom"] != white["chrom"] or variant["pos"] != white["pos"]:
        return False
    return partial_match(white["ref"], variant["ref"]) and partial_match(
        white["alt"], variant["alt"]
    )


def whitelist_variants(variants, whitelist: dict = {}):
    """Annotate and optionally rescue variants using whitelist info."""
    for key, v in variants.items():
        whitelist_reasons = []

        if onc_clnsig(v["info"]):
            whitelist_reasons.append("ClinvarPathogenicOncogenic")

        for white in whitelist.values():
            if matches_whitelist(v, white):
                whitelist_reasons.append("ClinicalList")
                break

        if whitelist_reasons:
            v["info"] += f";WhitelistStatus={','.join(whitelist_reasons)}"
            if v["filter"] != "PASS":
                v["info"] += f";WhitelistedFilters={v['filter']}"
                v["filter"] = "PASS"

    return variants


def write_info_headers(out_vcf):
    """Write INFO field definitions to the VCF header."""
    for info_id, desc in INFO_HEADERS.items():
        out_vcf.write(
            f'##INFO=<ID={info_id},Number=1,Type=String,Description="{desc}">\n'
        )


def write_vcf(output_path, input_path, variants):
    """Write updated VCF with whitelist annotations."""
    with gzip.open(input_path, "rt") as input_vcf, open(output_path, "w") as out_vcf:
        for line in input_vcf:
            if line.startswith("#CHROM"):
                write_info_headers(out_vcf)
                out_vcf.write(line)
                break
            else:
                out_vcf.write(line)

        for key, v in variants.items():
            fields = [
                v["chrom"],
                v["pos"],
                v["id"],
                v["ref"],
                v["alt"],
                v["qual"],
                v["filter"],
                v["info"],
                v["format"],
                *v["samples"],
            ]
            out_vcf.write("\t".join(fields) + "\n")


def read_variants(vcf: str) -> Dict[Tuple[str, str, str, str], Dict[str, str]]:
    """Parse VCF into a dictionary of variant entries."""
    variants = {}
    with gzip.open(vcf, "rt") as gz_vcf:
        for line in gz_vcf:
            if line.startswith("#"):
                continue
            fields: List[str] = line.strip().split("\t")
            chrom, pos, v_id, ref, alt, qual, filter_, info, format_, *samples = fields
            key = (chrom, pos, ref, alt)
            variants[key] = {
                "chrom": chrom,
                "pos": pos,
                "id": v_id,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": filter_,
                "info": info,
                "format": format_,
                "samples": samples,
            }
    return variants


def read_whitelist(csv_path):
    """Load whitelist CSV into dictionary indexed by ID."""
    whitelist = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            whitelist[int(row["id"])] = row
    return whitelist


# --- CLI Entry Point ---
@click.command()
@click.argument("vcf_path", type=click.Path(exists=True))
@click.option(
    "--whitelist-path",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to whitelist file.",
)
@click.argument("output_file", type=click.Path())
def main(vcf_path: str, whitelist_path: str, output_file: str) -> None:
    """VCF processor with whitelist annotation support."""
    variants = read_variants(vcf_path)

    if whitelist_path:
        whitelist = read_whitelist(whitelist_path)
        variants = whitelist_variants(variants, whitelist)
    else:
        variants = whitelist_variants(variants)

    write_vcf(output_file, vcf_path, variants)


if __name__ == "__main__":
    main()
