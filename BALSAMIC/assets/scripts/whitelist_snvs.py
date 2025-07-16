#!/usr/bin/env python
from typing import Dict, List, Tuple
import click
import gzip
import csv

def onc_clnsig(info: str) -> bool:
    for field in info.split(";"):
        if "=" not in field:
            continue
        key, value = field.split("=", 1)
        if key == "CLNSIG" and "pathogenic" in value.lower():
            return True
        if key == "ONC" and "oncogenic" in value.lower():
            return True
    return False

def matches_whitelist(variant, white):
    """Check if a variant matches a whitelist entry (with '...' support)."""
    if variant["chrom"] != white["chrom"] or variant["pos"] != white["pos"]:
        return False

    ref = variant["ref"]
    alt = variant["alt"]

    def partial_match(query, target):
        return query.strip("...") in target if "..." in query else query == target

    return partial_match(white["ref"], ref) and partial_match(white["alt"], alt)

def whitelist_variants(variants, whitelist):
    for key, v in variants.items():
        whitelist_reasons = []

        # Check CLNSIG / ONC evidence
        if onc_clnsig(v["info"]):
            whitelist_reasons.append("ClinvarPathogenicOncogenic")

        # Check against manual whitelist
        for white in whitelist.values():
            if matches_whitelist(v, white):
                whitelist_reasons.append("ClinicalList")
                break  # One match is sufficient

        if whitelist_reasons:
            v["info"] += f";WhitelistStatus={','.join(whitelist_reasons)}"

            if v["filter"] != "PASS":
                v["info"] += f";WhitelistedFilters={v['filter']}"
                v["filter"] = "PASS"
    return variants

def write_vcf(output_path, input_path, variants):
    with gzip.open(input_path, "rt") as input_vcf, open(output_path, "w") as out_vcf:
        for line in input_vcf:
            if line.startswith("#CHROM"):
                # Add new INFO field descriptions before column header line
                out_vcf.write('##INFO=<ID=WhitelistStatus,Number=1,Type=String,Description="Variant whitelisted based on CLNSIG/ONC or manually curated clinical list">\n')
                out_vcf.write('##INFO=<ID=WhitelistedFilters,Number=1,Type=String,Description="Original filters for whitelisted variants that were overridden">\n')
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
    """
    Reads variants from a VCF file, storing them in a dictionary keyed by variant identifiers.

    Parameters:
    - vcf (str): Path to the VCF file.

    Returns:
    - Dict[Tuple[str, str, str, str], Dict[str, str]]: Dictionary of variant data.
    """
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
    whitelist = {}
    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            whitelist[int(row["id"])] = row
    return whitelist


@click.command()
@click.argument("vcf_path", type=click.Path(exists=True))
@click.argument("whitelist_path", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
def main(vcf_path: str, whitelist_path: str, output_file: str) -> None:
    """
    Command-line interface to merge two VCF files with overlapping variants and headers,
    supporting bgzip format.

    Parameters:
    - vcf (str): Path to the first VCF file.
    - whitelist (str): Path to whitelist file.
    - output_file (str): Path to the output merged VCF file.
    """

    variants = read_variants(vcf_path)
    whitelist = read_whitelist(whitelist_path)


    whitelist_variants(variants, whitelist)
    write_vcf(output_file, vcf_path, variants)

if __name__ == "__main__":
    main()
