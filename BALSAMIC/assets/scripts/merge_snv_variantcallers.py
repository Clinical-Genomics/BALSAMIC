#!/usr/bin/env python
from typing import Dict, List, Tuple, Optional
import click
import gzip
import sys
import re


def collect_header(vcf: str) -> List[str]:
    """
    Collects the header lines from a VCF file.

    Parameters:
    - vcf (str): Path to the VCF file.

    Returns:
    - List[str]: List of header lines.
    """
    header = []
    with gzip.open(vcf, "rt") as gz_vcf:
        for line in gz_vcf:
            if line.startswith("#"):
                header.append(line.strip())
            else:
                break
    return header


def add_header_categories(header_categories, header):
    """
    Categorizes VCF headers by type (e.g., FILTER, INFO).

    Parameters:
    - header_categories (Dict[str, List[str]]): Existing header categories.
    - header (List[str]): List of header lines.

    Returns:
    - Dict[str, List[str]]: Updated header categories.
    """
    for hdr_row in header:
        hdr_row_split = hdr_row.split("=")
        cat = hdr_row_split[0].strip("#")
        if cat not in header_categories:
            header_categories[cat] = []
            header_categories[cat].append("=".join(hdr_row_split[1:]))
        else:
            header_categories[cat].append("=".join(hdr_row_split[1:]))
    return header_categories


def get_description_field(header_row: str) -> Optional[str]:
    """
    Extracts the Description field from a VCF header row.

    Parameters:
    - header_row (str): A header row.

    Returns:
    - Optional[str]: The description field, if present.
    """
    match = re.search(r'Description="(.*?)"', header_row)
    return match.group(1) if match else None


def update_description_text(row: str, new_description: str) -> str:
    """
    Updates the Description field in a VCF header row.

    Parameters:
    - row (str): A header row.
    - new_description (str): The new description text.

    Returns:
    - str: Updated header row.
    """
    return re.sub(r'Description="(.*?)"', f'Description="{new_description}"', row)


def get_number_field(header_row: str) -> Optional[str]:
    """
    Extracts the Number field from a VCF header row.

    Parameters:
    - header_row (str): A header row.

    Returns:
    - Optional[str]: The Number field, if present.
    """
    match = re.search(r"Number=([^,]+)", header_row)
    return match.group(1) if match else None


def update_number_text(row: str, new_number: str) -> str:
    """
    Updates the Number field in a VCF header row.

    Parameters:
    - row (str): A header row.
    - new_number (str): The new Number field value.

    Returns:
    - str: Updated header row.
    """
    return re.sub(r"Number=([^,]+)", f"Number={new_number}", row)


def merge_headers(vcf1, vcf2):
    """
    Merges headers from two VCF files, ensuring no duplicate or conflicting entries.

    Parameters:
    - vcf1 (str): Path to the first VCF file.
    - vcf2 (str): Path to the second VCF file.

    Returns:
    - List[str]: Merged header lines.
    """
    header1 = collect_header(vcf1)
    header2 = collect_header(vcf2)

    header_categories = {}
    header_categories = add_header_categories(header_categories, header1)
    header_categories = add_header_categories(header_categories, header2)

    # Handle variant header
    variant_header_string = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    variant_header = [key for key in header_categories if variant_header_string in key]

    count_unique_variant_headers = len(variant_header)

    if count_unique_variant_headers > 1:
        print("Error, variant headers in vcf1 and vcf2 do not match. Cannot merge.")
        sys.exit(1)

    for key in variant_header:
        header_categories.pop(key)

    # Merge identical rows
    merged_header_dict = {}
    for cat in header_categories:
        cat_lines = {}
        merged_header_dict[cat] = {}
        for line in header_categories[cat]:
            if line in cat_lines:
                # Line has already been added in this category, skip it
                continue
            else:
                try:
                    line_id = line.split("=")[1].split(",")[0]
                except:
                    # Line does not match standard VCF header format, just add it
                    cat_lines[line] = line
                    continue
            # Handle merging standard VCF header lines such as:
            # ##FILTER=<ID=balsamic_low_tumor_ad,Description="Set if not true: FORMAT/AD[0:1] >= 5.0">
            # ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
            if line_id not in cat_lines:
                cat_lines[line_id] = line
            else:
                # Merge lines with the same category and ID
                line1 = cat_lines[line_id]
                line2 = line
                cat_lines[line_id] = merge_header_row(line1, line2)
        merged_header_dict[cat] = cat_lines

    # Consolidate merged header into list to be written
    merged_header = []
    for cat in merged_header_dict:
        for key, line in merged_header_dict[cat].items():
            merged_header.append(f"##{cat}={line}")
    merged_header.append("#" + variant_header[0])
    merged_header = [line.strip("\n") for line in merged_header]
    return merged_header


def merge_header_row(line1: str, line2: str) -> str:
    """
    Merges two header rows with the same ID, resolving conflicts in Description and Number fields.

    Parameters:
    - line1 (str): First header row.
    - line2 (str): Second header row.

    Returns:
    - str: Merged header row.
    """
    updated_line = line1
    desc1 = get_description_field(line1)
    desc2 = get_description_field(line2)
    if desc1 and desc2 and desc1 != desc2:
        merged_desc = f"vcf1: {desc1} | vcf2: {desc2}"
        updated_line = update_description_text(updated_line, merged_desc)

    num1 = get_number_field(line1)
    num2 = get_number_field(line2)
    if num1 and num2 and num1 != num2:
        updated_line = update_number_text(updated_line, ".")
        print(
            f"Warning: Number fields differ. line1: {line1}, line2: {line2}",
            file=sys.stderr,
        )

    return updated_line


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
            fields = line.strip().split("\t")
            chrom, pos, id, ref, alt, qual, filter_, info, format_, *samples = fields
            key = (chrom, pos, ref, alt)
            variants[key] = {
                "chrom": chrom,
                "pos": pos,
                "id": id,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": filter_,
                "info": info,
                "format": format_,
                "samples": samples,
            }
    return variants


def merge_info_fields(info_fields: List[str]) -> str:
    """
    Merges multiple INFO fields into a single field, preserving all unique key-value pairs.

    Parameters:
    - info_fields (List[str]): List of INFO fields.

    Returns:
    - str: Merged INFO field.
    """
    merged_info = {}
    for field in info_fields:
        key, sep, value = field.partition("=")
        if sep:  # Key-value pair
            merged_info[key] = (
                merged_info.get(key, "") + ("," if key in merged_info else "") + value
            )
        else:  # Key only
            merged_info[key] = None
    return ";".join(f"{k}={v}" if v is not None else f"{k};" for k, v in merged_info.items())

def merge_variants(vcf1: str, vcf2: str) -> List[str]:
    """
    Merges variants from two VCF files, keeping all unique variants
    and merging INFO fields for shared variants, sorted in genomic order.

    Parameters:
    - vcf1 (str): Path to the first VCF file.
    - vcf2 (str): Path to the second VCF file.

    Returns:
    - List[str]: Merged variant lines, sorted in genomic order.
    """
    variants1 = read_variants(vcf1)
    variants2 = read_variants(vcf2)

    # Combine keys and sort by genomic order (chromosome, position)
    all_keys = sorted(
        set(variants1.keys()).union(variants2.keys()),
        key=lambda x: (parse_chromosome(x[0]), int(x[1])),
    )

    merged_variants = []

    for key in all_keys:
        variant1 = variants1.get(key)
        variant2 = variants2.get(key)

        if variant1 and variant2:
            # Variant exists in both VCFs; merge INFO fields
            id_ = variant1["id"] if variant1["id"] != "." else variant2["id"]
            merged_info = merge_info_fields(
                variant1["info"].split(";") + variant2["info"].split(";")
            )
            qual = variant1["qual"] if variant1["qual"] != "." else variant2["qual"]
            filter_ = variant1["filter"] if variant1["filter"] != "." else variant2["filter"]
            merged_variants.append(
                "\t".join(
                    [
                        variant1["chrom"],
                        variant1["pos"],
                        id_,
                        variant1["ref"],
                        variant1["alt"],
                        qual,
                        filter_,
                        merged_info,
                        variant1["format"],
                    ]
                    + variant1["samples"]
                )
            )
        elif variant1:
            # Unique to VCF1
            merged_variants.append(
                "\t".join(
                    [
                        variant1["chrom"],
                        variant1["pos"],
                        variant1["id"],
                        variant1["ref"],
                        variant1["alt"],
                        variant1["qual"],
                        variant1["filter"],
                        variant1["info"],
                        variant1["format"],
                    ]
                    + variant1["samples"]
                )
            )
        elif variant2:
            # Unique to VCF2
            merged_variants.append(
                "\t".join(
                    [
                        variant2["chrom"],
                        variant2["pos"],
                        variant2["id"],
                        variant2["ref"],
                        variant2["alt"],
                        variant2["qual"],
                        variant2["filter"],
                        variant2["info"],
                        variant2["format"],
                    ]
                    + variant2["samples"]
                )
            )

    return merged_variants


def parse_chromosome(chrom: str) -> Tuple[int, str]:
    """
    Parses chromosome names to ensure proper sorting.

    Parameters:
    - chrom (str): Chromosome name.

    Returns:
    - Tuple[int, str]: Parsed chromosome for sorting (numeric or special).
    """
    # Handle numeric chromosomes as integers for proper sorting
    if chrom.isdigit():
        return (int(chrom), "")
    # Handle special chromosomes like X, Y, M, etc.
    chrom_order = {"X": 23, "Y": 24, "M": 25}
    return (chrom_order.get(chrom.upper(), 99), chrom)


def write_merged_vcf(output_file: str, header: List[str], variants: List[str]) -> None:
    """
    Writes the merged VCF to a file.

    Parameters:
    - output_file (str): Path to the output file.
    - header (List[str]): Merged header lines.
    - variants (List[str]): Merged variant lines.
    """
    with open(output_file, "w") as f:
        f.write("\n".join(header) + "\n")
        f.write("\n".join(variants) + "\n")


@click.command()
@click.argument("vcf1", type=click.Path(exists=True))
@click.argument("vcf2", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
def main(vcf1: str, vcf2: str, output_file: str) -> None:
    """
    Command-line interface to merge two VCF files with overlapping variants and headers,
    supporting bgzip format.

    Parameters:
    - vcf1 (str): Path to the first VCF file.
    - vcf2 (str): Path to the second VCF file.
    - output_file (str): Path to the output merged VCF file.
    """
    # Merge header
    merged_header = merge_headers(vcf1, vcf2)

    # Merge variants based on the same chrom, pos, ref, alt, and merge INFO fields
    merged_variants = merge_variants(vcf1, vcf2)

    # Write the merged VCF to output file
    write_merged_vcf(output_file, merged_header, merged_variants)


if __name__ == "__main__":
    main()
