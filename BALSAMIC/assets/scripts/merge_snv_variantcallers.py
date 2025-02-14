#!/usr/bin/env python
from typing import Dict, List, Tuple, Optional
from datetime import datetime
import click
import gzip
import sys
import os
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


def add_header_categories(
    header_categories: Dict[str, List[str]], header: List[str]
) -> Dict[str, List[str]]:
    """
    Categorizes VCF headers by type (e.g., FILTER, INFO).

    Parameters:
    - header_categories (Dict[str, List[str]]): Existing header categories.
    - header (List[str]): List of header lines.

    Returns:
    - Dict[str, List[str]]: Updated header categories.
    """
    for hdr_row in header:
        # Example hdr_row:
        # ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
        hdr_row_split: List[str] = hdr_row.split("=")
        # Example hdr_row_split
        # ["##FORMAT", "<ID", "AD,Number", "R,Type", "Integer,Description", "Allelic depths for the ref and alt alleles in the order listed">
        category: str = hdr_row_split[0].strip("#")
        # Example category:
        # FORMAT
        if category not in header_categories:
            header_categories[category] = []
            header_categories[category].append("=".join(hdr_row_split[1:]))
        else:
            header_categories[category].append("=".join(hdr_row_split[1:]))
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


class VCFHeaderMergeError(Exception):
    """Exception raised for errors in merging VCF headers."""

    pass


def merge_headers(vcf1: str, vcf2: str) -> List[str]:
    """Merges headers from two VCF files, ensuring no duplicate or conflicting entries."""

    def extract_line_id(line: str) -> str:
        """Extracts the ID from a VCF header line."""
        try:
            return line.split("=")[1].split(",")[0]
        except IndexError:
            return ""

    def process_category_lines(lines: List[str], category: str) -> Dict[str, str]:
        """Processes header lines within a category, merging duplicate IDs."""
        cat_lines = {}
        for line in lines:
            line_id = extract_line_id(line)
            if not line_id:  # Add malformed lines directly
                cat_lines[line] = line
            else:
                cat_lines[line_id] = merge_header_row(
                    cat_lines.get(line_id, line), line
                )
        if category == "INFO":
            cat_lines[
                "AF_LIST"
            ] = '<ID=AF_LIST,Number=.,Type=Float,Description="Allele Frequency list from both variant callers of a merged variant, in positional argument order">'
            cat_lines[
                "DP_LIST"
            ] = '<ID=DP_LIST,Number=.,Type=Integer,Description="Total Depth list from both variant callers of a merged variant, in positional argument order">'
        return cat_lines

    vcf1_name, vcf2_name = map(os.path.basename, [vcf1, vcf2])
    formatted_date = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    header1, header2 = collect_header(vcf1), collect_header(vcf2)

    header_categories: Dict[str, List[str]] = add_header_categories({}, header1)
    header_categories: Dict[str, List[str]] = add_header_categories(
        header_categories, header2
    )

    # Validate variant headers
    variant_header = [
        key
        for key in header_categories
        if "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" in key
    ]
    if len(variant_header) > 1:
        raise VCFHeaderMergeError(
            f"Error: Variant headers in {vcf1_name} and {vcf2_name} do not match. Cannot merge."
        )

    if variant_header:
        header_categories.pop(variant_header[0], None)

    # Process and merge headers
    merged_header_dict = {
        category: process_category_lines(lines, category)
        for category, lines in header_categories.items()
    }

    # Assemble merged header
    merged_header = [
        f"##{category}={line}"
        for category, cat_lines in merged_header_dict.items()
        for line in cat_lines.values()
    ]
    merged_header.extend(
        [
            f"##merge_snv_variantcallers=merge_snv_variantcallers.py {vcf1_name} {vcf2_name} --output output_merged.vcf",
            f"##merge_snv_variantcallers_processing_time={formatted_date}",
            f"##merge_snv_variantcallers=values in merged INFO fields are listed in the order of the input files: first from {vcf1_name}, then from {vcf2_name}",
        ]
    )
    if variant_header:
        merged_header.append("#" + variant_header[0])

    return [line.strip("\n") for line in merged_header]


def merge_header_row(line1: str, line2: str) -> str:
    """
    Merges two header rows with the same ID, resolving conflicts in Description and Number fields.

    Parameters:
    - line1 (str): First header row.
    - line2 (str): Second header row.

    Returns:
    - str: Merged header row.
    """
    updated_line: str = line1

    # Merge Description field for two header_rows with the same category and ID
    desc1: str = get_description_field(line1)
    desc2: str = get_description_field(line2)
    if desc1 and desc2 and desc1 != desc2:
        merged_desc: str = f"vcf1: {desc1} | vcf2: {desc2}"
        updated_line: str = update_description_text(updated_line, merged_desc)

    num_field1: str = get_number_field(line1)
    num_field2: str = get_number_field(line2)
    if num_field1 and num_field2 and num_field1 != num_field2:
        # If the number field is different for the same ID, set the num_field to "."
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


def merge_info_fields(info_fields: List[str]) -> str:
    """
    Merges multiple INFO fields into a single field, preserving all unique key-value pairs.

    Parameters:
    - info_fields (List[str]): List of INFO fields, where each field is a string.

    Returns:
    - str: Merged INFO field as a semicolon-separated string.
    """
    unique_fields = {"AF", "DP"}  # Set for faster lookups
    merged_info = {}

    def parse_info_field(field: str) -> Tuple[str, Optional[str]]:
        """Parses an individual INFO field into a key-value pair."""
        key, sep, value = field.partition("=")
        return key, value if sep else None

    def add_to_merged_info(key: str, value: Optional[str]) -> None:
        """Handles merging logic for INFO fields."""
        if value is None:
            merged_info[key] = None
        elif key in merged_info:
            merged_info[key] += f",{value}"
        else:
            merged_info[key] = value

    # Parse and merge INFO fields
    for field in info_fields:
        key, value = parse_info_field(field)
        add_to_merged_info(key, value)

    def process_unique_fields() -> None:
        """Ensures single-value fields retain only the first occurrence."""
        for key in unique_fields & merged_info.keys():
            values = merged_info[key].split(",")
            if len(values) > 1:
                merged_info[f"{key}_LIST"] = merged_info[key]
                merged_info[key] = values[0]

    process_unique_fields()

    # Construct the final INFO string
    return ";".join(
        f"{key}={value}" if value is not None else key
        for key, value in merged_info.items()
    )


def merge_filters(filter1: str, filter2: str) -> str:
    """
    Merges FILTER fields from two variants according to the following rules:
    - If both filters are "PASS", the result is "PASS".
    - If one or both filters contain non-"PASS" values, merge those values into a semicolon-separated string.
    - Remove "PASS" if it is present in either filter when merging non-"PASS" values.

    Parameters:
    - filter1 (str): The FILTER field of the first variant.
    - filter2 (str): The FILTER field of the second variant.

    Returns:
    - str: Merged FILTER field.
    """
    # Handle missing or placeholder values (".")
    filters1 = set(filter1.split(";")) if filter1 != "." else set()
    filters2 = set(filter2.split(";")) if filter2 != "." else set()

    # Combine filters
    merged_filters = filters1.union(filters2)

    # Remove "PASS" if there are any other filters
    if "PASS" in merged_filters and len(merged_filters) > 1:
        merged_filters.discard("PASS")

    # If the result is empty, default to "."
    return ";".join(sorted(merged_filters)) if merged_filters else "."


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
    variants1: Dict[Tuple[str, str, str, str], Dict[str, str]] = read_variants(vcf1)
    variants2: Dict[Tuple[str, str, str, str], Dict[str, str]] = read_variants(vcf2)

    # Combine keys and sort by genomic order (chromosome, position)
    all_keys: List[Tuple[str, str, str, str]] = sorted(
        set(variants1.keys()).union(variants2.keys()),
        key=lambda x: (parse_chromosome(x[0]), int(x[1])),
    )

    merged_variants = []

    for key in all_keys:
        variant1: Dict[str, str] = variants1.get(key)
        variant2: Dict[str, str] = variants2.get(key)

        if variant1 and variant2:
            # Variant exists in both VCFs; merge INFO fields

            id_ = variant1["id"] if variant1["id"] != "." else variant2["id"]

            # Merge info fields, creating a list of values for shared IDs (AD=123,124)
            merged_info: str = merge_info_fields(
                variant1["info"].split(";") + variant2["info"].split(";")
            )
            qual = variant1["qual"] if variant1["qual"] != "." else variant2["qual"]

            # Merge filter columns
            filter_ = merge_filters(variant1["filter"], variant2["filter"])

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
    merged_header: List[str] = merge_headers(vcf1, vcf2)

    # Merge variants based on the same chrom, pos, ref, alt, and merge INFO fields
    merged_variants: List[str] = merge_variants(vcf1, vcf2)

    # Write the merged VCF to output file
    write_merged_vcf(output_file, merged_header, merged_variants)


if __name__ == "__main__":
    main()
