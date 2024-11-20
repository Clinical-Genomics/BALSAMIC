#!/usr/bin/env python
import gzip
import sys
import re
import click

def collect_header(vcf):
    header = []
    with gzip.open(vcf, 'rt') as gz_vcf:
        for line in gz_vcf:
            if line.startswith("#"):
                header.append(line)
            else:
                break
    return header

def add_header_categories(header_categories, header):
    for hdr_row in header:
        hdr_row_split = hdr_row.split("=")
        cat = hdr_row_split[0].strip("#")
        if cat not in header_categories:
            header_categories[cat] = []
            header_categories[cat].append("=".join(hdr_row_split[1:]))
        else:
            header_categories[cat].append("=".join(hdr_row_split[1:]))
    return header_categories

def get_description_field(header_row):
    pattern = r'Description="(.*?)"'
    match = re.search(pattern, header_row)
    if match:
        return match.group(1)
    return None

def update_description_text(row, new_description):
    return re.sub(r'Description="(.*?)"', f'Description="{new_description}"', row)

def get_number_field(header_row):
    pattern = r'Number=([^,]+)'
    match = re.search(pattern, header_row)
    if match:
        return match.group(1)
    return None

def update_number_text(row, new_number):
    return re.sub(r'Number=([^,]+)', f'Number="{new_number}"', row)

def merge_headers(vcf1, vcf2):
    """
    Combines all header lines from two VCF files, adding unique lines from each.
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
    merged_header.append(variant_header[0])
    merged_header = [line.strip("\n") for line in merged_header]
    return merged_header


def merge_header_row(line1, line2):
    updated_line = line1

    # If Descriptions are different, merge them
    description1 = get_description_field(line1)
    description2 = get_description_field(line2)
    if description1 and description2 and description1 != description2:
        # Merge Descriptions
        merged_description = f"vcf1: {description1} | vcf2: {description2}"
        updated_line = update_description_text(line1, merged_description)

    # If Number is different, set to "."
    number1 = get_number_field(line1)
    number2 = get_number_field(line2)
    if number1 and number2 and number1 != number2:
        # Conflicting Number type, setting to "."
        # Print warning to standard error
        print(
            f"Warning: Number fields differ for header lines. line1: {line1}, line2: {line2}",
            file=sys.stderr
        )
        updated_line = update_number_text(updated_line, ".")

    return updated_line


def read_variants(vcf):
    variants = {}

    with gzip.open(vcf, 'rt') as gz_vcf:
        for line in gz_vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            fields = line.strip().split("\t")

            normal = False
            if len(fields) > 10:  # TUMOR + NORMAL
                chrom, pos, id, ref, alt, qual, filter, info, format, tumor, normal = fields[:11]
            else:  # TUMOR ONLY
                chrom, pos, id, ref, alt, qual, filter, info, format, tumor = fields[:11]

                # Use a tuple of (chrom, pos, ref, alt) as the key
            key = (chrom, pos, ref, alt)

            # Store the variant along with its ID, QUAL, and FILTER, and INFO field
            variants[key] = {}
            variants[key]["chrom"] = chrom
            variants[key]["pos"] = pos
            variants[key]["id"] = id
            variants[key]["ref"] = ref
            variants[key]["alt"] = alt
            variants[key]["qual"] = qual
            variants[key]["filter"] = filter
            variants[key]["info"] = info
            variants[key]["format"] = format
            variants[key]["tumor"] = tumor
            if normal:
                variants[key]["normal"] = normal

    return variants


def merge_info_fields(info_fields):
    """
    Merges INFO fields by combining all unique key-value pairs from a list of INFO fields.
    """
    merged_info = {}
    merged_info_single_fields = {}
    for info in info_fields:
        info_key_value = info.split("=")
        if len(info_key_value) < 2:
            key = info_key_value[0]
            merged_info_single_fields[key] = key
            continue

        key = info_key_value[0]
        value = info_key_value[1]

        if key in merged_info:
            existing_val = merged_info[key]
            merged_info[key] = ",".join([existing_val, value])
        else:
            merged_info[key] = value

    # Recreate INFO field with merged key-value pairs
    merged_info_str = ";".join([f"{key}" for key in merged_info_single_fields])
    merged_info_str = merged_info_str + ";".join([f"{key}={value}" for key, value in merged_info.items()])
    return merged_info_str


def merge_variants(vcf1, vcf2):
    """
    Merges variants from two VCF files based on the same chromosome, position, and alteration.
    INFO fields are merged for variants with the same keys. ID, QUAL, and FILTER values
    from vcf1 are prioritized if variants overlap.
    """
    variants_vcf1 = read_variants(vcf1)
    variants_vcf2 = read_variants(vcf2)

    # Now, merge the variants with the same key
    merged_variants = []
    for variant in variants_vcf1:
        v_dict = variants_vcf1[variant]
        chrom = v_dict["chrom"]
        pos = v_dict["pos"]
        id = v_dict["id"]
        ref = v_dict["ref"]
        alt = v_dict["alt"]
        qual = v_dict["qual"]
        filter = v_dict["filter"]
        info_1 = v_dict["info"].split(";")
        format = v_dict["format"]
        tumor = v_dict["tumor"]
        if "normal" in v_dict:
            normal = v_dict["normal"]

        if variant in variants_vcf2:
            # Choose ID and merge INFO field
            id_2 = variants_vcf2[variant]["id"]
            info_2 = variants_vcf2[variant]["info"].split(";")

        if id != id_2 and id == ".":  # Switch ID to ID from VCF2
            id = id_2

        info_fields = []
        for info in info_1:
            info_fields.append(info)
        for info in info_2:
            info_fields.append(info)

        merged_info = merge_info_fields(info_fields)

        # Combine the variant fields (ID, QUAL, FILTER prioritized from vcf1)
        if "normal" in variants_vcf1[variant]:
            merged_variant = f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{merged_info}\t{format}\t{tumor}\t{normal}"
        else:
            merged_variant = f"{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{merged_info}\t{format}\t{tumor}"
        merged_variants.append(merged_variant)

    return merged_variants

def write_merged_vcf(output_file, header, merged_variants):
    """
    Writes the merged variants along with the header to the output VCF file.
    """
    with (
        gzip.open(output_file, "wt")
        if output_file.endswith(".gz")
        else open(output_file, "w")
    ) as f:
        # Write the header
        for line in header:
            f.write(line + "\n")
        # Write merged variants
        for variant in merged_variants:
            f.write(variant + "\n")

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