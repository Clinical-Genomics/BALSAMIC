#!/usr/bin/env python
import vcfpy
import gzip
import click
from typing import Dict, Tuple, Any


def merge_headers(reader1: vcfpy.Reader, reader2: vcfpy.Reader) -> vcfpy.Header:
    """
    Combines all header lines from two VCF files, adding unique lines from each.

    Parameters:
    - reader1 (vcfpy.Reader): First VCF reader.
    - reader2 (vcfpy.Reader): Second VCF reader.

    Returns:
    - vcfpy.Header: Merged header containing unique lines from both readers.
    """
    merged_header = reader1.header.copy()

    # Define the types of header lines we want to merge
    header_line_types = ["INFO", "FILTER", "FORMAT", "ALT"]

    # Loop over each line type and merge lines from reader2 that are not in merged_header
    for line_type in header_line_types:
        # Check if the appropriate method for each header line type is available
        method_name = f"{line_type.lower()}_ids"
        if hasattr(merged_header, method_name):
            # Get the list of IDs for the header line type
            line_ids = getattr(merged_header, method_name)()
        else:
            # Handle cases where the specific method is not available
            line_ids = [line.id for line in merged_header.get_lines(line_type)]

        for line in reader2.header.get_lines(line_type):
            if line.id not in line_ids:
                merged_header.add_line(line)

    return merged_header


def merge_vcf_records(record1: vcfpy.Record, record2: vcfpy.Record) -> Dict[str, Any]:
    """
    Merges INFO fields of two matching VCF records.

    Parameters:
    - record1 (vcfpy.Record): First VCF record.
    - record2 (vcfpy.Record): Second VCF record.

    Returns:
    - Dict[str, Any]: Dictionary of merged INFO fields.
    """
    merged_info = record1.INFO.copy()

    for key, value in record2.INFO.items():
        if key in merged_info:
            if isinstance(merged_info[key], list):
                if isinstance(value, list):
                    merged_info[key].extend(
                        [v for v in value if v not in merged_info[key]]
                    )
                else:
                    merged_info[key].append(value)
            else:
                merged_info[key] = [merged_info[key]]
                if isinstance(value, list):
                    merged_info[key].extend(value)
                else:
                    merged_info[key].append(value)
        else:
            merged_info[key] = value

    return merged_info


def process_vcf_records(
    writer: vcfpy.Writer, reader1: vcfpy.Reader, reader2: vcfpy.Reader
) -> None:
    """
    Processes and writes merged VCF records from two VCF readers.

    Parameters:
    - writer (vcfpy.Writer): Writer to output the merged records.
    - reader1 (vcfpy.Reader): First VCF reader.
    - reader2 (vcfpy.Reader): Second VCF reader.
    """

    def get_variant_key(record: vcfpy.Record) -> Tuple[str, int, str, Tuple[str, ...]]:
        """Creates a unique key for each variant based on CHROM, POS, REF, and ALT."""
        return (
            record.CHROM,
            record.POS,
            record.REF,
            tuple(str(alt) for alt in record.ALT),
        )

    variants1 = {get_variant_key(record): record for record in reader1}
    variants2 = {get_variant_key(record): record for record in reader2}

    processed = set()

    for key, record1 in variants1.items():
        if key in variants2:
            record2 = variants2[key]
            merged_info = merge_vcf_records(record1, record2)
            merged_record = vcfpy.Record(
                CHROM=record1.CHROM,
                POS=record1.POS,
                ID=record1.ID,
                REF=record1.REF,
                ALT=record1.ALT,
                QUAL=record1.QUAL,
                FILTER=record1.FILTER,
                INFO=merged_info,
                FORMAT=record1.FORMAT,
                calls=record1.calls,
            )
            writer.write_record(merged_record)
            processed.add(key)
        else:
            writer.write_record(record1)

    for key, record2 in variants2.items():
        if key not in processed:
            writer.write_record(record2)

    writer.close()


def merge_vcf_files(vcf1: str, vcf2: str, output_file: str) -> None:
    """
    Merges two VCF files, combining headers and overlapping records.

    Parameters:
    - vcf1 (str): Path to the first VCF file.
    - vcf2 (str): Path to the second VCF file.
    - output_file (str): Path to the output merged VCF file.
    """
    reader1 = vcfpy.Reader.from_path(vcf1)
    reader2 = vcfpy.Reader.from_path(vcf2)
    merged_header = merge_headers(reader1, reader2)

    with (
        gzip.open(output_file, "wt")
        if output_file.endswith(".gz")
        else open(output_file, "w")
    ) as f:
        writer = vcfpy.Writer.from_stream(f, merged_header)
        process_vcf_records(writer, reader1, reader2)


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
    merge_vcf_files(vcf1, vcf2, output_file)
    print(f"Merged VCF file created at: {output_file}")


if __name__ == "__main__":
    main()
