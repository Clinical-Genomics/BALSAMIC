#!/usr/bin/env python
import vcfpy
import gzip
import click


def merge_headers(reader1, reader2):
    """Combines all header lines from two VCF files, adding unique lines from each."""
    merged_header = reader1.header.copy()

    # Define the types of header lines we want to merge
    header_line_types = ["INFO", "FILTER", "FORMAT", "ALT"]

    # Loop over each line type and merge lines from reader2 that are not in merged_header
    for line_type in header_line_types:
        for line in reader2.header.get_lines(line_type):
            if line.id not in getattr(merged_header, f"{line_type.lower()}_ids")():
                merged_header.add_line(line)

    return merged_header


def merge_vcf_records(record1, record2):
    """Merges INFO fields of two matching VCF records."""
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


def merge_vcf_files(file1, file2, output_file):
    reader1 = vcfpy.Reader.from_path(file1)
    reader2 = vcfpy.Reader.from_path(file2)

    # Merge headers from both VCF files
    merged_header = merge_headers(reader1, reader2)

    # Open the output file in gzip mode if it has a .gz extension
    if output_file.endswith(".gz"):
        with gzip.open(output_file, "wt") as f:
            writer = vcfpy.Writer.from_stream(f, merged_header)
            process_vcf_records(writer, reader1, reader2)
    else:
        with open(output_file, "w") as f:
            writer = vcfpy.Writer.from_stream(f, merged_header)
            process_vcf_records(writer, reader1, reader2)


def process_vcf_records(writer, reader1, reader2):
    """Processes and writes merged VCF records."""
    variants1 = {
        (
            record.CHROM,
            record.POS,
            record.REF,
            tuple(str(alt) for alt in record.ALT),
        ): record
        for record in reader1
    }
    variants2 = {
        (
            record.CHROM,
            record.POS,
            record.REF,
            tuple(str(alt) for alt in record.ALT),
        ): record
        for record in reader2
    }

    processed = set()

    # Merge overlapping records
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


@click.command()
@click.argument("file1", type=click.Path(exists=True))
@click.argument("file2", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
def main(file1, file2, output_file):
    """Merge two VCF files with overlapping variants and headers, supporting bgzip format."""
    merge_vcf_files(file1, file2, output_file)
    print(f"Merged VCF file created at: {output_file}")


if __name__ == "__main__":
    main()
