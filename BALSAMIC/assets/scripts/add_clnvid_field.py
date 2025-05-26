import sys
import gzip
from typing import TextIO


def add_clnvid_header(output_handle: TextIO) -> None:
    """
    Writes the INFO header line for the CLNVID field to the output VCF.

    Parameters:
        output_handle (TextIO): Writable handle to the output VCF.
    """
    header_line = '##INFO=<ID=CLNVID,Number=1,Type=Integer,Description="ClinVar Variation ID">'
    output_handle.write(f"{header_line}\n")


def process_vcf(input_path: str, output_path: str) -> None:
    """
    Processes a gzipped VCF file, adds the CLNVID INFO field based on the ID column,
    and writes to a new gzipped VCF file.

    Parameters:
        input_path (str): Path to the input VCF file (gzip compressed).
        output_path (str): Path to the output VCF file (gzip compressed).
    """
    with gzip.open(input_path, 'rt') as infile, gzip.open(output_path, 'wt') as outfile:
        for line in infile:
            if line.startswith('##'):
                outfile.write(line)
            elif line.startswith('#'):
                add_clnvid_header(outfile)
                outfile.write(line)
            else:
                fields = line.strip().split('\t')
                vcf_id = fields[2]
                if fields[7] == ".":
                    fields[7] = f"CLNVID={vcf_id}"
                else:
                    fields[7] += f";CLNVID={vcf_id}"
                outfile.write('\t'.join(fields) + '\n')


def main() -> None:
    """
    Main entry point for the script. Expects two command-line arguments:
    input path and output path.
    """
    if len(sys.argv) != 3:
        print("Usage: add_clnvid_field.py <input.vcf.gz> <output.vcf.gz>", file=sys.stderr)
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    process_vcf(input_vcf, output_vcf)