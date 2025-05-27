import click
from typing import TextIO


def add_clnvid_header(output_handle: TextIO) -> None:
    """
    Writes the INFO header line for the CLNVID field to the output VCF.

    Parameters:
        output_handle (TextIO): Writable handle to the output VCF.
    """
    header_line = (
        '##INFO=<ID=CLNVID,Number=1,Type=Integer,Description="ClinVar Variation ID">'
    )
    output_handle.write(f"{header_line}\n")


def process_vcf(input_path: str, output_path: str) -> None:
    """
    Processes a plain text VCF file, adds the CLNVID INFO field based on the ID column,
    and writes to a new plain text VCF file.

    Parameters:
        input_path (str): Path to the input VCF file.
        output_path (str): Path to the output VCF file.
    """
    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        for line in infile:
            if line.startswith("##"):
                outfile.write(line)
            elif line.startswith("#"):
                add_clnvid_header(outfile)
                outfile.write(line)
            else:
                fields = line.strip().split("\t")
                vcf_id = fields[2]
                if fields[7] == ".":
                    fields[7] = f"CLNVID={vcf_id}"
                else:
                    fields[7] += f";CLNVID={vcf_id}"
                outfile.write("\t".join(fields) + "\n")


@click.command()
@click.argument(
    "input_path", type=click.Path(exists=True, readable=True, dir_okay=False)
)
@click.argument("output_path", type=click.Path(writable=True, dir_okay=False))
def main(input_path: str, output_path: str) -> None:
    """
    Adds a CLNVID INFO field to each record in a VCF file based on the ID column.

    INPUT_PATH: Path to the input (non-gzipped) VCF file.
    OUTPUT_PATH: Path to the output (non-gzipped) VCF file.
    """
    process_vcf(input_path, output_path)


if __name__ == "__main__":
    main()
