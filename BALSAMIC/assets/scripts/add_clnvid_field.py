import pysam
import click


def add_clnvid_header(output_handle) -> None:
    """
    Writes the INFO header line for the CLNVID field to the output VCF.
    """
    header_line = '##INFO=<ID=CLNVID,Number=1,Type=Integer,Description="ClinVar Variation ID">'
    output_handle.write(f"{header_line}\n".encode("utf-8"))


def process_vcf(input_path: str, output_path: str) -> None:
    """
    Processes a bgzipped VCF file using pysam, adds the CLNVID INFO field based on the ID column,
    and writes to a new bgzipped VCF file in a tabix-compatible format.
    """
    with pysam.BGZFile(input_path, "r") as infile, pysam.BGZFile(
        output_path, "w"
    ) as outfile:
        for raw_line in infile:
            line = raw_line.decode("utf-8")
            if line.startswith("##"):
                outfile.write(raw_line)
            elif line.startswith("#"):
                add_clnvid_header(outfile)
                outfile.write(raw_line)
            else:
                fields = line.strip().split("\t")
                vcf_id = fields[2]
                if fields[7] == ".":
                    fields[7] = f"CLNVID={vcf_id}"
                else:
                    fields[7] += f";CLNVID={vcf_id}"
                modified_line = "\t".join(fields) + "\n"
                outfile.write(modified_line.encode("utf-8"))


@click.command()
@click.argument(
    "input_path", type=click.Path(exists=True, readable=True, dir_okay=False)
)
@click.argument("output_path", type=click.Path(writable=True, dir_okay=False))
def main(input_path: str, output_path: str) -> None:
    """
    Adds a CLNVID INFO field to each record in a bgzipped VCF file based on the ID column.

    INPUT_PATH: Path to the input VCF file (.vcf.gz).
    OUTPUT_PATH: Path to the output VCF file (.vcf.gz).
    """
    process_vcf(input_path, output_path)


if __name__ == "__main__":
    main()
