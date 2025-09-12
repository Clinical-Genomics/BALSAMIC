import pysam
import click


def add_clnvid_header(output_handle) -> None:
    """
    Writes the INFO header line for the CLNVID field to the output VCF.
    """
    # Use String to be safe; ID column can contain non-numeric identifiers.
    header_line = '##INFO=<ID=CLNVID,Number=1,Type=String,Description="ClinVar Variation ID taken from the VCF ID column">'
    output_handle.write(f"{header_line}\n".encode("utf-8"))


def process_vcf(input_path: str, output_path: str) -> None:
    """
    Processes a bgzipped VCF file using pysam, adds the CLNVID INFO field based on the ID column,
    and writes to a new bgzipped VCF file in a tabix-compatible format.
    """
    saw_clnvid_header = False

    with pysam.BGZFile(input_path, "r") as infile, pysam.BGZFile(
        output_path, "w"
    ) as outfile:
        for raw_line in infile:
            # Pass through meta headers; track if CLNVID header already exists
            if raw_line.startswith(b"##"):
                if b"##INFO=<ID=CLNVID" in raw_line:
                    saw_clnvid_header = True
                outfile.write(raw_line)
                continue

            # Column header line: add CLNVID header if missing, then write
            if raw_line.startswith(b"#"):
                if not saw_clnvid_header:
                    add_clnvid_header(outfile)
                    saw_clnvid_header = True
                outfile.write(raw_line)
                continue

            # Variant line
            line = raw_line.decode("utf-8").rstrip("\n")
            fields = line.split("\t")

            # Ensure we have at least up to INFO column
            if len(fields) < 8:
                # Malformed line, write back unchanged
                outfile.write((line + "\n").encode("utf-8"))
                continue

            vcf_id = fields[2]
            info = fields[7]

            # Only add CLNVID when ID column is not '.'
            if vcf_id != ".":
                if info == "." or info == "":
                    info = f"CLNVID={vcf_id}"
                elif "CLNVID=" not in info:
                    info = f"{info};CLNVID={vcf_id}"

            # Write updated INFO back to fields
            fields[7] = info
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
