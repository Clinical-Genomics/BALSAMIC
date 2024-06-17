import click
import pandas as pd


@click.command()
@click.option(
    "-o",
    "--output-file",
    required=True,
    type=click.Path(exists=False),
    help="Name of output-file.",
)
@click.option(
    "-c",
    "--normalised-coverage-path",
    required=True,
    type=click.Path(exists=True),
    help="Input CNVkit cnr result.",
)
def create_gens_cov_file(output_file, normalised_coverage_path):
    """
    Post-processes the CNVkit cnr output for upload to GENS.
    Removing Antitarget regions and outputting the coverages in multiple resolution-formats.

    :param output_file: Path to GENS output.cov file
    :param normalised_coverage_path: Path to input CNVkit cnr file.
    """
    # Process CNVkit file
    log2_data = []
    cnvkit_df = pd.read_csv(normalised_coverage_path, sep="\t")
    for index, row in cnvkit_df.iterrows():
        if row["gene"] == "Antitarget":
            continue
        midpoint = row["start"] + int((row["end"] - row["start"]) / 2)
        log2_data.append(
            f"{row['chromosome']}\t{midpoint-1}\t{midpoint}\t{row['log2']}"
        )

    # Write log2 data to output file
    with open(output_file, "w") as log2_out:
        for resolution in ["o", "a", "b", "c", "d"]:
            for line in log2_data:
                log2_out.write(f"{resolution}_{line}\n")


if __name__ == "__main__":
    create_gens_cov_file()
