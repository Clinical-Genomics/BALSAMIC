import click
from BALSAMIC.utils.io import read_csv
import numpy as np
import warnings


def calculate_log2_ratio(purity, log2_ratio, ploidy):
    # Ensure that the inputs are within valid ranges
    if not (0 <= purity <= 1):
        raise ValueError("Purity must be between 0 and 1")

    if ploidy <= 0:
        raise ValueError("Ploidy must be a positive number")

    # Ploidy and purity adjustment formula reference to PureCN issue: https://github.com/lima1/PureCN/issues/40
    log2_adjusted = (
        purity * log2_ratio * ploidy + 2 * (1 - purity) * log2_ratio - 2 * (1 - purity)
    ) / (purity * ploidy)

    return log2_adjusted


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
    help="Input CNVkit cnr. result.",
)
@click.option(
    "-p",
    "--tumor-purity-path",
    required=True,
    type=click.Path(exists=True),
    help="Tumor purity file from PureCN",
)
def create_gens_cov_file(
    output_file: str, normalised_coverage_path: str, tumor_purity_path: str
):
    """Post-processes the CNVkit .cnr output for upload to GENS.

    Removes Antitarget regions, adjusts for purity and ploidy and outputs the coverages in multiple resolution-formats.

    Args:
        output_file: Path to GENS output.cov file
        normalised_coverage_path: Path to input CNVkit cnr file.
        tumor_purity_path: Path to PureCN purity estimate csv file
    """
    log2_data = []

    # Process CNVkit file
    cnr_dict_list: list[dict] = read_csv(
        csv_path=normalised_coverage_path, delimeter="\t"
    )

    # Process PureCN purity file
    purecn_dict_list: list[dict] = read_csv(csv_path=tumor_purity_path, delimeter=",")
    purity = float(purecn_dict_list[0]["Purity"])
    ploidy = float(purecn_dict_list[0]["Ploidy"])

    for row in cnr_dict_list:
        if row["gene"] == "Antitarget":
            continue

        # find midpoint
        start: int = int(row["start"])
        end: int = int(row["end"])
        region_size: int = end - start
        midpoint: int = start + int(region_size / 2)

        # adjust log2 ratio
        log2: float = float(row["log2"])
        log2: float = calculate_log2_ratio(purity, log2, ploidy)
        log2: float = round(log2, 4)

        # store values in list
        log2_data.append(f"{row['chromosome']}\t{midpoint - 1}\t{midpoint}\t{log2}")

    # Write log2 data to output file
    with open(output_file, "w") as log2_out:
        for resolution in ["o", "a", "b", "c", "d"]:
            for line in log2_data:
                log2_out.write(f"{resolution}_{line}\n")


if __name__ == "__main__":
    create_gens_cov_file()