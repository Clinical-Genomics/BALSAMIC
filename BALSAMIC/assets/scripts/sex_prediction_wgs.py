import click
import json
from typing import List
import numpy as np

def read_txt_file(filepath: str) -> List[List[str]]:
    """Read ascat statistics file and return rows as a list of lists."""
    with open(filepath, "r") as rf:
        rows = rf.readlines()
        return [r.strip("\n").split(" ") for r in rows]

def get_ascat_sex_prediction(ascat_sample_statistics_path):

    sample_stats = read_txt_file(ascat_sample_statistics_path)
    for line in sample_stats:
        if line[0] == "GenderChrFound":
            male_y_or_n = line[1]

    if male_y_or_n == "Y":
        predicted_sex = "male"
    else:
        predicted_sex = "female"

    return {"case_sex": predicted_sex}

def predict_sex_from_y_and_x_cov(tumor_y_coverage_path: str, tumor_x_coverage_path: str) -> dict:
    """Compute statistics on a column of numbers in a file.

    Args:
        file_path (str): Path to the input file containing a column of numbers.

    Returns:
        dict: A dictionary containing mean, median, standard deviation, min, and max.
    """
    # Load the file
    tumor_y_cov = np.loadtxt(tumor_y_coverage_path)

    # Compute median
    tumor_median_y_cov = np.median(tumor_y_cov)

    # Load the file
    tumor_x_cov = np.loadtxt(tumor_x_coverage_path)

    # Compute median
    tumor_median_x_cov = np.median(tumor_x_cov)

    # Compute fraction
    tumor_y_x_median_frac = round(tumor_median_y_cov / tumor_median_x_cov, 5)

    # Decide sex
    if tumor_y_x_median_frac > 0.1:
        predicted_sex = "male"
    elif tumor_y_x_median_frac < 0.08:
        predicted_sex = "female"
    else:
        predicted_sex = "unknown"

    return {"case_sex": predicted_sex,
            "tumor_median_x_coverage": tumor_median_x_cov,
            "tumor_median_y_coverage": tumor_median_y_cov,
            "tumor_y_x_median_frac": tumor_y_x_median_frac
            }

def write_json(json_obj: dict, path: str) -> None:
    """Write JSON format data to an output file."""
    try:
        with open(path, "w") as fn:
            json.dump(json_obj, fn, indent=4)
    except OSError as error:
        raise OSError(f"Error while writing JSON file: {path}, error: {error}")


@click.command()
@click.option(
    "--case-ascat-statistics",
    type=click.Path(exists=True),
    required=False,
    help="Optional path to the ascat statistics file (for TN cases)",
)
@click.option(
    "--sample-y-coverage",
    type=click.Path(exists=True),
    required=False,
    help="Optional path to file with coverage per base in Y chromosome (for T only cases)",
)
@click.option(
    "--sample-x-coverage",
    type=click.Path(exists=True),
    required=False,
    help="Optional path to file with coverage per base in X chromosome (for T only cases)",
)
@click.option(
    "--output",
    type=click.Path(writable=True),
    required=True,
    help="Path to the output file to be created.",
)
def sex_check(
    case_ascat_statistics,
    sample_y_coverage,
    sample_x_coverage,
    output,
):
    if case_ascat_statistics:
        predicted_sex = get_ascat_sex_prediction(case_ascat_statistics)
    elif sample_y_coverage and sample_x_coverage:
        predicted_sex = predict_sex_from_y_and_x_cov(sample_y_coverage, sample_x_coverage)
    else:
        print("Missing input files to predict sample sex, exiting")
        return

    # Write json report
    write_json(predicted_sex, output)


if __name__ == "__main__":
    sex_check()
