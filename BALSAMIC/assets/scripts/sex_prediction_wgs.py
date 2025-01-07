import click
import json
from typing import List, Dict
import numpy as np

def read_txt_file(filepath: str) -> List[List[str]]:
    """
    Reads a text file containing space-separated values and returns the data as a list of lists.

    Args:
        filepath (str): Path to the input text file.

    Returns:
        List[List[str]]: A list of lists where each inner list represents a row of the file.
    """
    with open(filepath, "r") as rf:
        rows = rf.readlines()
    return [r.strip("\n").split(" ") for r in rows]

def get_ascat_sex_prediction(ascat_sample_statistics_path: str) -> Dict[str, str]:
    """
    Predicts the sex of a sample based on ASCAT statistics.

    Args:
        ascat_sample_statistics_path (str): Path to the ASCAT statistics file.

    Returns:
        Dict[str, str]: A dictionary containing the predicted sex of the case.
    """
    sample_stats = read_txt_file(ascat_sample_statistics_path)
    male_y_or_n = "N"  # Default value if "GenderChrFound" line is missing
    for line in sample_stats:
        if line[0] == "GenderChrFound":
            male_y_or_n = line[1]
            break

    predicted_sex = "male" if male_y_or_n == "Y" else "female"
    return {"case_sex": predicted_sex}

def predict_sex_from_y_and_x_cov(
    tumor_y_coverage_path: str, tumor_x_coverage_path: str
) -> Dict[str, float]:
    """
    Predicts the sex of a sample based on the median coverage of Y and X chromosomes.

    Args:
        tumor_y_coverage_path (str): Path to the file containing Y chromosome coverage values.
        tumor_x_coverage_path (str): Path to the file containing X chromosome coverage values.

    Returns:
        Dict[str, float]: A dictionary containing predicted sex, median coverages, and their ratio.
    """
    tumor_y_cov = np.loadtxt(tumor_y_coverage_path)
    tumor_x_cov = np.loadtxt(tumor_x_coverage_path)

    tumor_median_y_cov = np.median(tumor_y_cov)
    tumor_median_x_cov = np.median(tumor_x_cov)
    tumor_y_x_median_frac = round(tumor_median_y_cov / tumor_median_x_cov, 5)

    if tumor_y_x_median_frac > 0.1:
        predicted_sex = "male"
    elif tumor_y_x_median_frac < 0.08:
        predicted_sex = "female"
    else:
        predicted_sex = "unknown"

    return {
        "case_sex": predicted_sex,
        "tumor_median_x_coverage": tumor_median_x_cov,
        "tumor_median_y_coverage": tumor_median_y_cov,
        "tumor_y_x_median_frac": tumor_y_x_median_frac,
    }

def write_json(json_obj: Dict, path: str) -> None:
    """
    Writes a dictionary to a JSON file.

    Args:
        json_obj (Dict): The data to be written as JSON.
        path (str): Path to the output file.

    Raises:
        OSError: If there is an error writing the file.
    """
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
    help="Optional path to the ASCAT statistics file (for TN cases)",
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
    case_ascat_statistics: str,
    sample_y_coverage: str,
    sample_x_coverage: str,
    output: str,
) -> None:
    """
    Determines the sex of a sample based on ASCAT statistics or chromosome coverage.

    Args:
        case_ascat_statistics (str): Path to the ASCAT statistics file.
        sample_y_coverage (str): Path to the Y chromosome coverage file.
        sample_x_coverage (str): Path to the X chromosome coverage file.
        output (str): Path to the output file for the results.
    """
    if case_ascat_statistics:
        predicted_sex = get_ascat_sex_prediction(case_ascat_statistics)
    elif sample_y_coverage and sample_x_coverage:
        predicted_sex = predict_sex_from_y_and_x_cov(
            sample_y_coverage, sample_x_coverage
        )
    else:
        click.echo("Missing input files to predict sample sex, exiting", err=True)
        return

    write_json(predicted_sex, output)

if __name__ == "__main__":
    sex_check()
