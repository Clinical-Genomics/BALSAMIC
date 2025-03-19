import click
import json
from typing import List, Dict
import numpy as np


def predict_sex_from_y_and_x_cov(
    y_coverage_path: str, x_coverage_path: str, sample_type: str
) -> Dict[str, float]:
    """
    Predicts the sex of a sample based on the median coverage of Y and X chromosomes.

    Args:
        y_coverage_path (str): Path to the file containing Y chromosome coverage values.
        x_coverage_path (str): Path to the file containing X chromosome coverage values.

    Returns:
        Dict[str, float]: A dictionary containing predicted sex, median coverages, and their ratio.
    """
    y_cov = np.loadtxt(y_coverage_path)
    x_cov = np.loadtxt(x_coverage_path)

    median_y_cov = np.median(y_cov)
    median_x_cov = np.median(x_cov)
    y_x_median_frac = round(median_y_cov / median_x_cov, 5)

    if y_x_median_frac > 0.1:
        predicted_sex = "male"
    elif y_x_median_frac < 0.08:
        predicted_sex = "female"
    else:
        predicted_sex = "unknown"

    return {
        sample_type: {
            "predicted_sex": predicted_sex,
            "median_x_coverage": median_x_cov,
            "median_y_coverage": median_y_cov,
            "y_x_median_frac": y_x_median_frac,
        }
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
    "--tumor-y-coverage",
    type=click.Path(exists=True),
    required=True,
    help="Optional path to file with coverage per base in Y chromosome (for T only cases)",
)
@click.option(
    "--tumor-x-coverage",
    type=click.Path(exists=True),
    required=True,
    help="Optional path to file with coverage per base in X chromosome (for T only cases)",
)
@click.option(
    "--normal-y-coverage",
    type=click.Path(exists=True),
    required=False,
    help="Optional path to file with coverage per base in Y chromosome (for T only cases)",
)
@click.option(
    "--normal-x-coverage",
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
def predict_sex_wgs(
    tumor_y_coverage: str,
    tumor_x_coverage: str,
    normal_y_coverage: str,
    normal_x_coverage: str,
    output: str,
) -> None:
    """
    Determines the sex of a sample based on ASCAT statistics or chromosome coverage.

    Args:
        tumor_y_coverage (str): Path to the Y chromosome coverage file for the tumor sample.
        tumor_x_coverage (str): Path to the X chromosome coverage file for the tumor sample.
        normal_y_coverage (str): Path to the Y chromosome coverage file for the normal sample.
        normal_x_coverage (str): Path to the X chromosome coverage file for the normal sample.
        output (str): Path to the output file for the results.
    """
    predicted_sex = predict_sex_from_y_and_x_cov(
        tumor_y_coverage, tumor_x_coverage, "tumor"
    )

    if normal_y_coverage and normal_x_coverage:
        predicted_normal_sex = predict_sex_from_y_and_x_cov(
            normal_y_coverage, normal_x_coverage, "normal"
        )
        predicted_sex.update(predicted_normal_sex)

    write_json(predicted_sex, output)


if __name__ == "__main__":
    predict_sex_wgs()
