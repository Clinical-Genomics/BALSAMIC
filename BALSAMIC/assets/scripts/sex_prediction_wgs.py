import click
import json
from typing import List

def read_txt_file(filepath: str) -> List[List[str]]:
    """Read ascat statistics file and return rows as a list of lists."""
    with open(filepath, "r") as rf:
        rows = rf.readlines()
        return [r.strip("\n").split(" ") for r in rows]

def get_ascat_sex_prediction(ascat_sample_statistics_path, sample_type):

    sample_stats = read_txt_file(ascat_sample_statistics_path)
    for line in sample_stats:
        if line[0] == "GenderChrFound":
            male_y_or_n = line[1]

    if male_y_or_n == "Y":
        predicted_sex = "male"
    else:
        predicted_sex = "female"

    return {sample_type: {"predicted_sex": predicted_sex}}


def case_sex_prediction(predicted_sex):
    predicted_sex["case_sex"] = {}
    tumor_sex = predicted_sex["tumor"]["predicted_sex"]
    if "normal" in predicted_sex:
        normal_sex = predicted_sex["normal"]["predicted_sex"]
        if tumor_sex != normal_sex:
            case_sex = "conflicting"
        else:
            case_sex = tumor_sex
    else:
        case_sex = tumor_sex

    predicted_sex["case_sex"] = case_sex
    return predicted_sex


def write_json(json_obj: dict, path: str) -> None:
    """Write JSON format data to an output file."""
    try:
        with open(path, "w") as fn:
            json.dump(json_obj, fn, indent=4)
    except OSError as error:
        raise OSError(f"Error while writing JSON file: {path}, error: {error}")


@click.command()
@click.option(
    "--tumor-ascat-statistics",
    type=click.Path(exists=True),
    required=True,
    help="Path to the ascat statistics tumor file.",
)
@click.option(
    "--normal-ascat-statistics",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to the ascat statistics normal file.",
)
@click.option(
    "--output",
    type=click.Path(writable=True),
    required=True,
    help="Path to the output file to be created.",
)
def sex_check(
    tumor_ascat_statistics,
    normal_ascat_statistics,
    output,
):

    predicted_sex = get_ascat_sex_prediction(tumor_ascat_statistics, "tumor")

    if normal_ascat_statistics:
        normal_sex_prediction = get_ascat_sex_prediction(normal_ascat_statistics, "normal")
        predicted_sex.update(normal_sex_prediction)

    # Create case-level prediction (compare tumor and normal sex)
    predicted_sex = case_sex_prediction(predicted_sex)

    # Write json report
    write_json(predicted_sex, output)


if __name__ == "__main__":
    sex_check()
