import click
import json
from typing import List

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
    required=True,
    help="Path to the ascat statistics file.",
)
@click.option(
    "--output",
    type=click.Path(writable=True),
    required=True,
    help="Path to the output file to be created.",
)
def sex_check(
    case_ascat_statistics,
    output,
):
    predicted_sex = get_ascat_sex_prediction(case_ascat_statistics)

    # Write json report
    write_json(predicted_sex, output)


if __name__ == "__main__":
    sex_check()
