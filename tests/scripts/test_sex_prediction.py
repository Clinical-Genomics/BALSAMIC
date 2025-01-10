from pathlib import Path
from BALSAMIC.assets.scripts.sex_prediction_tga import predict_sex_main
from BALSAMIC.assets.scripts.sex_prediction_wgs import predict_sex_wgs
from BALSAMIC.utils.io import read_json


def test_tga_male_sex_prediction(
    male_target_cnn_file: str, male_antitarget_cnn_file: str, tmp_path, cli_runner
):
    """Ensure that TGA sex prediction is working using cnvkit cnn files for male sample."""
    # GIVEN the output path, and input cnn files
    output_path = str(tmp_path / "male_sex_prediction.json")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        predict_sex_main,
        [
            "--target-cnn-tumor",
            male_target_cnn_file,
            "--antitarget-cnn-tumor",
            male_antitarget_cnn_file,
            "--output",
            output_path,
        ],
    )

    # THEN check if the JSON is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()

    # THEN check that male sex has been correctly predicted
    male_sex_prediction: dict = read_json(output_path)
    assert male_sex_prediction["case_sex"] == "male"


def test_tga_female_sex_prediction(
    female_target_cnn_file: str, female_antitarget_cnn_file: str, tmp_path, cli_runner
):
    """Ensure that TGA sex prediction is working using cnvkit cnn files for female sample."""
    # GIVEN the output path, and input cnn files
    output_path = str(tmp_path / "female_sex_prediction.json")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        predict_sex_main,
        [
            "--target-cnn-tumor",
            female_target_cnn_file,
            "--antitarget-cnn-tumor",
            female_antitarget_cnn_file,
            "--output",
            output_path,
        ],
    )

    # THEN check if the JSON is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()

    # THEN check that female sex has been correctly predicted
    female_sex_prediction: dict = read_json(output_path)
    assert female_sex_prediction["case_sex"] == "female"


def test_tga_conflicting_sex_prediction(
    male_target_cnn_file: str,
    male_antitarget_cnn_file: str,
    female_target_cnn_file: str,
    female_antitarget_cnn_file: str,
    tmp_path,
    cli_runner,
):
    """Ensure that TGA sex is reported as conflicting when using mixed up samples of opposite sexes."""
    # GIVEN the output path, and input cnn files
    output_path = str(tmp_path / "conflicting_sex_prediction.json")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        predict_sex_main,
        [
            "--target-cnn-tumor",
            male_target_cnn_file,
            "--antitarget-cnn-tumor",
            male_antitarget_cnn_file,
            "--target-cnn-normal",
            female_target_cnn_file,
            "--antitarget-cnn-normal",
            female_antitarget_cnn_file,
            "--output",
            output_path,
        ],
    )

    # THEN check if the JSON is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()

    # THEN check that conflicting sex has been correctly predicted
    conflicting_sex_prediction: dict = read_json(output_path)
    assert conflicting_sex_prediction["case_sex"] == "conflicting"


def test_wgs_tumor_only_male_sex_prediction(
    male_200k_x_coverage: str, male_200k_y_coverage: str, tmp_path, cli_runner
):
    """Ensure that WGS sex prediction is working for tumor only using per base X and Y coverages for male sample."""
    # GIVEN the output path, and input cnn files
    output_path = str(tmp_path / "male_wgs_sex_prediction.json")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        predict_sex_wgs,
        [
            "--sample-x-coverage",
            male_200k_x_coverage,
            "--sample-y-coverage",
            male_200k_y_coverage,
            "--output",
            output_path,
        ],
    )

    # THEN check if the JSON is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()

    # THEN check that male sex has been correctly predicted
    sex_prediction: dict = read_json(output_path)
    assert sex_prediction["case_sex"] == "male"


def test_wgs_tumor_only_female_sex_prediction(
    female_200k_x_coverage: str, female_200k_y_coverage: str, tmp_path, cli_runner
):
    """Ensure that WGS sex prediction is working for tumor only using per base X and Y coverages for female sample."""
    # GIVEN the output path, and input cnn files
    output_path = str(tmp_path / "female_wgs_sex_prediction.json")

    # WHEN invoking the python script
    result = cli_runner.invoke(
        predict_sex_wgs,
        [
            "--sample-x-coverage",
            female_200k_x_coverage,
            "--sample-y-coverage",
            female_200k_y_coverage,
            "--output",
            output_path,
        ],
    )

    # THEN check if the JSON is correctly created and there are no errors
    assert result.exit_code == 0
    assert Path(output_path).exists()

    # THEN check that female sex has been correctly predicted
    sex_prediction: dict = read_json(output_path)
    assert sex_prediction["case_sex"] == "female"
