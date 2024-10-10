"""Test extending bedfile."""
from pathlib import Path

from click.testing import CliRunner, Result

from BALSAMIC.assets.scripts.postprocess_gens_cnvkit import create_gens_cov_file


def test_create_gens_cov_file(
    gens_dummy_cnvkit_cnr,
    gens_dummy_purecn_purity,
    gens_dummy_cov_bed_expected,
    tmp_path: Path,
    cli_runner: CliRunner,
) -> None:
    """Test postprocess CNVkit TGA output for GENS script."""

    # GIVEN an input cnvkit cnr file and a pureCN purity csv

    # GIVEN a path to an outputfile
    gens_cnvkit_cov: Path = Path(tmp_path, "gens_tga.cov.bed")

    # WHEN running the postprocess gens cnvkit script
    result: Result = cli_runner.invoke(
        create_gens_cov_file,
        [
            "--output-file",
            gens_cnvkit_cov.as_posix(),
            "--normalised-coverage-path",
            gens_dummy_cnvkit_cnr,
            "--tumor-purity-path",
            gens_dummy_purecn_purity,
        ],
    )

    # THEN the output bedfile should exist
    assert result.exit_code == 0
    assert gens_cnvkit_cov.is_file()

    # WHEN reading produced output file and expected output file
    with open(gens_cnvkit_cov, "r") as actual_file:
        test_output = actual_file.read()

    with open(gens_dummy_cov_bed_expected, "r") as expected_file:
        expected_output = expected_file.read()

    # THEN test file and expected test file should be identical
    assert (
        test_output == expected_output
    ), "The expected and produced TGA GENS files do not match."
