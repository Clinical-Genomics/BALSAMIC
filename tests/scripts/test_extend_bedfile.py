"""Test extending bedfile."""
from pathlib import Path

from click.testing import CliRunner, Result

from BALSAMIC.assets.scripts.extend_bedfile import extend_bedfile
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS


def test_extend_bedfile(
    bedfile_path: Path, tmp_path: Path, cli_runner: CliRunner
) -> None:
    """Test extended bedfile with regions shorter than minimum_region_size used in CNV analysis."""

    # GIVEN an input bedfile with regions smaller than minimum region size

    # GIVEN a minimum region size
    minimum_region_size = WORKFLOW_PARAMS["expand_short_bedregions"][
        "minimum_region_size"
    ]

    # GIVEN an output bedfile
    extended_bed_path: Path = Path(tmp_path, "test_extended_bed.bed")

    # WHEN running the extend bedfile script
    result: Result = cli_runner.invoke(
        extend_bedfile,
        [
            "--min-region-size",
            minimum_region_size,
            bedfile_path.as_posix(),
            extended_bed_path.as_posix(),
        ],
    )

    # THEN the output bedfile should exist
    assert result.exit_code == 0
    assert extended_bed_path.is_file()

    with open(bedfile_path, "r") as input_bedfile:
        lines_before = []
        for line in input_bedfile:
            chrom_start_end = line.strip().split("\t")
            lines_before.append(chrom_start_end)

    with open(extended_bed_path, "r") as output_bedfile:
        lines_after = []
        for line in output_bedfile:
            chrom_start_end = line.strip().split("\t")
            lines_after.append(chrom_start_end)

    # THEN the first row should be unchanged
    assert lines_before[0] == lines_after[0]

    # THEN the second row should be extended by 99 bases
    assert lines_before[1] == [1, 10005000, 10005001]
    assert lines_after[1] == [1, 10004950, 10005050]
    region_size = lines_after[1][2] - lines_after[1][1]
    assert region_size == minimum_region_size
