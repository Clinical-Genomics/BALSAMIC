from pathlib import Path

from BALSAMIC.assets.scripts.edit_vcf_info import edit_vcf_info


def test_edit_vcf_info(tmp_path, cli_runner):
    """test cli command is properly executed  and returns output file"""

    # GIVEN input VCF file path, path to output file and varaintcaller name
    input_vcf_path = (
        "tests/test_data/references/vep/test.vardict.all.filtered.pass.vcf.gz"
    )
    output_vcf_path = tmp_path / "edited_vcf.gz"
    varcaller = "newcaller"

    # WHEN calling edit_vcf_info script
    result = cli_runner.invoke(
        edit_vcf_info,
        [
            "--input_vcf",
            input_vcf_path,
            "--output_vcf",
            str(output_vcf_path),
            "--variant_caller",
            varcaller,
        ],
    )

    # THEN it should return exit code and checks for output file paths
    assert result.exit_code == 0
    assert Path(output_vcf_path).exists()
