from pathlib import Path

from BALSAMIC.assets.scripts.edit_vcf_info import edit_vcf_info

def test_edit_vcf_info(tmp_path, cli_runner):
    input_vcf_path = "tests/test_data/references/vep/test.vardict.all.filtered.pass.vcf.gz"
    output_vcf_path = tmp_path / "edited_vcf.gz"
    varcaller = "newcaller"

    result = cli_runner.invoke(
        edit_vcf_info, 
	 [
	 "--input_vcf",
	 input_vcf_path,
	 "--output_vcf", 
	 str(output_vcf_path),
	 "--variant_caller",
	 varcaller
	 ]
    )

    assert result.exit_code == 0
    assert Path(output_vcf_path).exists()

