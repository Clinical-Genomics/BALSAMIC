from pathlib import Path
from unittest import mock

import logging


def test_deliver_tumor_only_panel(
    invoke_cli,
    environ,
    tumor_only_config,
    helpers,
    sentieon_install_dir,
    sentieon_license,
    caplog,
):
    # GIVEN a tumor-normal config file
    helpers.read_config(tumor_only_config)
    actual_delivery_report = Path(helpers.delivery_dir, helpers.case_id + ".hk")
    multiqc_data_dir: Path = Path(helpers.result_dir, "qc", "multiqc_data")
    multiqc_data_dir.mkdir(parents=True, exist_ok=True)
    multiqc_data_json: Path = Path(multiqc_data_dir, "multiqc_data.json")
    multiqc_data_json.touch()

    with mock.patch.dict(
        environ,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ), caplog.at_level(logging.DEBUG):
        # WHEN running analysis
        result = invoke_cli(
            [
                "report",
                "deliver",
                "--sample-config",
                tumor_only_config,
                "--disable-variant-caller",
                "cnvkit",
            ]
        )

        # THEN it should run without any error
        assert result.exit_code == 0
        assert actual_delivery_report.is_file()
        assert "following" in caplog.text


def test_deliver_tumor_normal_panel(
    invoke_cli,
    environ,
    tumor_normal_config,
    helpers,
    sentieon_install_dir,
    sentieon_license,
    caplog,
):
    # GIVEN a tumor-normal config file
    helpers.read_config(tumor_normal_config)
    actual_delivery_report = Path(helpers.delivery_dir, helpers.case_id + ".hk")

    # Actual delivery files dummies with and without index
    cnv_result_dir = Path(helpers.result_dir, "cnv")
    cnv_result_dir.mkdir(parents=True, exist_ok=True)
    actual_delivery_file = Path(cnv_result_dir, "tumor.merged.cns")
    actual_delivery_file.touch()

    vep_result_dir = Path(helpers.result_dir, "vep")
    vep_result_dir.mkdir(parents=True, exist_ok=True)
    touch_vcf_delivery_file = Path(
        vep_result_dir, "SNV.somatic." + helpers.case_id + ".vardict.vcf.gz"
    )
    touch_vcf_delivery_file.touch()
    touch_vcf_delivery_file_index = Path(
        vep_result_dir, "SNV.somatic." + helpers.case_id + ".vardict.vcf.gz.tbi"
    )
    touch_vcf_delivery_file_index.touch()
    multiqc_data_dir: Path = Path(helpers.result_dir, "qc", "multiqc_data")
    multiqc_data_dir.mkdir(parents=True, exist_ok=True)
    multiqc_data_json: Path = Path(multiqc_data_dir, "multiqc_data.json")
    multiqc_data_json.touch()

    # Temporary files to be ignored by delivery
    vcf_result_dir = Path(helpers.result_dir, "vcf")
    vcf_result_dir.mkdir(parents=True, exist_ok=True)
    touch_temp_no_delivery_file = Path(
        vcf_result_dir, "CNV.somatic." + helpers.case_id + ".cnvkit.vcf.gz"
    )
    touch_temp_no_delivery_file.touch()

    with mock.patch.dict(
        environ,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ), caplog.at_level(logging.DEBUG):
        # WHEN running analysis
        result = invoke_cli(
            ["report", "deliver", "--sample-config", tumor_normal_config]
        )

        # THEN it should run without any error
        assert result.exit_code == 0
        assert actual_delivery_report.is_file()
        assert "following" in caplog.text
