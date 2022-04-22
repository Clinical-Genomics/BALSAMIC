import json
import os
from pathlib import Path
from unittest import mock
from BALSAMIC import __version__ as balsamic_version

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


def test_deliver_tumor_only_panel_files(
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

    # GIVEN a list of expected files to be included in the report
    delivery_files = [
        "sample_tumor_only.json",
        "sample_tumor_only_report.html",
        "sample_tumor_only_BALSAMIC_" + balsamic_version + "_graph.pdf",
        "multiqc_data.json",
        "sample_tumor_only_metrics_deliverables.yaml",
    ]

    # GIVEN a compound tag that should be added to the report
    compound_tags = [
        ["balsamic-config"],
        ["balsamic-report"],
        ["balsamic-dag"],
        ["multiqc-json", "json"],
        ["qc-metrics-yaml", "yaml"],
    ]

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

        # THEN it should run without any errors
        assert result.exit_code == 0
        assert actual_delivery_report.is_file()
        with open(actual_delivery_report, "r") as report:
            report = json.load(report)["files"]
            assert report

        # THEN all the expected files and tags should be included in the report
        for file in report:
            assert os.path.basename(file["path"]) in delivery_files
            assert file["tag"] in compound_tags


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
