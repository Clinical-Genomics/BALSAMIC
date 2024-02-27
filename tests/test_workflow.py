from unittest import mock
import logging
import snakemake

from BALSAMIC.constants.analysis import AnalysisWorkflow
from BALSAMIC.utils.cli import get_snakefile

MOCKED_OS_ENVIRON = "os.environ"


def test_workflow_tumor_only_tga_hg19(
    tumor_only_config,
    sentieon_install_dir,
    sentieon_license,
    caplog,
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "single"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_only_config
    caplog.set_level(logging.INFO)

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for TGA, hg19-tumor-only should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)

        print(caplog.text)

        # THEN the following rules should not be included
        # assert "igh_dux4_detection_tumor_only" not in caplog.text
        assert "cnvpytor_tumor_only" not in caplog.text


def test_workflow_tumor_normal_tga_hg19(
    tumor_normal_config,
    sentieon_install_dir,
    sentieon_license,
    caplog,
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "paired"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_normal_config
    caplog.set_level(logging.INFO)

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for TGA, hg19-tumor-normal should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)

        # THEN the following rules should not be included
        # assert "igh_dux4_detection_tumor_normal" not in caplog.text


def test_workflow_tumor_only_wgs_hg19(
    tumor_only_wgs_config,
    sentieon_install_dir,
    sentieon_license,
    caplog,
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "single"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_only_wgs_config
    caplog.set_level(logging.INFO)

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for WGS, hg19-tumor-only should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)

        # THEN the following rules should be included
        # assert "igh_dux4_detection_tumor_only" in caplog.text


def test_workflow_tumor_normal_wgs_hg19(
    tumor_normal_wgs_config,
    sentieon_install_dir,
    sentieon_license,
    caplog,
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "paired"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_normal_wgs_config
    caplog.set_level(logging.INFO)

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for WGS, hg19-tumor-normal should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)

        # THEN the following rules should be included
        # assert "igh_dux4_detection_tumor_normal" in caplog.text


def test_workflow_qc_tumor_only_hg19(
    tumor_only_config_qc, sentieon_install_dir, sentieon_license
):
    # GIVEN a sample config dict and a snakefile
    workflow = "single"
    snakefile = get_snakefile(workflow, AnalysisWorkflow.BALSAMIC_QC)
    config_json = tumor_only_config_qc

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for QC, hg19-tumor-only should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)


def test_workflow_qc_tumor_normal_hg19(
    tumor_normal_config_qc, sentieon_install_dir, sentieon_license
):
    # GIVEN a sample config dict and a snakefile
    workflow = "paired"
    snakefile = get_snakefile(workflow, AnalysisWorkflow.BALSAMIC_QC)
    config_json = tumor_normal_config_qc

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for QC, hg19-tumor-normal should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)
