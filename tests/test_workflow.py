from unittest import mock
import snakemake

from BALSAMIC.utils.cli import get_snakefile

MOCKED_OS_ENVIRON = "os.environ"


def test_workflow_tumor_only_tga_hg19(
    tumor_only_config, sentieon_install_dir, sentieon_license
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "single"
    reference_genome = "hg19"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_only_config

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


def test_workflow_tumor_normal_tga_hg19(
    tumor_normal_config, sentieon_install_dir, sentieon_license
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "paired"
    reference_genome = "hg19"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_normal_config

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


def test_workflow_tumor_only_wgs_hg19(
    tumor_only_wgs_config, sentieon_install_dir, sentieon_license
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "single"
    reference_genome = "hg19"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_only_wgs_config

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


def test_workflow_tumor_normal_wgs_hg19(
    tumor_normal_wgs_config, sentieon_install_dir, sentieon_license
):
    # GIVEN a sample config dict and a snakefile
    analysis_type = "paired"
    reference_genome = "hg19"
    analysis_workflow = "balsamic"
    snakefile = get_snakefile(analysis_type, analysis_workflow)
    config_json = tumor_normal_wgs_config

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


def test_workflow_qc_tumor_only_hg19(tumor_only_config):
    # GIVEN a sample config dict and a snakefile
    workflow = "paired"
    reference_genome = "hg19"
    analysis_workflow = "balsamic-qc"
    snakefile = get_snakefile(workflow, analysis_workflow)
    config_json = tumor_only_config

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for QC, hg19-tumor-only should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)


def test_workflow_qc_tumor_normal_hg19(tumor_normal_config):
    # GIVEN a sample config dict and a snakefile
    workflow = "paired"
    reference_genome = "hg19"
    analysis_workflow = "balsamic-qc"
    snakefile = get_snakefile(workflow, analysis_workflow)
    config_json = tumor_normal_config

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for QC, hg19-tumor-normal should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)


def test_workflow_qc_tumor_only_canfam3(tumor_only_config):
    # GIVEN a sample config dict and a snakefile
    workflow = "single"
    reference_genome = "canfam3"
    analysis_workflow = "balsamic-qc"
    snakefile = get_snakefile(workflow, analysis_workflow)
    config_json = tumor_only_config

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for QC, canfam3-tumor-only should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)


def test_workflow_qc_tumor_normal_canfam3(tumor_normal_config):
    # GIVEN a sample config dict and a snakefile
    workflow = "paired"
    reference_genome = "canfam3"
    analysis_workflow = "balsamic-qc"
    snakefile = get_snakefile(workflow, analysis_workflow)
    config_json = tumor_normal_config

    # WHEN invoking snakemake module with dry run option
    # THEN the snakemake workflow for QC, canfam3-tumor-normal should run successfully.
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
    ):
        assert snakemake.snakemake(snakefile, configfiles=[config_json], dryrun=True)
