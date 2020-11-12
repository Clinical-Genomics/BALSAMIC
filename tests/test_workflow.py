from unittest import mock
import snakemake

from BALSAMIC.utils.cli import get_snakefile

MOCKED_OS_ENVIRON = 'os.environ'


def test_workflow_tumor_normal(tumor_normal_config, sentieon_install_dir,
                               sentieon_license):
    # GIVEN a sample config dict and snakefile
    workflow = 'paired'
    snakefile = get_snakefile(workflow)
    config_json = tumor_normal_config

    # WHEN invoking snakemake module with dryrun option
    # THEN it should return true
    with mock.patch.dict(
            MOCKED_OS_ENVIRON, {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        assert snakemake.snakemake(snakefile,
                                   configfiles=[config_json],
                                   dryrun=True)


def test_workflow_tumor_only(tumor_only_config, sentieon_install_dir,
                             sentieon_license):
    # GIVEN a sample config dict and snakefile
    workflow = 'single'
    snakefile = get_snakefile(workflow)
    config_json = tumor_only_config

    # WHEN invoking snakemake module with dryrun option
    # THEN it should return true
    with mock.patch.dict(
            MOCKED_OS_ENVIRON, {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        assert snakemake.snakemake(snakefile,
                                   configfiles=[config_json],
                                   dryrun=True)


def test_workflow_qc(tumor_normal_config, tumor_only_config,
                     sentieon_install_dir, sentieon_license):
    # GIVEN a sample config dict and snakefile
    workflow = 'qc'
    snakefile = get_snakefile(workflow)

    # WHEN invoking snakemake module with dryrun option
    # THEN it should return true
    with mock.patch.dict(
            MOCKED_OS_ENVIRON, {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        for config_json in (tumor_normal_config, tumor_only_config):
            assert snakemake.snakemake(snakefile,
                                       configfiles=[config_json],
                                       dryrun=True)


def test_workflow_sentieon(tumor_normal_wgs_config, tumor_only_wgs_config,
                           sentieon_install_dir, sentieon_license):
    # GIVEN a sample config dict and snakefile
    workflows = [('single', tumor_only_wgs_config),
                 ('paired', tumor_normal_wgs_config)]
    sequencing_type = "wgs"

    # WHEN invoking snakemake module with dryrun option
    # THEN it should return true
    with mock.patch.dict(
            MOCKED_OS_ENVIRON, {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        for workflow in workflows:
            analysis_type = workflow[0]
            config = workflow[1]
            snakefile = get_snakefile(analysis_type, sequencing_type)
            assert snakemake.snakemake(snakefile,
                                       configfiles=[config],
                                       dryrun=True)


def test_umiworkflow_tumor_only(tumor_only_umi_config):
    # GIVEN a sample config dict and snakefile
    workflow = 'umi'
    snakefile = get_snakefile(workflow)
    config_json = tumor_only_umi_config
    config_tumorlod = "tests/test_data/references/tumorlod.json"

    # WHEN invoking snakemake module with dryrun option
    # THEN it should return true
    with mock.patch.dict(
        MOCKED_OS_ENVIRON, {
            'SENTIEON_LICENSE': sentieon_license,
            'SENTIEON_INSTALL_DIR': sentieon_install_dir
        }):
        assert snakemake.snakemake(snakefile,
                                   configfiles=[config_json, config_tumorlod],
                                   dryrun=True)
