from pathlib import Path
from unittest import mock


def test_status_tumor_only_panel(invoke_cli, tumor_only_config,
                                 sentieon_install_dir, sentieon_license):
    # GIVEN a tumor-only config file
    # WHEN running analysis
    with mock.patch.dict(
            'os.environ', {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        result = invoke_cli([
            'report', 'status', '--show-only-missing', '--sample-config',
            tumor_only_config
        ])

        # THEN it should run without any error
        assert result.exit_code == 0


def test_status_tumor_normal_panel(invoke_cli, tumor_normal_config, helpers,
                                   sentieon_install_dir, sentieon_license):
    # GIVEN a tumor-normal config file
    # WHEN running analysis with three actual delivery files
    # Actual delivery files dummies with and without index
    helpers.read_config(tumor_normal_config)
    normal_bam_result_dir = Path(helpers.result_dir, "bam")
    normal_bam_result_dir.mkdir(parents=True, exist_ok=True)
    normal_bam_delivery_file = Path(normal_bam_result_dir, "normal.merged.bam")
    normal_bam_delivery_file.touch()

    tumor_bam_result_dir = Path(helpers.result_dir, "bam")
    tumor_bam_result_dir.mkdir(parents=True, exist_ok=True)
    tumor_bam_delivery_file = Path(tumor_bam_result_dir, "tumor.merged.bam")
    tumor_bam_delivery_file.touch()

    with mock.patch.dict(
            'os.environ', {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        result = invoke_cli([
            'report', 'status', '--print-files', '--sample-config',
            tumor_normal_config
        ])

        # THEN it should run without any error
        assert result.exit_code == 0
