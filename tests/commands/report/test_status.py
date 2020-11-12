from unittest import mock


def test_status_tumor_only_panel(invoke_cli, tumor_only_config,
                                 sentieon_install_dir, sentieon_license):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    with mock.patch.dict(
            'os.environ', {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        result = invoke_cli(
            ['report', 'status', '--sample-config', tumor_only_config])

        # THEN it should run without any error
        assert result.exit_code == 0


def test_status_tumor_normal_panel(invoke_cli, tumor_normal_config,
                                   sentieon_install_dir, sentieon_license):
    # GIVEN a tumor-normal config file
    # WHEN running analysis
    with mock.patch.dict(
            'os.environ', {
                'SENTIEON_LICENSE': sentieon_license,
                'SENTIEON_INSTALL_DIR': sentieon_install_dir
            }):
        result = invoke_cli(
            ['report', 'status', '--sample-config', tumor_normal_config])

        # THEN it should run without any error
        assert result.exit_code == 0
