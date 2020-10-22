def test_init_container(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    # WHEN creating config.json in reference dir
    result = invoke_cli(
        ['init', 'container', '--dry', '--out-dir',
         str(test_new_dir)])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0


def test_init_container_force(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()

    # WHEN creating config.json in reference dir
    result = invoke_cli([
        'init', 'container', '--force', '--dry', '--out-dir',
        str(test_new_dir)
    ])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0


def test_init_container_specific_tag(invoke_cli, tmp_path):
    # Given a dummy path
    test_new_dir = tmp_path / "test_container_dir"
    test_new_dir.mkdir()
    dummy_tag = "cool_new_feature"

    # WHEN creating config.json in reference dir
    result = invoke_cli([
        'init', 'container', '--force', '--dry', '--container-version',
        dummy_tag, '--out-dir',
        str(test_new_dir)
    ])

    # THEN output config and pdf file generate and command exit code 0
    assert result.exit_code == 0
