def test_run_reference(invoke_cli, tmp_path, singularity_container):
    # Given test_reference.json
    test_new_dir = tmp_path / "test_reference_dir_with_run"
    test_new_dir.mkdir()

    test_output_reference_config = test_new_dir / "config.json"
    test_output_reference_pdf = test_new_dir / "generate_ref_worflow_graph.pdf"

    result_config = invoke_cli([
        'config', 'reference', '-c', 'secret_key', '--singularity',
        singularity_container, '-o',
        str(test_new_dir)
    ])

    # WHEN creating config.json in reference dir
    result_run = invoke_cli(
        ['run', 'reference', '-c',
         str(test_output_reference_config)])

    # THEN output config, pdf file generation, and reference dry run command exit code 0
    assert result_run.exit_code == 0
