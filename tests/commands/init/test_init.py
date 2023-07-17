"""Test Balsamic init command."""
from functools import partial
from pathlib import Path
from unittest import mock

from click.testing import Result
from graphviz import Source

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.analysis import RunMode
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.cluster import ClusterAccount
from BALSAMIC.constants.constants import EXIT_SUCCESS, EXIT_FAIL


def test_init_hg(
    invoke_cli: partial,
    tmp_path: Path,
    cosmic_key: str,
    config_json: str,
    reference_graph: str,
):
    """Test Balsamic init command."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.HG19,
            "--cosmic-key",
            cosmic_key,
        ]
    )

    # THEN the human reference generation workflow should have successfully started
    assert Path(tmp_path, balsamic_version, GenomeVersion.HG19, config_json).exists()
    assert Path(
        tmp_path, balsamic_version, GenomeVersion.HG19, reference_graph
    ).exists()
    assert result.exit_code == EXIT_SUCCESS


def test_init_canfam(
    invoke_cli: partial,
    tmp_path: Path,
    cosmic_key: str,
    config_json: str,
    reference_graph: str,
):
    """Test Balsamic canine workflow init command."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.CanFam3,
        ]
    )

    # THEN the canine reference generation workflow should have successfully started
    assert Path(tmp_path, balsamic_version, GenomeVersion.CanFam3, config_json).exists()
    assert Path(
        tmp_path, balsamic_version, GenomeVersion.CanFam3, reference_graph
    ).exists()
    assert result.exit_code == EXIT_SUCCESS


def test_init_hg_no_cosmic_key(invoke_cli: partial, tmp_path: Path, cosmic_key: str):
    """Test Balsamic init command when a COSMIC key is not provided."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.HG19,
        ]
    )

    # THEN an exception should have been raised
    assert (
        f"No COSMIC authentication key specified. It is required when using {GenomeVersion.HG19} reference"
        in result.output
    )
    assert result.exit_code == EXIT_FAIL


def test_init_hg_run_analysis(
    invoke_cli: partial,
    tmp_path: Path,
    cosmic_key: str,
    config_json: str,
    reference_graph: str,
):
    """Test Balsamic init command when actually running the analysis."""

    # GIVEN a temporary output directory, a cluster account, and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.HG19,
            "--cosmic-key",
            cosmic_key,
            "--run-mode",
            RunMode.CLUSTER,
            "--account",
            ClusterAccount.DEVELOPMENT,
            "--run-analysis",
        ]
    )

    # THEN the human reference generation workflow should have successfully started
    assert Path(tmp_path, balsamic_version, GenomeVersion.HG19, config_json).exists()
    assert Path(
        tmp_path, balsamic_version, GenomeVersion.HG19, reference_graph
    ).exists()
    assert result.exit_code == EXIT_SUCCESS


def test_init_hg_run_analysis_no_account(
    invoke_cli: partial, tmp_path: Path, cosmic_key: str
):
    """Test Balsamic init command when actually running the analysis without specifying a cluster account."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    result: Result = invoke_cli(
        [
            "init",
            "--out-dir",
            tmp_path.as_posix(),
            "--genome-version",
            GenomeVersion.HG19,
            "--cosmic-key",
            cosmic_key,
            "--run-mode",
            RunMode.CLUSTER,
            "--run-analysis",
        ]
    )

    # THEN an exception should have been raised
    assert "A cluster account is required for cluster run mode" in result.output
    assert result.exit_code == EXIT_FAIL


def test_init_hg_graph_exception(
    invoke_cli: partial,
    tmp_path: Path,
    cosmic_key: str,
    config_json: str,
    reference_graph: str,
):
    """Test Balsamic init command with a graphviz exception."""

    # GIVEN a temporary output directory and a COSMIC key

    # WHEN invoking the init command
    with mock.patch.object(Source, "render"):
        result: Result = invoke_cli(
            [
                "init",
                "--out-dir",
                tmp_path.as_posix(),
                "--genome-version",
                GenomeVersion.HG19,
                "--cosmic-key",
                cosmic_key,
            ]
        )

    # THEN the human reference generation workflow should have successfully started
    assert "Workflow graph generation failed" in result.output
    assert Path(tmp_path, balsamic_version, GenomeVersion.HG19, config_json).exists()
    assert not Path(
        tmp_path, balsamic_version, GenomeVersion.HG19, reference_graph
    ).exists()
    assert result.exit_code == EXIT_FAIL
