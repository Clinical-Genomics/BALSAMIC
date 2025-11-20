import importlib
from pathlib import Path
from click.testing import CliRunner
import subprocess
import logging
import pytest


MODULE_NAME = "BALSAMIC.assets.scripts.scan_finished_jobid_status"


@pytest.fixture(scope="module")
def m():
    # Import the target module once for all tests
    return importlib.import_module(MODULE_NAME)


# ---------- find_job_logs ----------


def test_find_job_logs_filters_numeric(tmp_path: Path, m, caplog):
    caplog.set_level(logging.DEBUG)
    # numeric names
    (tmp_path / "111.log").write_text("X")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "222.log").write_text("Y")
    # non-numeric or wrong suffix
    (tmp_path / "abc.log").write_text("Z")
    (tmp_path / "333.txt").write_text("nope")

    logs = m.find_job_logs(tmp_path)
    assert logs == {
        "111": tmp_path / "111.log",
        "222": tmp_path / "sub" / "222.log",
    }
    # One debug for the non-job log
    assert any("Skipping non-job log file" in rec.message for rec in caplog.records)
    assert any("Discovered 2 job logs" in rec.message for rec in caplog.records)


# ---------- get_job_state (per-job sacct wrapper) ----------


def test_get_job_state_success(monkeypatch, m, caplog):
    caplog.set_level(logging.DEBUG)

    class Dummy:
        # Simulate:
        # sacct --parsable2 -o JobID,State
        # 42|COMPLETED
        stdout = "42|COMPLETED\n"

    def fake_run(*args, **kwargs):
        return Dummy()

    monkeypatch.setattr(subprocess, "run", fake_run)

    out = m.get_job_state("42")
    # We expect the normalized state string
    assert out == "COMPLETED"
    # And a debug log about normalization
    assert any(
        "Job 42 sacct raw state='COMPLETED' -> normalized='COMPLETED'" in rec.message
        for rec in caplog.records
    )


def test_get_job_state_prefers_parent_and_normalizes(monkeypatch, m, caplog):
    caplog.set_level(logging.DEBUG)

    class Dummy:
        # Here we have multiple records:
        # - parent: 42|FAILED node_fail
        # - step:   42.0|CANCELLED+
        stdout = "42|FAILED node_fail\n42.0|CANCELLED+\n"

    def fake_run(*args, **kwargs):
        return Dummy()

    monkeypatch.setattr(subprocess, "run", fake_run)

    out = m.get_job_state("42")
    # Parent record is chosen, and state is normalized to the first token without '+'
    assert out == "FAILED"
    assert any(
        "Job 42 sacct raw state='FAILED node_fail' -> normalized='FAILED'"
        in rec.message
        for rec in caplog.records
    )


def test_get_job_state_not_found_executable(monkeypatch, m, caplog):
    caplog.set_level(logging.ERROR)

    def boom(*args, **kwargs):
        raise FileNotFoundError("no sacct")

    monkeypatch.setattr(subprocess, "run", boom)

    assert m.get_job_state("7") is None
    assert any("sacct executable not found" in rec.message for rec in caplog.records)


def test_get_job_state_calledprocesserror(monkeypatch, m, caplog):
    caplog.set_level(logging.DEBUG)

    def boom(*args, **kwargs):
        raise subprocess.CalledProcessError(
            returncode=1, cmd=["sacct"], stderr="bad things"
        )

    monkeypatch.setattr(subprocess, "run", boom)

    assert m.get_job_state("99") is None
    assert any(
        "Could not check job 99 via sacct" in rec.message for rec in caplog.records
    )
    assert any(
        "sacct stderr for 99: bad things" in rec.message for rec in caplog.records
    )


# ---------- write_results ----------


def test_write_results_all_sections(tmp_path: Path, m):
    out = tmp_path / "result.txt"

    failed = [
        ("101", Path("/logs/101.log")),
        ("202", Path("/logs/202.log")),
    ]
    cancelled = [
        ("303", Path("/logs/303.log")),
    ]
    unknown = ["404", "505"]

    # New section in write_results: jobs that failed but later succeeded
    resolved_failures = [
        ("909", Path("/logs/909.log"), "FAILED", "rule_x"),
    ]

    m.write_results(out, failed, cancelled, unknown, resolved_failures)
    txt = out.read_text()

    # Timestamp header present (we don't assert the exact time)
    assert "=== Job status check at " in txt

    # Failed block (note: new heading text)
    assert "FAILED JOBS (no successful retry):" in txt
    assert "101\t/logs/101.log" in txt
    assert "202\t/logs/202.log" in txt

    # Cancelled block (note: new heading text)
    assert "CANCELLED JOBS (no successful retry):" in txt
    assert "303\t/logs/303.log" in txt

    # Resolved failures block
    assert "Some jobs failed but succeeded on retry" in txt
    # Example job
    assert "909\t/logs/909.log\tFAILED" in txt

    # Unknown block
    assert "Unknown status jobs:" in txt
    assert "404\tNA" in txt
    assert "505\tNA" in txt

    # No SUCCESSFUL marker because there are unresolved failures/cancellations
    assert "SUCCESSFUL (no unresolved failed/cancelled jobs)" not in txt


def test_write_results_success_only(tmp_path: Path, m):
    out = tmp_path / "ok.txt"

    # All lists empty — means everything succeeded
    m.write_results(
        output_file=out,
        failed=[],
        cancelled=[],
        unknown=[],
        resolved_failures=[],
    )

    txt = out.read_text()
    assert "SUCCESSFUL" in txt


# ---------- CLI: check_failed_jobs ----------


def test_cli_no_logs(monkeypatch, m, tmp_path: Path, caplog):
    caplog.set_level(logging.WARNING)
    # Force find_job_logs -> {}
    monkeypatch.setattr(m, "find_job_logs", lambda p: {})
    runner = CliRunner()
    out_file = tmp_path / "out.txt"
    result = runner.invoke(
        m.check_failed_jobs,
        [str(tmp_path), "-o", str(out_file), "--log-level", "INFO"],
    )
    # Should not crash, should warn, and not create output
    assert result.exit_code == 0
    assert any("No job logs found" in rec.message for rec in caplog.records)
    assert not out_file.exists()


def test_cli_classifies_and_writes(monkeypatch, m, tmp_path: Path, caplog):
    caplog.set_level(logging.DEBUG)

    # Provide a fake mapping {jobid: path}
    fake_map = {
        "100": Path("/fake/100.log"),
        "200": Path("/fake/200.log"),
        "300": Path("/fake/300.log"),
        "400": Path("/fake/400.log"),
    }
    monkeypatch.setattr(m, "find_job_logs", lambda p: fake_map)

    # get_job_state behavior per jobid – already normalized states
    def fake_get(jobid: str):
        table = {
            "100": "FAILED",
            "200": "CANCELLED",
            "300": "COMPLETED",
            "400": None,  # e.g., missing sacct info
        }
        return table[jobid]

    monkeypatch.setattr(m, "get_job_state", fake_get)

    out = tmp_path / "report.txt"
    runner = CliRunner()
    result = runner.invoke(
        m.check_failed_jobs,
        [str(tmp_path), "-o", str(out), "--log-level", "INFO"],
    )
    assert result.exit_code == 0
    txt = out.read_text()

    # We now expect the "resolved failures" section, because the rule dir /fake
    # has at least one COMPLETED job (300), and 100/200 are failure-like states.
    assert "Some jobs failed but succeeded on retry:" in txt
    assert "100\t/fake/100.log\tFAILED" in txt
    assert "200\t/fake/200.log\tCANCELLED" in txt

    # Unknown section and job 400
    assert "Unknown status jobs:" in txt
    assert "400\tNA" in txt

    # No unresolved failed/cancelled jobs -> SUCCESS marker present
    assert "SUCCESSFUL" in txt

    # Completed job should not appear at all in error/resolution sections
    assert "300" not in txt
