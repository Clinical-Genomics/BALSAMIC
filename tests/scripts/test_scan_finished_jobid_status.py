import importlib
import types
from pathlib import Path
from click.testing import CliRunner
import subprocess
import textwrap
import re
import logging
import pytest

# 👇 Change this to your actual module filename (no .py)
MODULE_NAME = "BALSAMIC.assets.scripts.scan_finished_jobid_status"


@pytest.fixture(scope="module")
def m():
    # Import the target module once for all tests
    return importlib.import_module(MODULE_NAME)


# ---------- parse_state ----------


def test_parse_state_ok(m):
    out = "JobId=1234 JobState=FAILED Reason=NonZeroExitCode"
    assert m.parse_state(out) == "FAILED"


def test_parse_state_missing(m, caplog):
    caplog.set_level(logging.DEBUG)
    out = "JobId=5678 SomeOtherField=foo"
    assert m.parse_state(out) is None
    assert any("JobState not found" in rec.message for rec in caplog.records)


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


# ---------- get_job_state ----------


def test_get_job_state_success(monkeypatch, m):
    class Dummy:
        stdout = "JobId=42 JobState=COMPLETED"

    def fake_run(*args, **kwargs):
        return Dummy()

    monkeypatch.setattr(subprocess, "run", fake_run)

    out = m.get_job_state("42")
    assert "JobState=COMPLETED" in out


def test_get_job_state_not_found(monkeypatch, m, caplog):
    caplog.set_level(logging.ERROR)

    def boom(*args, **kwargs):
        raise FileNotFoundError("no scontrol")

    monkeypatch.setattr(subprocess, "run", boom)

    assert m.get_job_state("7") is None
    assert any("scontrol executable not found" in rec.message for rec in caplog.records)


def test_get_job_state_calledprocesserror(monkeypatch, m, caplog):
    caplog.set_level(logging.DEBUG)

    def boom(*args, **kwargs):
        raise subprocess.CalledProcessError(
            returncode=1, cmd=["scontrol"], stderr="bad things"
        )

    monkeypatch.setattr(subprocess, "run", boom)

    assert m.get_job_state("99") is None
    assert any("Could not check job 99" in rec.message for rec in caplog.records)
    assert any(
        "scontrol stderr for 99 bad things" in rec.message for rec in caplog.records
    )


# ---------- write_results ----------


def test_write_results_all_sections(tmp_path: Path, m):
    out = tmp_path / "result.txt"
    failed = [("101", Path("/logs/101.log")), ("202", Path("/logs/202.log"))]
    cancelled = [("303", Path("/logs/303.log"))]
    unknown = ["404", "505"]

    m.write_results(out, failed, cancelled, unknown)
    txt = out.read_text()

    # Timestamp header present (we don't assert the exact time)
    assert "=== Job status check at " in txt

    # Failed block
    assert "Failed jobs:" in txt
    assert "101\t/logs/101.log" in txt
    assert "202\t/logs/202.log" in txt

    # Cancelled block
    assert "Cancelled jobs:" in txt
    assert "303\t/logs/303.log" in txt

    # Unknown block
    assert "Unknown status jobs:" in txt
    assert "404\tNA" in txt and "505\tNA" in txt

    # No SUCCESSFUL marker because there are failures/cancellations
    assert "SUCCESSFUL" not in txt


def test_write_results_success_only(tmp_path: Path, m):
    out = tmp_path / "ok.txt"
    m.write_results(out, failed=[], cancelled=[], unknown=[])
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

    # get_job_state behavior per jobid
    def fake_get(jobid: str):
        table = {
            "100": "JobId=100 JobState=FAILED",
            "200": "JobId=200 JobState=CANCELLED",
            "300": "JobId=300 JobState=COMPLETED",
            "400": None,  # e.g., missing/failed scontrol
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

    # Present sections
    assert "Failed jobs:" in txt
    assert "100\t/fake/100.log" in txt

    assert "Cancelled jobs:" in txt
    assert "200\t/fake/200.log" in txt

    assert "Unknown status jobs:" in txt
    assert "400\tNA" in txt

    # Completed job should not appear in any section
    assert "300" not in re.findall(r"^(?:Failed|Cancelled|Unknown).*$", txt, flags=re.M)
