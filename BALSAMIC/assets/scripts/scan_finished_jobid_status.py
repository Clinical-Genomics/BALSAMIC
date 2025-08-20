#!/usr/bin/env python3
from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import click


LOG = logging.getLogger(__name__)


def find_job_logs(log_root: Path) -> Dict[str, Path]:
    """
    Recursively find *.log files whose basename is a numeric jobid.
    Returns {jobid -> log_path}.
    """
    job_logs: Dict[str, Path] = {}
    for p in log_root.rglob("*.log"):
        if p.stem.isdigit():  # e.g. "9727982.log" -> "9727982"
            job_logs[p.stem] = p
        else:
            LOG.debug("Skipping non-job log file: %s", p)
    LOG.info("Discovered %d job logs under %s", len(job_logs), log_root)
    return job_logs


def get_job_state(jobid: str, scontrol: str = "scontrol") -> Optional[str]:
    """
    Return raw output of `scontrol show job JOBID`, or None if the query fails.
    """
    try:
        LOG.debug("Running %s show job %s", scontrol, jobid)
        result = subprocess.run(
            [scontrol, "show", "job", jobid],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout
    except FileNotFoundError:
        LOG.error("scontrol executable not found: %s", scontrol)
        return None
    except subprocess.CalledProcessError as e:
        LOG.warning("Could not check job %s (may not exist). rc=%s", jobid, e.returncode)
        LOG.debug("scontrol stderr for %s: %s", jobid, e.stderr)
        return None


def parse_state(scontrol_output: str) -> Optional[str]:
    """
    Extract JobState from scontrol text, e.g. 'JobState=FAILED'.
    Returns the state string (e.g. 'FAILED') or None if not found.
    """
    m = re.search(r"JobState=(\S+)", scontrol_output)
    state = m.group(1) if m else None
    if state is None:
        LOG.debug("JobState not found in scontrol output")
    return state


def write_results(
    output_file: Path,
    failed: List[Tuple[str, Path]],
    cancelled: List[Tuple[str, Path]],
) -> None:
    """
    Write both failed and cancelled sections (if any). If neither exist, write SUCCESS.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open("w") as out_f:
        if failed:
            out_f.write("Failed jobs:\n")
            for jobid, log_path in failed:
                out_f.write(f"{jobid}\t{log_path}\n")
            out_f.write("\n")

        if cancelled:
            out_f.write("Cancelled jobs:\n")
            for jobid, log_path in cancelled:
                out_f.write(f"{jobid}\t{log_path}\n")
            out_f.write("\n")

        if not failed and not cancelled:
            out_f.write("SUCCESS\n")

    LOG.info(
        "Results written to %s (failed=%d, cancelled=%d)",
        output_file, len(failed), len(cancelled)
    )


@click.command()
@click.argument("log_dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(writable=True, path_type=Path),
    help="Path to output file for results (FAILED/CANCELLED or SUCCESS).",
)
@click.option(
    "--scontrol",
    default="scontrol",
    show_default=True,
    help="Path to the scontrol executable.",
)
@click.option(
    "--log-level",
    default="DEBUG",
    show_default=True,
    type=click.Choice(["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], case_sensitive=False),
    help="Logging verbosity.",
)
def check_failed_jobs(log_dir: Path, output: Path, scontrol: str, log_level: str) -> None:
    """
    Recursively scan LOG_DIR for SLURM *.log files (stdout+stderr combined),
    extract job IDs from filenames, and check their states via `scontrol show job JOBID`.
    """
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    LOG.info("Scanning logs under: %s", log_dir)
    job_logs = find_job_logs(log_dir)

    failed: List[Tuple[str, Path]] = []
    cancelled: List[Tuple[str, Path]] = []

    if not job_logs:
        LOG.warning("No job logs found (no files matching '*.log')")
        return

    for jobid in sorted(job_logs.keys(), key=int):
        out_text = get_job_state(jobid, scontrol=scontrol)
        if not out_text:
            # Can't classify without job info; skip but note it.
            LOG.debug("Skipping job %s due to missing scontrol output", jobid)
            continue

        state = parse_state(out_text)
        if state == "FAILED":
            failed.append((jobid, job_logs[jobid]))
        elif state == "CANCELLED":
            cancelled.append((jobid, job_logs[jobid]))
        else:
            LOG.debug("Job %s state is %s", jobid, state)

    write_results(output, failed, cancelled)


if __name__ == "__main__":
    check_failed_jobs()