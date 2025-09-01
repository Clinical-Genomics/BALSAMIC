#!/usr/bin/env python3
from __future__ import annotations

import logging
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime
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
            LOG.debug(f"Skipping non-job log file: {p}")
    LOG.info(f"Discovered {len(job_logs)} job logs under {log_root}")
    return job_logs


def get_job_state(jobid: str) -> Optional[str]:
    """
    Return raw output of `scontrol show job JOBID`, or None if the query fails.
    """
    try:
        LOG.debug(f"Running show job scontrol {jobid}")
        result = subprocess.run(
            ["scontrol", "show", "job", jobid],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout
    except FileNotFoundError:
        LOG.error("scontrol executable not found: scontrol")
        return None
    except subprocess.CalledProcessError as e:
        LOG.warning(f"Could not check job {jobid} (may not exist). rc={e.returncode}")
        LOG.debug(f"scontrol stderr for {jobid} {e.stderr}")
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
    unknown: List[str],
) -> None:
    """
    Append job results to output_file.
    Each run is prefixed with a timestamp header.
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with output_file.open("a") as out_f:
        out_f.write(f"=== Job status check at {timestamp} ===\n")

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

        if unknown:
            out_f.write("Unknown status jobs:\n")
            for jobid in unknown:
                out_f.write(f"{jobid}\tNA\n")
            out_f.write("\n")

        if not failed and not cancelled:
            out_f.write("SUCCESSFUL\n\n")

    LOG.info(
        f"Appended results to {output_file} (failed={len(failed)}, cancelled={len(cancelled)} unknown={len(unknown)})"
    )


@click.command()
@click.argument(
    "log_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(writable=True, path_type=Path),
    help="Path to output file for results (FAILED/CANCELLED or SUCCESS).",
)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(
        ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], case_sensitive=False
    ),
    help="Logging verbosity.",
)
def check_failed_jobs(log_dir: Path, output: Path, log_level: str) -> None:
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
    unknown: List[str] = []

    if not job_logs:
        LOG.warning("No job logs found (no files matching '*.log')")
        return

    for jobid in sorted(job_logs.keys(), key=int):
        out_text = get_job_state(jobid)
        if not out_text:
            # Can't classify without job info; skip but note it.
            LOG.warning(
                f"Missing scontrol output for job {jobid} -- setting status UNKNOWN"
            )
            unknown.append(jobid)
            continue

        state = parse_state(out_text)
        if state == "FAILED":
            failed.append((jobid, job_logs[jobid]))
        elif state == "CANCELLED":
            cancelled.append((jobid, job_logs[jobid]))
        else:
            LOG.debug(f"Job {jobid} state is {state}")

    write_results(output, failed, cancelled, unknown)


if __name__ == "__main__":
    check_failed_jobs()
