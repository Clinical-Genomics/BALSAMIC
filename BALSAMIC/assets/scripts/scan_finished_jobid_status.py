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

FAILURE_LIKE_STATES = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY"}
SUCCESS_STATES = {"COMPLETED"}


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
    Look up the final job state via `sacct`.

    We prefer the top-level job record (e.g. '10683002')
    over step records (e.g. '10683002.batch', '10683002.0').

    Returns a normalized state string like 'COMPLETED', 'FAILED',
    'CANCELLED', etc., or None if not found.
    """
    cmd = [
        "sacct",
        "-j",
        jobid,
        "--noheader",  # no column headers
        "--parsable2",  # '|' separator, stable columns
        "-o",
        "JobID,State",
    ]

    try:
        LOG.debug("Running: %s", " ".join(cmd))
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        LOG.error("sacct executable not found: sacct")
        return None
    except subprocess.CalledProcessError as e:
        LOG.warning("Could not check job %s via sacct (rc=%s)", jobid, e.returncode)
        LOG.debug("sacct stderr for %s: %s", jobid, e.stderr)
        return None

    lines = [ln.strip() for ln in result.stdout.splitlines() if ln.strip()]
    if not lines:
        LOG.debug("No sacct records returned for job %s", jobid)
        return None

    # Each line looks like: "10683002|FAILED" or "10683002.0|CANCELLED+"
    records = [ln.split("|") for ln in lines]

    # Prefer the exact jobid (no step suffix)
    parent_record = next((r for r in records if r[0] == jobid), None)
    chosen = parent_record or records[0]

    raw_state = chosen[1]
    # Normalize things like "CANCELLED+" or "FAILED node_fail"
    state = raw_state.split()[0].rstrip("+")
    LOG.debug("Job %s sacct raw state=%r -> normalized=%r", jobid, raw_state, state)

    return state


def write_results(
    output_file: Path,
    failed: List[Tuple[str, Path]],
    cancelled: List[Tuple[str, Path]],
    unknown: List[str],
    resolved_failures: List[Tuple[str, Path, str, str]],
) -> None:
    """
    Append job results to output_file.
    Each run is prefixed with a timestamp header.

    resolved_failures items are (jobid, log_path, state, rule_key).
    """
    output_file.parent.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with output_file.open("a") as out_f:
        out_f.write(f"=== Job status check at {timestamp} ===\n\n")

        # If there are no *unresolved* failures/cancellations, consider run successful.
        if not failed and not cancelled:
            out_f.write("SUCCESSFUL\n\n")

        if failed:
            out_f.write("FAILED JOBS (no successful retry):\n")
            for jobid, log_path in failed:
                out_f.write(f"{jobid}\t{log_path}\n")
            out_f.write("\n")

        if cancelled:
            out_f.write("CANCELLED JOBS (no successful retry):\n")
            for jobid, log_path in cancelled:
                out_f.write(f"{jobid}\t{log_path}\n")
            out_f.write("\n")

        if resolved_failures:
            out_f.write("\n")
            out_f.write(
                "NOTE:\n"
                "Some jobs failed but succeeded on retry:\n"
                "(jobid\tlog_path\toriginal_state)\n"
            )
            for jobid, log_path, state, rule_key in resolved_failures:
                out_f.write(f"{jobid}\t{log_path}\t{state}\n")
            out_f.write("\n")

        if unknown:
            out_f.write("Unknown status jobs:\n")
            for jobid in unknown:
                out_f.write(f"{jobid}\tNA\n")
            out_f.write("\n")

    LOG.info(
        "Appended results to %s (failed=%d, cancelled=%d, resolved_failures=%d, unknown=%d)",
        output_file,
        len(failed),
        len(cancelled),
        len(resolved_failures),
        len(unknown),
    )


def derive_rule_key(log_root: Path, log_path: Path) -> str:
    """
    Derive a "rule key" from the log path.

    Prefer path relative to log_root; fall back to absolute parent directory.
    """
    try:
        return str(log_path.parent.relative_to(log_root))
    except ValueError:
        return str(log_path.parent)


def group_jobs_by_rule(
    log_dir: Path, job_logs: Dict[str, Path]
) -> Tuple[Dict[str, List[Tuple[str, Path, Optional[str]]]], List[str]]:
    """
    Query scontrol for each job and group them per rule directory.

    Returns:
        rule_to_jobs: {rule_key -> [(jobid, log_path, state), ...]}
        unknown: list of jobids with missing scontrol info
    """
    rule_to_jobs: Dict[str, List[Tuple[str, Path, Optional[str]]]] = {}
    unknown: List[str] = []

    for jobid in sorted(job_logs.keys(), key=int):
        log_path = job_logs[jobid]

        state = get_job_state(jobid)
        if state is None:
            LOG.warning(
                "Missing sacct state for job %s -- setting status UNKNOWN", jobid
            )
            unknown.append(jobid)
            continue
        rule_key = derive_rule_key(log_dir, log_path)

        LOG.debug("Job %s in rule dir %s has state %s", jobid, rule_key, state)

        rule_to_jobs.setdefault(rule_key, []).append((jobid, log_path, state))

    return rule_to_jobs, unknown


def classify_jobs(
    rule_to_jobs: Dict[str, List[Tuple[str, Path, Optional[str]]]]
) -> Tuple[
    List[Tuple[str, Path]],
    List[Tuple[str, Path]],
    List[Tuple[str, Path, str, str]],
]:
    """
    Classify jobs into:
      - failed (no successful retry for that rule)
      - cancelled (no successful retry for that rule)
      - resolved_failures (failed/cancelled but with a successful retry)
    """
    failed: List[Tuple[str, Path]] = []
    cancelled: List[Tuple[str, Path]] = []
    resolved_failures: List[Tuple[str, Path, str, str]] = []

    for rule_key, jobs in rule_to_jobs.items():
        has_success = any(state in SUCCESS_STATES for _, _, state in jobs)

        for jobid, log_path, state in jobs:
            if state not in FAILURE_LIKE_STATES:
                continue

            if has_success:
                # Rule eventually succeeded; treat as resolved failure.
                resolved_failures.append((jobid, log_path, state, rule_key))
            elif state == "FAILED":
                failed.append((jobid, log_path))
            else:
                cancelled.append((jobid, log_path))

    return failed, cancelled, resolved_failures


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

    If multiple jobs share the same rule log directory and at least one of them
    completes successfully, earlier failures in that directory are reported under
    a separate heading as "Failed jobs with successful retry".
    """
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    LOG.info("Scanning logs under: %s", log_dir)
    job_logs = find_job_logs(log_dir)

    if not job_logs:
        LOG.warning("No job logs found (no files matching '*.log')")
        return

    rule_to_jobs, unknown = group_jobs_by_rule(log_dir, job_logs)
    failed, cancelled, resolved_failures = classify_jobs(rule_to_jobs)

    write_results(output, failed, cancelled, unknown, resolved_failures)


if __name__ == "__main__":
    check_failed_jobs()
