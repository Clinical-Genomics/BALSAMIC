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
        out_f.write(f"=== Job status check at {timestamp} ===\n")

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
                "(jobid\tlog_path\toriginal_state\n"
            )
            for jobid, log_path, state, rule_key in resolved_failures:
                out_f.write(f"{jobid}\t{log_path}\t{state}\n")
            out_f.write("\n")

        if unknown:
            out_f.write("Unknown status jobs:\n")
            for jobid in unknown:
                out_f.write(f"{jobid}\tNA\n")
            out_f.write("\n")

        # If there are no *unresolved* failures/cancellations, consider run successful.
        if not failed and not cancelled:
            out_f.write("SUCCESSFUL (no unresolved failed/cancelled jobs)\n\n")

    LOG.info(
        "Appended results to %s (failed=%d, cancelled=%d, resolved_failures=%d, unknown=%d)",
        output_file,
        len(failed),
        len(cancelled),
        len(resolved_failures),
        len(unknown),
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

    unknown: List[str] = []

    # rule_key -> list of (jobid, log_path, state)
    rule_to_jobs: Dict[str, List[Tuple[str, Path, Optional[str]]]] = {}

    if not job_logs:
        LOG.warning("No job logs found (no files matching '*.log')")
        return

    for jobid in sorted(job_logs.keys(), key=int):
        log_path = job_logs[jobid]

        out_text = get_job_state(jobid)
        if not out_text:
            LOG.warning(
                "Missing scontrol output for job %s -- setting status UNKNOWN", jobid
            )
            unknown.append(jobid)
            continue

        state = parse_state(out_text)
        # Derive "rule key" from parent directory relative to log_dir.
        try:
            rule_key = str(log_path.parent.relative_to(log_dir))
        except ValueError:
            # Fallback if for some reason the path isn't under log_dir
            rule_key = str(log_path.parent)

        LOG.debug("Job %s in rule dir %s has state %s", jobid, rule_key, state)

        rule_to_jobs.setdefault(rule_key, []).append((jobid, log_path, state))

    failed: List[Tuple[str, Path]] = []
    cancelled: List[Tuple[str, Path]] = []
    resolved_failures: List[Tuple[str, Path, str, str]] = []

    failure_like_states = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY"}
    success_states = {
        "COMPLETED"
    }  # extend if you want to treat other states as success

    # Now classify based on per-rule success/failure
    for rule_key, jobs in rule_to_jobs.items():
        has_success = any(state in success_states for _, _, state in jobs)

        for jobid, log_path, state in jobs:
            if state not in failure_like_states:
                # Don't treat non-failure states as errors
                continue

            if has_success:
                # This rule eventually succeeded; earlier failures go to the "resolved" section
                resolved_failures.append((jobid, log_path, state, rule_key))
            else:
                # No successful retry in this rule dir; keep behavior as before
                if state == "FAILED":
                    failed.append((jobid, log_path))
                else:
                    cancelled.append((jobid, log_path))

    write_results(output, failed, cancelled, unknown, resolved_failures)


if __name__ == "__main__":
    check_failed_jobs()
