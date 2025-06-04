import click
import re
import subprocess
from pathlib import Path


def extract_job_ids(log_dir):
    jobid_pattern = re.compile(r"\.(\d+)\.(?:out|err)$")
    return {
        match.group(1)
        for file in log_dir.iterdir()
        if (match := jobid_pattern.search(file.name))
    }


def get_job_state(jobid):
    try:
        result = subprocess.run(
            ["/usr/bin/scontrol", "show", "job", jobid],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout
    except subprocess.CalledProcessError:
        click.echo(f"Could not check job {jobid} (may not exist).")
        return None


def extract_stderr_path(output_text):
    match = re.search(r"StdErr=(\S+)", output_text)
    return match.group(1) if match else "N/A"


def categorize_job(jobid, output_text, failed, cancelled):
    stderr_path = extract_stderr_path(output_text)
    if "JobState=FAILED" in output_text:
        failed.append((jobid, stderr_path))
    elif "JobState=CANCELLED" in output_text:
        cancelled.append((jobid, stderr_path))


def write_results(output_file, failed, cancelled):
    with output_file.open("a") as out_f:
        if failed:
            out_f.write("Failed jobs:\n")
            for jobid, stderr in failed:
                out_f.write(f"{jobid}\t{stderr}\n")
        elif cancelled:
            out_f.write("Cancelled jobs:\n")
            for jobid, stderr in cancelled:
                out_f.write(f"{jobid}\t{stderr}\n")
        else:
            click.echo("All jobs completed successfully.")
    click.echo(f"Results written to {output_file}")


@click.command()
@click.argument(
    "log_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(writable=True, path_type=Path),
    help="Path to output file for failed jobs.",
)
def check_failed_jobs(log_dir: Path, output: Path):
    """
    Scan LOG_DIR for SLURM log files (*.out, *.err), extract job IDs,
    and check if any jobs have failed or been cancelled using `scontrol show job JOBID`.

    If --output is provided, results are written to a file.
    """
    job_ids = extract_job_ids(log_dir)

    if not job_ids:
        click.echo("No job IDs found in log filenames.")
        return

    failed_jobs_info = []
    cancelled_jobs_info = []

    for jobid in sorted(job_ids):
        output_text = get_job_state(jobid)
        if output_text:
            categorize_job(jobid, output_text, failed_jobs_info, cancelled_jobs_info)

    write_results(output, failed_jobs_info, cancelled_jobs_info)


if __name__ == "__main__":
    check_failed_jobs()
