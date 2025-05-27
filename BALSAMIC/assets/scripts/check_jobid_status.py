import click
import re
import subprocess
from pathlib import Path


@click.command()
@click.argument(
    "log_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "--output",
    "-o",
    type=click.Path(writable=True, path_type=Path),
    help="Path to output file for failed jobs.",
)
def check_failed_jobs(log_dir, output):
    """
    Scan LOG_DIR for SLURM log files (*.out, *.err), extract job IDs,
    and check if any jobs have failed using `scontrol show job JOBID`.

    If --output is provided, results are written to a file.
    """
    jobid_pattern = re.compile(r"\.(\d+)\.(?:out|err)$")
    stderr_pattern = re.compile(r"StdErr=(\S+)")
    job_ids = set()

    # Extract job IDs from filenames
    for file in log_dir.iterdir():
        match = jobid_pattern.search(file.name)
        if match:
            job_ids.add(match.group(1))

    if not job_ids:
        click.echo("No job IDs found in log filenames.")
        return

    failed_jobs_info = []

    # Check SLURM job state
    for jobid in sorted(job_ids):
        try:
            result = subprocess.run(
                ["scontrol", "show", "job", jobid],
                capture_output=True,
                text=True,
                check=True,
            )
            output_text = result.stdout
            if "JobState=FAILED" in output_text:
                stderr_match = stderr_pattern.search(output_text)
                stderr_path = stderr_match.group(1) if stderr_match else "N/A"
                failed_jobs_info.append((jobid, stderr_path))
        except subprocess.CalledProcessError:
            click.echo(f"Could not check job {jobid} (may not exist).")

    # Output results
    if output:
        with output.open("w") as out_f:
            if failed_jobs_info:
                out_f.write("Failed jobs:\n")
                for jid, stderr in failed_jobs_info:
                    out_f.write(f"{jid}\t{stderr}\n")
            else:
                click.echo("All jobs completed successfully.\n")
        click.echo(f"Results written to {output}")
    else:
        if failed_jobs_info:
            click.echo("Failed jobs:")
            for jid, stderr in failed_jobs_info:
                click.echo(f"  - {jid}\t{stderr}")
        else:
            click.echo("All jobs completed successfully.")


if __name__ == "__main__":
    check_failed_jobs()
