import subprocess
import json
import pytest

from unittest import mock

from BALSAMIC.commands.run.scheduler import SbatchScheduler
from BALSAMIC.commands.run.scheduler import QsubScheduler
from BALSAMIC.commands.run.scheduler import read_sample_config
from BALSAMIC.commands.run.scheduler import write_sacct_file
from BALSAMIC.commands.run.scheduler import submit_job
from BALSAMIC.commands.run.scheduler import main as scheduler_main
from BALSAMIC.utils.cli import createDir


def test_scheduler_slurm_py(snakemake_job_script, tumor_normal_config, capsys):
    # GIVEN a jobscript, dependencies, joutput job id, and sample comamnd
    test_jobid = '999999999999'
    test_return_value = 'Submitted batch job ' + test_jobid
    scheduler_args = [
        '9000', '9001', '9002', snakemake_job_script['snakescript']
    ]
    scheduler_profile_slurm = 'slurm'
    with open(tumor_normal_config, 'r') as input_config:
        sample_config = json.load(input_config)

    # Create directory for log and script
    script_dir = createDir(sample_config['analysis']['script'])
    log_dir = createDir(sample_config['analysis']['log'])

    # Construct scheduler's cmd
    scheduler_cmd = [
        "--sample-config", tumor_normal_config, "--profile",
        scheduler_profile_slurm, "--qos", "low", "--account", "development",
        "--log-dir", log_dir, "--script-dir", script_dir, "--result-dir",
        sample_config['analysis']['result']
    ]
    scheduler_cmd.extend(scheduler_args)

    # WHEN calling scheduler_main with mocked subprocess
    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = test_return_value.encode('utf-8')
        scheduler_main(scheduler_cmd)

    # THEN sacct file should be written with the job id(s)
    with open(log_dir + '/sample_tumor_normal.sacct', 'r') as fin:
        assert fin.read() == test_jobid + "\n"

    # THEN captured output is job id
    captured = capsys.readouterr()
    assert captured.out == test_jobid + "\n"


def test_scheduler_qsub_py(snakemake_job_script, tumor_normal_config, capsys):
    # GIVEN a jobscript, dependencies, joutput job id, and sample comamnd
    test_jobname = 'script.sh'
    test_return_value = f'Your job 31415 ("{test_jobname}") has been submitted'
    scheduler_args = [
        '1000', '1001', '1002', snakemake_job_script['snakescript']
    ]
    scheduler_profile_qsub = 'qsub'
    with open(tumor_normal_config, 'r') as input_config:
        sample_config = json.load(input_config)

    # Create directory for log and script
    script_dir = createDir(sample_config['analysis']['script'])
    log_dir = createDir(sample_config['analysis']['log'])

    # Construct scheduler's cmd
    scheduler_cmd = [
        "--sample-config", tumor_normal_config, "--profile",
        scheduler_profile_qsub, "--qos", "low", "--account", "development",
        "--log-dir", log_dir, "--script-dir", script_dir, "--result-dir",
        sample_config['analysis']['result']
    ]
    scheduler_cmd.extend(scheduler_args)

    # WHEN calling scheduler_main with mocked subprocess
    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = test_return_value.encode('utf-8')
        scheduler_main(scheduler_cmd)

    # THEN sacct file should be written with the job id(s)
    with open(log_dir + '/sample_tumor_normal.sacct', 'r') as fin:
        assert fin.read() == test_jobname + "\n"

    # THEN captured output is job id
    captured = capsys.readouterr()
    assert captured.out == test_jobname + "\n"


def test_submit_job_slurm(snakemake_job_script):
    # GIVEN a jobid
    test_jobid = '1234'
    test_return_value = 'Submitted batch job ' + test_jobid

    # WHEN getting jobid for slurm
    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = test_return_value.encode('utf-8')
        actual_jobid = submit_job(['random_command'], 'slurm')

    # THEN output jobid should match
    assert actual_jobid == test_jobid


def test_submit_job_qsub(snakemake_job_script):
    # GIVEN a jobid
    test_jobname = 'script.sh'
    test_return_value = f'Your job 31415 ("{test_jobname}") has been submitted'

    # WHEN getting jobid for slurm
    with mock.patch.object(subprocess, 'run') as mocked:
        mocked.return_value.stdout = test_return_value.encode('utf-8')
        actual_jobname = submit_job(['random_command'], 'qsub')

    # THEN output jobid should match
    assert actual_jobname == test_jobname


def test_SbatchScheduler():
    # GIVEN values for sbatch command
    sbatch_cmd = SbatchScheduler()
    sbatch_cmd.account = "development"
    sbatch_cmd.dependency = "afterok:12345"
    sbatch_cmd.error = "test_job.err"
    sbatch_cmd.output = "test_job.out"
    sbatch_cmd.mail_type = "FAIL"
    sbatch_cmd.mail_user = "john.doe@example.com"
    sbatch_cmd.ntasks = "2"
    sbatch_cmd.qos = "low"
    sbatch_cmd.time = "01:00:00"
    sbatch_cmd.script = "example_script.sh"
    sbatch_cmd.partition = "dummy_partition"

    # WHEN sbatch command is built
    sbatch_cmd = sbatch_cmd.build_cmd()

    # THEN sbatch command string is constructed
    assert isinstance(sbatch_cmd, str)
    assert sbatch_cmd == (
        'sbatch --account "development" --dependency "afterok:12345" --error "test_job.err" '
        '--output "test_job.out" --mail-type "FAIL" --mail-user "john.doe@example.com" '
        '--ntasks "2" --qos "low" --time "01:00:00" --partition dummy_partition example_script.sh'
    )


def test_qsub_scheduler():
    # GIVEN values for qsub command
    qsub_cmd = QsubScheduler()
    qsub_cmd.account = "development"
    qsub_cmd.dependency = ['test_jobname.sh']
    qsub_cmd.error = "test_job.err"
    qsub_cmd.output = "test_job.out"
    qsub_cmd.mail_type = "FAIL"
    qsub_cmd.mail_user = "john.doe@example.com"
    qsub_cmd.ntasks = "2"
    qsub_cmd.qos = "low"
    qsub_cmd.time = "01:00:00"
    qsub_cmd.script = "example_script.sh"

    # WHEN qsub command is built
    qsub_cmd = qsub_cmd.build_cmd()

    # THEN qsub command should be constructed
    assert isinstance(qsub_cmd, str)
    assert qsub_cmd == (
        'qsub -V -S /bin/bash  -q development  -e test_job.err  -o test_job.out  -m s   -M '
        'john.doe@example.com  -p low  -l excl=1  -pe mpi 2   -hold_jid test_jobname.sh  example_script.sh '
    )


def test_read_sample_config_err(config_files):
    with pytest.raises(Exception):
        # GIVEN a bed file instead of json file
        bed_file = config_files['panel_bed_file']

        # WHEN calling read_sample_config
        # THEN It should raise the exception error
        assert read_sample_config(bed_file)


def test_write_sacct_file_err():
    with pytest.raises(FileNotFoundError):
        # GIVEN a non-existing file path and jobid
        dummy_file_path = "dummy/dummy_fname"
        dummy_jobid = "12345"

        # WHEN calling write_sacct_file
        # THEN It should raise the exception
        assert write_sacct_file(dummy_file_path, dummy_jobid)


def test_submit_job_err():
    with pytest.raises(subprocess.CalledProcessError):
        # GIVEN a wrong command
        sbatch_cmd = "SBATCH jobscript.sh"
        profile = 'slurm'

        # WHEN calling submit_job function
        # THEN it should return the exit code 1 and raise the subprocess error
        assert submit_job(sbatch_cmd, profile)
