#!/usr/bin/env python3

from BALSAMIC.commands.run.sbatch import SbatchScheduler, QsubScheduler


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

    # WHEN sbatch command is built
    sbatch_cmd = sbatch_cmd.build_cmd()

    # THEN sbatch command string is constructed
    assert isinstance(sbatch_cmd, str)
    assert sbatch_cmd == 'sbatch --account "development" --dependency "afterok:12345" --error "test_job.err" --output "test_job.out" --mail-type "FAIL" --mail-user "john.doe@example.com" --ntasks "2" --qos "low" --time "01:00:00" example_script.sh'


def test_qsub_scheduler():
    # GIVEN values for qsub command
    qsub_cmd = QsubScheduler()
    qsub_cmd.account = "development"
    qsub_cmd.dependency = ['12345']
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
    assert qsub_cmd == 'qsub  -A development  -e test_job.err  -o test_job.out  -m FAIL  -M john.doe@example.com  -p low  -l "walltime=01:00:00,nodes=1:ppn=2"   -W "depend=afterok:12345"  example_script.sh '
