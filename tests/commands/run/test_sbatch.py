

from BALSAMIC.

def test_sbatch():
    # GIVEN values for sbatch command
    sbatch_cmd = sbatch()
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
