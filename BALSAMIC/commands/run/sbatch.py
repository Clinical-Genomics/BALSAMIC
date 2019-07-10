#!/usr/bin/env python3
import sys
import os
import re
import subprocess
import json
import argparse
import shutil
from BALSAMIC.utils.cli import sbatch as sbatch_cmd
from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(description='''
    This is an internal script and should be invoked independently.
    This script gets a list of arugments (see --help) and submits a job to slurm as
    afterok dependency.
    ''')
parser.add_argument("dependencies",
                    nargs="*",
                    help="{{dependencies}} from snakemake")
parser.add_argument("snakescript", help="Snakemake script")
parser.add_argument("--sample-config", help="balsamic config sample output")
parser.add_argument('--slurm-account', required=True, help='SLURM account name')
parser.add_argument('--slurm-qos', default='low', help='SLURM job QOS')
parser.add_argument('--slurm-mail-type', help='SLURM mail type')
parser.add_argument('--slurm-mail-user', help='SLURM mail user')
parser.add_argument("--dir-log", help="Log directory")
parser.add_argument("--dir-result", help="Result directory")
parser.add_argument("--dir-script", help="Script directory")

args = parser.parse_args()

with open(args.sample_config) as f:
  sample_config = json.load(f)

logpath = args.dir_log
scriptpath = args.dir_script
resultpath = args.dir_result
jobscript = args.snakescript
job_properties = read_job_properties(jobscript)

shutil.copy2(jobscript, scriptpath)
jobscript = os.path.join(scriptpath, os.path.basename(jobscript))
log_error = os.path.join(logpath, os.path.basename(jobscript) + "_%j.err")
log_output = os.path.join(logpath, os.path.basename(jobscript) + "_%j.out")

sacct_file = os.path.join(logpath, sample_config["analysis"]["sample_id"] + ".sacct")

time = job_properties["cluster"]["time"]
cpu = job_properties["cluster"]["n"]

if not args.slurm_mail_type:
    mail_type = job_properties["cluster"]["mail_type"]

balsamic_status = os.getenv("BALSAMIC_STATUS","conda")
if "BALSAMIC_STATUS" == "container":
  if "BALSAMIC_BIND_PATH" not in os.environ:
      raise ValueError("BALSAMIC_BIND_PATH environment variable was not found")
  else:
      bind_path = os.getenv("BALSAMIC_BIND_PATH")

  if "BALSAMIC_MAIN_ENV" not in os.environ:
      raise ValueError("BALSAMIC_MAIN_ENV environment variable was not found")
  else:
      main_env = os.getenv("BALSAMIC_MAIN_ENV")

  if "BALSAMIC_CONTAINER" not in os.environ:
      raise ValueError("BALSAMIC_CONTAINER environment variable was not found")
  else:
      container = os.getenv("BALSAMIC_CONTAINER")

  sbatch_script = os.path.join(scriptpath, "sbatch." + os.path.basename(jobscript))

  with open(sbatch_script, 'a') as f:
      f.write("#!/bin/bash" + "\n")
      if balsamic_status == "container":
         f.write(f"function balsamic_run {{ singularity exec -B {bind_path} --app {main_env} {container} $@; }}" + "\n")
         f.write(f"# Snakemake original script {jobscript}" + "\n")
         f.write(f"balsamic_run bash {jobscript}" + "\n")
  
  sbatch_file = os.path.join(logpath, sample_config["analysis"]["sample_id"] + ".sbatch")

sbatch = sbatch_cmd()
sbatch.account = args.slurm_account
if args.dependencies:
    sbatch.dependency = ','.join(["afterok:%s" % d for d in args.dependencies])
sbatch.error = log_error 
sbatch.output = log_output
sbatch.mail_type = mail_type
sbatch.mail_user = args.slurm_mail_user
sbatch.ntasks = cpu
sbatch.time = time

if "BALSAMIC_STATUS" == "container":
  sbatch.script = sbatch_script
else:
  sbatch.script = jobscript 

# run sbatch cmd
try:
    res = subprocess.run(sbatch.build_cmd(), check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

# Get jobid
res = res.stdout.decode()
try:
    m = re.search("Submitted batch job (\d+)", res)
    jobid = m.group(1)
    print(jobid)
except Exception as e:
    print(e)
    raise

if "BALSAMIC_STATUS" == "container":
  with open(sbatch_file, 'a') as f:
      f.write(sbatch.build_cmd() + "\n")
      f.write(sys.executable + "\n")

with open(sacct_file, 'a') as f:
    f.write(jobid + "\n")
