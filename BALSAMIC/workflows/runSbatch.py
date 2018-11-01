#!/usr/bin/env python3
import sys
import re
import os
import subprocess
import json
import argparse
from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(
    description=
    'This is an internal script and should be invoked independently. This script gets a list of arugments (see --help) and submits a job to slurm as afterok dependency. All inputs are mandatory!'
)
parser.add_argument(
    "dependencies",
    nargs="*",
    help="{{dependencies}} string given by snakemake\n")
parser.add_argument(
    "snakescript",
    help=
    "Snakemake generated shell script with commands to execute snakemake rule\n"
)
parser.add_argument(
    "--sample-config", help="Config file to read sbatch settings from.")
parser.add_argument("--dir-log", help="Log directory")
parser.add_argument("--dir-result", help="Result directory")
parser.add_argument("--dir-script", help="Script directory")
parser.add_argument(
    "--qos", default="low", help="QOS for sbatch jobs. [Default: low]")
args = parser.parse_args()

sample_config = json.load(open(args.sample_config))

jobscript = args.snakescript
job_properties = read_job_properties(jobscript)
logpath = args.dir_log
scriptpath = args.dir_script
resultpath = args.dir_result

time = job_properties["cluster"]["time"]
cpu = job_properties["cluster"]["n"]
account_slurm = job_properties["cluster"]["account"]
mail_type = job_properties["cluster"]["mail_type"]
mail_user = job_properties["cluster"]["mail_user"]

subprocess.call('cp ' + jobscript + ' ' + scriptpath + '/', shell=True)

scriptname = jobscript.split("/")
scriptname = scriptname[-1]
jobscript = os.path.join(scriptpath, scriptname)
sacct_file = os.path.join(logpath, sample_config["analysis"]["sample_id"] + ".sacct")

output_log = os.path.join(logpath, scriptname + "_%j.out")
error_log = os.path.join(logpath, scriptname + "_%j.err")
cmdline = 'sbatch -A {account} -n {n} -t {time} --qos={qos} -o {output_log} -e {error_log} --mail-type {mail_type} --mail-user {mail_user}'.format(
    n=cpu,
    time=time,
    qos=args.qos,
    output_log=output_log,
    error_log=error_log,
    account=account_slurm,
    mail_type=mail_type,
    mail_user=mail_user)

cmdline += " "

dependencies = args.dependencies
if dependencies:
    cmdline += '--dependency=' + \
        ','.join(["afterok:%s" % d for d in dependencies])

cmdline += " " + jobscript + " | cut -d' ' -f 4"

cmdline += " " + " >> " + sacct_file

subprocess.call(cmdline, shell=True)
subprocess.call("tail -n1 " +  sacct_file, shell=True)
