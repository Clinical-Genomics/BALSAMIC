#!/usr/bin/env python3
"""
Submit this clustering script for sbatch to snakemake with:
    snakemake --forceall --snakefile Snakefile \
      --cluster-config config_cluster.json \
      -j 9999 \
      --cluster 'python3 runSbatch.py --config config_cluster.json {dependencies}' \
      --immediate-submit --notemp
"""
import sys
import re
import os
import json
import argparse
from snakemake.utils import read_job_properties

parser = argparse.ArgumentParser(description='Process script and dependecies.')
parser.add_argument(
    "dependencies",
    nargs="*",
    help="{{dependencies}} string given by snakemake\n")
parser.add_argument(
    "snakescript",
    help="Snakemake generated shell script with commands to execute snakemake rule\n")
# parser.add_argument(
#    "--cluster-config",
#    help="Config file to read sbatch settings from.")
parser.add_argument(
    "--sample-config",
    help="Config file to read sbatch settings from.")
args = parser.parse_args()

sample_config = json.load(open(args.sample_config))
jobscript = args.snakescript
job_properties = read_job_properties(jobscript)

time = job_properties["cluster"]["time"]
cpu = job_properties["cluster"]["n"]
account_slurm = job_properties["cluster"]["account"]
mail_type = job_properties["cluster"]["mail_type"]
mail_user = job_properties["cluster"]["mail_user"]

logpath = os.path.join(
    sample_config["analysis"]["analysis_dir"],
    sample_config["analysis"]["sample_id"],
    sample_config["analysis"]["log"])
scriptpath = os.path.join(
    sample_config["analysis"]["analysis_dir"],
    sample_config["analysis"]["sample_id"],
    sample_config["analysis"]["script"])
resultpath = os.path.join(
    sample_config["analysis"]["analysis_dir"],
    sample_config["analysis"]["sample_id"],
    sample_config["analysis"]["result"])

os.system('mkdir -p ' + logpath)
os.system('mkdir -p ' + scriptpath)
os.system('mkdir -p ' + resultpath)

os.system('cp ' + jobscript + ' ' + scriptpath)
scriptname = jobscript.split("/")
scriptname = scriptname[-1]
jobscript = os.path.join(scriptpath, scriptname)

output_log = logpath + scriptname + ".out"
error_log = logpath + scriptname + ".err"
cmdline = 'sbatch -A {account} -n {n} -t {time} --qos=low -o {output_log} -e {error_log} --mail-type {mail_type} --mail-user {mail_user}'.format(
    n=cpu, time=time, output_log=output_log, error_log=error_log, account=account_slurm, mail_type=mail_type, mail_user=mail_user)

cmdline += " "

dependencies = args.dependencies
if dependencies:
    cmdline += '--dependency=' + \
        ','.join(["afterok:%s" % d for d in dependencies])

cmdline += " " + jobscript + " | cut -d' ' -f 4"

os.system(cmdline)
