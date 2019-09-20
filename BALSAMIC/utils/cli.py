import os
import re
import json
import yaml
import sys
from io import StringIO
from pathlib import Path
from itertools import chain


class CaptureStdout(list):
    '''
    Captures stdout.
    '''

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


class SnakeMake:
    """
    To build a snakemake command using cli options

    Params:
    ------
    case_name       - analysis case name
    working_dir     - working directory for snakemake
    configfile      - sample configuration file (json) output of -
                      balsamic-config-sample
    run_mode        - run mode - cluster or local shell run
    cluster_config  - cluster config json file
    scheduler       - slurm command constructor
    log_path        - log file path
    script_path     - file path for slurm scripts
    result_path     - result directory
    qos             - QOS for sbatch jobs
    account         - scheduler(e.g. slurm) account
    mail_user       - email to account to send job run status
    forceall        - To add '--forceall' option for snakemake
    run_analysis    - To run pipeline
    use_singularity - To use singularity
    singularity_bind- Singularity bind path
    singularity_arg - Singularity arguments to pass to snakemake
    sm_opt          - snakemake additional options
    """

    def __init__(self):
        self.case_name = None
        self.working_dir = None
        self.snakefile = None
        self.configfile = None
        self.run_mode = None
        self.cluster_config = None
        self.scheduler = None
        self.log_path = None
        self.script_path = None
        self.result_path = None
        self.qos = None
        self.account = None
        self.mail_type = None
        self.mail_user = None
        self.forceall = False
        self.run_analysis = False
        self.use_singularity = True
        self.singularity_bind = None
        self.singularity_arg = str()
        self.sm_opt = None

    def build_cmd(self):
        forceall = ''
        sm_opt = ''
        cluster_cmd = ''
        dryrun = ''

        if self.forceall:
            forceall = " --forceall "

        if self.sm_opt:
            sm_opt = " ".join(self.sm_opt)

        if not self.run_analysis:
            dryrun = " --dryrun "

        if self.use_singularity:
            self.singularity_arg = " --use-singularity --singularity-args '"
            for bind_path in self.singularity_bind:
                self.singularity_arg += " --bind {}:{}".format(bind_path, bind_path)
            self.singularity_arg += "' "

        if self.run_mode == 'slurm':
            sbatch_cmd = " 'python3 {} ".format(self.scheduler) + \
                " --sample-config " + self.configfile + \
                " --slurm-account " + self.account + \
                " --slurm-qos " + self.qos + \
                " --log-dir " + self.log_path + \
                " --script-dir " + self.script_path + \
                " --result-dir " + self.result_path

            if self.mail_user:
                sbatch_cmd += " --slurm-mail-user " + self.mail_user

            if self.mail_type:
                sbatch_cmd += " --slurm-mail-type " + self.mail_type

            sbatch_cmd += " {dependencies} '"

            cluster_cmd = " --immediate-submit -j 999 " + \
                " --jobname BALSAMIC." + self.case_name + ".{rulename}.{jobid}.sh" + \
                " --cluster-config " + self.cluster_config + \
                " --cluster " + sbatch_cmd

        sm_cmd = " snakemake --notemp -p " + \
            " --directory " + self.working_dir + \
            " --snakefile " + self.snakefile + \
            " --configfile " + self.configfile + \
            self.singularity_arg + \
            " " + forceall + " " + dryrun + \
            " " + cluster_cmd + " " + sm_opt

        return sm_cmd


def add_doc(docstring):
    """
    A decorator for adding docstring. Taken shamelessly from stackexchange.
    """

    def document(func):
        func.__doc__ = docstring
        return func

    return document


def createDir(path, interm_path=[]):
    '''
    Creates directories by recursively checking if it exists,
    otherwise increments the number
    '''
    if os.path.isdir(os.path.abspath(path)):
        basepath = os.path.basename(os.path.abspath(path))
        basepath_number = 0
        if "." in basepath:
            basepath_number = int(basepath.split(".")[1])
        basepath_string = basepath.split(".")[0]
        basepath_number += 1
        path = os.path.join(os.path.dirname(os.path.abspath(path)),
                            ".".join([basepath_string,
                                      str(basepath_number)]))
        interm_path.append(path)
        createDir(path, interm_path)
        return os.path.abspath(interm_path[-1])
    else:
        os.makedirs(os.path.abspath(path), exist_ok=True)
        return os.path.abspath(path)


def get_packages(yaml_file):
    '''
    return packages found in a conda yaml file

    input: conda yaml file path
    output: list of packages
    '''
    try:
        with open(yaml_file, 'r') as f:
            pkgs = yaml.safe_load(f)['dependencies']
    except OSError as error:
        raise error

    return (pkgs)


def write_json(json_out, output_config):

    try:
        with open(output_config, "w") as fn:
            json.dump(json_out, fn, indent=4)
    except OSError as error:
        raise error


def get_package_split(condas):
    '''
    Get a list of conda env files, and extract pacakges

    input: conda env files
    output: dict of packages and their version
    '''

    pkgs = [
        "bwa", "bcftools", "cutadapt", "fastqc", "gatk", "manta", "picard",
        "sambamba", "strelka", "samtools", "tabix", "vardic"
    ]

    pkgs = dict(
        [[y.split("=")[0], y.split("=")[1]]
         for y in set(chain.from_iterable([get_packages(s) for s in condas]))
         if y.split("=")[0] in pkgs])

    return (pkgs)


def iterdict(dic):
    """ dictionary iteration - returns generator"""
    for key, value in dic.items():
        if isinstance(value, dict):
            yield from iterdict(value)
        else:
            yield key, value


def get_sbatchpy():
    """
    Returns a string path for runSbatch.py
    """

    p = Path(__file__).parents[1]
    sbatch = str(Path(p, 'commands/run/sbatch.py'))

    return sbatch


def get_snakefile(analysis_type, sequencing_type = "targeted"):
    """
    Return a string path for variant calling snakefile.
    """

    p = Path(__file__).parents[1]
    if analysis_type == "qc":
        snakefile = Path(p, 'workflows', 'Alignment')
    elif analysis_type in ["single", "paired"]:
        snakefile = Path(p, 'workflows', 'VariantCalling')
        if sequencing_type == "wgs":
            snakefile = Path(p, 'workflows', 'VariantCalling_sentieon')
    elif analysis_type == "generate_ref":
        snakefile = Path(p, 'workflows', 'GenerateRef')

    return str(snakefile)


def get_config(config_name):
    """
    Return a string path for config file.
    """

    p = Path(__file__).parents[1]
    config_file = str(Path(p, 'config', config_name + ".json"))
    if Path(config_file).exists():
        return config_file
    else:
        raise FileNotFoundError(f'Config for {config_name} was not found.')


def get_ref_path(input_json):
    """
    Set full path to reference files
    Input: reference config file
    Return: json file with abspath
    """
    with open(input_json) as fh:
        ref_json = json.load(fh)
        for k, v in ref_json['reference'].items():
            ref_json['reference'][k] = os.path.abspath(v)

    return ref_json
