import os
import json
import yaml
import sys
import collections
import BALSAMIC
import snakemake
import re
import shutil
import logging
import click
import graphviz


from pathlib import Path
from colorclass import Color
from io import StringIO
from itertools import chain
from collections import defaultdict
from BALSAMIC.utils.constants import CONDA_ENV_PATH

LOG = logging.getLogger(__name__)

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
    
    case_name       - analysis case name
    working_dir     - working directory for snakemake
    configfile      - sample configuration file (json) output of balsamic-config-sample
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
        self.profile = None
        self.cluster_config = str() 
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
                self.singularity_arg += " --bind {}:{}".format(
                    bind_path, bind_path)
            self.singularity_arg += "' "

        if self.run_mode == 'cluster':
            sbatch_cmd = " '{} {} ".format(sys.executable, self.scheduler) + \
                " --sample-config " + self.configfile + \
                " --profile " + self.profile + \
                " --account " + self.account + \
                " --qos " + self.qos + \
                " --log-dir " + self.log_path + \
                " --script-dir " + self.script_path + \
                " --result-dir " + self.result_path

            if self.mail_user:
                sbatch_cmd += " --mail-user " + self.mail_user

            if self.mail_type:
                sbatch_cmd += " --mail-type " + self.mail_type

            sbatch_cmd += " {dependencies} '"

            cluster_cmd = " --immediate-submit -j 999 " + \
                " --jobname BALSAMIC." + self.case_name + ".{rulename}.{jobid}.sh" + \
                " --cluster-config " + self.cluster_config + \
                " --cluster " + sbatch_cmd

        sm_cmd = " snakemake --notemp -p " + \
            " --directory " + self.working_dir + \
            " --snakefile " + self.snakefile + \
            " --configfiles " + self.configfile + " " + self.cluster_config + \
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


def get_schedulerpy():
    """
    Returns a string path for scheduler.py
    """

    p = Path(__file__).parents[1]
    scheduler = str(Path(p, 'commands/run/scheduler.py'))

    return scheduler


def get_snakefile(analysis_type, sequencing_type="targeted"):
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

def recursive_default_dict():
    '''
    Recursivly create defaultdict.
    '''
    return collections.defaultdict(recursive_default_dict)


def convert_defaultdict_to_regular_dict(inputdict: dict):
    '''
    Recursively convert defaultdict to dict.
    '''
    if isinstance(inputdict, collections.defaultdict):
        inputdict = {
            key: convert_defaultdict_to_regular_dict(value)
            for key, value in inputdict.items()
        }
    return inputdict


def merge_dict_on_key(dict_1, dict_2, by_key):
    '''
    Merge two list of dictionaries based on key
    '''
    merged_dict = defaultdict(dict)
    for interm_list in (dict_1, dict_2):
        for item in interm_list:
            merged_dict[item[by_key]].update(item)
    merged_dict_list = merged_dict.values()
    return merged_dict_list


def find_file_index(file_path):
    indexible_files = {
        ".bam": [".bam.bai", ".bai"],
        ".cram": [".cram.cai", ".cai"],
        ".vcf.gz": [".vcf.gz.tbi"],
        ".vcf": [".vcf.tbi"],
    }

    file_path_index = set()
    for file_extension, file_index_extensions in indexible_files.items():
        if file_path.endswith(file_extension):
            for file_index_extension in file_index_extensions:
                new_file_path = file_path.replace(
                    file_extension, file_index_extension
                )
                if os.path.isfile(new_file_path):
                    file_path_index.add(new_file_path)

    return list(file_path_index)

def get_file_extension(file_path):
    known_multi_extensions = ['.vcf.gz', '.vcf.gz.tbi', '.vcf.tbi', '.fastq.gz']
    file_extension = ""
    for known_ext in known_multi_extensions:
        if file_path.endswith(known_ext):
            file_extension = known_ext
            break

    if not file_extension:
        file_name, file_extension = os.path.splitext(file_path)

    return file_extension

def get_from_two_key(input_dict, from_key, by_key, by_value, default=None):
    '''
    Given two keys with list of values of same length, find matching index of by_value in from_key from by_key.
    
    from_key and by_key should both exist
    '''

    matching_value = default
    if from_key in input_dict and by_key in input_dict and by_value in input_dict[from_key]:
            idx = input_dict[from_key].index(by_value)
            matching_value = input_dict[by_key][idx]

    return matching_value


def get_file_status_string(file_to_check):
    """
      Checks if file exsits. and returns a string with checkmark or redcorss mark
      if it exists or doesn't exist respectively.
      Always assume file doesn't exist, unless proven otherwise.
    """
    return_str = Color(u"[{red}\u2717{/red}] File missing: ") + file_to_check

    file_status = os.path.isfile(file_to_check)
    if file_status: 
        return_str = Color(u"[{green}\u2713{/green}] Found: ") + file_to_check
    
    return return_str, file_status

def merge_json(*args):
    """
    Take a list of json files and merges them together
    Input: list of json file
    Output: dictionary of merged json
    """

    json_out = dict()
    for json_file in args:
        try:
            if isinstance(json_file, dict):
                json_out = {**json_out, **json_file}
            else:
                with open(json_file) as fn:
                    json_out = {**json_out, **json.load(fn)}
        except OSError as error:
            raise error

    return json_out

def validate_fastq_pattern(sample):
    fq_pattern = re.compile(r"R_[12]" + ".fastq.gz$")
    sample_basename = Path(sample).name

    file_str = sample_basename[0:(
        fq_pattern.search(sample_basename).span()[0] + 1)]
    return file_str


def get_panel_chrom(panel_bed) -> list:
    lines = [line.rstrip('\n') for line in open(panel_bed, 'r')]
    return {s.split('\t')[0] for s in lines}


def get_bioinfo_tools_list(conda_env_path) -> dict:
    bioinfo_tools = {}
    for yaml_file in Path(conda_env_path).rglob('*.yaml'):
        with open(yaml_file, "r") as f:
            packages = yaml.safe_load(f).get("dependencies")
            for p in packages:
                try:
                    name, version = p.split("=")
                except ValueError:
                    name, version = p, None
                finally:
                    bioinfo_tools[name] = version
    return bioinfo_tools


def get_sample_dict(tumor, normal) -> dict:
    samples = {}
    if normal:
        for sample in normal:
            key, val = get_sample_names(sample, "normal")
            samples[key] = val

    for sample in tumor:
        key, val = get_sample_names(sample, "tumor")
        samples[key] = val
    return samples


def get_sample_names(file, sample_type):
    file_str = validate_fastq_pattern(file)
    if file_str:
        return file_str, {
            "file_prefix": file_str,
            "type": sample_type,
            "readpair_suffix": ["1", "2"]
        }


def create_fastq_symlink(casefiles, symlink_dir: Path):
    for filename in casefiles:
        parent_dir = Path(filename).parents[0]
        file_str = validate_fastq_pattern(filename)
        for f in parent_dir.rglob(f'*{file_str}*.fastq.gz'):
            try:
                LOG.info(
                    f"Creating symlink {f} -> {Path(symlink_dir, f.name)}")
                Path(symlink_dir, f.name).symlink_to(f)
            except FileExistsError:
                LOG.info(f"File {symlink_dir / f.name} exists, skipping")


def generate_graph(config_collection_dict, config_path):
    with CaptureStdout() as graph_dot:
        snakemake.snakemake(snakefile=get_snakefile(
            analysis_type=config_collection_dict["analysis"]["analysis_type"],
            sequencing_type=config_collection_dict["analysis"]
            ["sequencing_type"]),
                            dryrun=True,
                            configfiles=[config_path],
                            printrulegraph=True)

    graph_title = "_".join([
        'BALSAMIC', BALSAMIC.__version__,
        config_collection_dict["analysis"]["case_id"]
    ])
    graph_dot = "".join(graph_dot).replace(
        'snakemake_dag {',
        'BALSAMIC { label="' + graph_title + '";labelloc="t";')
    graph_obj = graphviz.Source(
        graph_dot,
        filename=".".join(
            config_collection_dict["analysis"]["dag"].split(".")[:-1]),
        format="pdf",
        engine="dot")
    graph_obj.render(cleanup=True)
