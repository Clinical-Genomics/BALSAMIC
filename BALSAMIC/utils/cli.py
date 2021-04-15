import os
import json
import shutil
import logging
import sys
import collections
import re
import subprocess
from pathlib import Path
from io import StringIO
from distutils.spawn import find_executable

import yaml
import snakemake
import graphviz
from colorclass import Color

import BALSAMIC
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


class CaptureStdout(list):
    """
    Captures stdout.
    """

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
    quiet           - Quiet mode for snakemake
    singularity_arg - Singularity arguments to pass to snakemake
    sm_opt          - snakemake additional options
    disable_variant_caller - Disable variant caller
    dragen          - enable/disable dragen suite
    slurm_profiler  - enable slurm profiler
    """

    def __init__(self):
        self.case_name = str()
        self.working_dir = str()
        self.snakefile = str()
        self.configfile = str()
        self.run_mode = str()
        self.profile = str()
        self.cluster_config = str()
        self.scheduler = str()
        self.log_path = str()
        self.script_path = str()
        self.result_path = str()
        self.qos = str()
        self.account = str()
        self.mail_type = str()
        self.mail_user = str()
        self.forceall = False
        self.run_analysis = False
        self.quiet = False
        self.report = str()
        self.use_singularity = True
        self.singularity_bind = str()
        self.singularity_arg = str()
        self.sm_opt = str()
        self.disable_variant_caller = str()
        self.dragen = False
        self.slurm_profiler = str()

    def build_cmd(self):
        forceall = str()
        quiet_mode = str()
        sm_opt = str()
        cluster_cmd = str()
        dryrun = str()
        report = str()
        snakemake_config_key_value = list()

        if self.forceall:
            forceall = "--forceall"

        if self.report:
            report = "--report {}".format(self.report)

        if self.quiet:
            quiet_mode = " --quiet "

        if self.sm_opt:
            sm_opt = " ".join(self.sm_opt)

        if not self.run_analysis:
            dryrun = "--dryrun"

        if self.disable_variant_caller:
            snakemake_config_key_value.append(
                f'disable_variant_caller={self.disable_variant_caller}')

        if self.dragen:
            snakemake_config_key_value.append('dragen=True')

        if snakemake_config_key_value:
            snakemake_config_key_value.insert(0, "--config")

        if self.use_singularity:
            self.singularity_arg = "--use-singularity --singularity-args ' --cleanenv "
            for bind_path in self.singularity_bind:
                self.singularity_arg += " --bind {}:{}".format(
                    bind_path, bind_path)
            self.singularity_arg += "' "

        if self.run_mode == "cluster":
            sbatch_cmd = (" '{} {} "
                          " --sample-config {} --profile {} "
                          " --account {} --qos {} "
                          " --log-dir {} --script-dir {} "
                          " --result-dir {} ".format(
                              sys.executable,
                              self.scheduler,
                              self.configfile,
                              self.profile,
                              self.account,
                              self.qos,
                              self.log_path,
                              self.script_path,
                              self.result_path,
                          ))

            if self.slurm_profiler:
                sbatch_cmd += " --slurm-profiler {}".format(
                    self.slurm_profiler)

            if self.mail_user:
                sbatch_cmd += " --mail-user {} ".format(self.mail_user)

            if self.mail_type:
                sbatch_cmd += " --mail-type {} ".format(self.mail_type)

            sbatch_cmd += " {dependencies} '"

            cluster_cmd = (" --immediate-submit -j 999 "
                           "--jobname BALSAMIC.{}.{{rulename}}.{{jobid}}.sh "
                           "--cluster-config {} --cluster {} ".format(
                               self.case_name, self.cluster_config,
                               sbatch_cmd))

        # Merge snakmake config key value list
        snakemake_config_key_value = " ".join(snakemake_config_key_value)

        sm_cmd = (
            f" snakemake --notemp -p "
            f" --directory {self.working_dir} --snakefile {self.snakefile} --configfiles {self.configfile} "
            f" {self.cluster_config} {self.singularity_arg} {quiet_mode} "
            f" {forceall} {dryrun} {cluster_cmd} "
            f" {report} {snakemake_config_key_value} {sm_opt}")
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
    """
    Creates directories by recursively checking if it exists,
    otherwise increments the number
    """
    if os.path.isdir(os.path.abspath(path)):
        basepath = os.path.basename(os.path.abspath(path))
        basepath_number = 0
        if "." in basepath:
            basepath_number = int(basepath.split(".")[1])
        basepath_string = basepath.split(".")[0]
        basepath_number += 1
        path = os.path.join(
            os.path.dirname(os.path.abspath(path)),
            ".".join([basepath_string, str(basepath_number)]),
        )
        interm_path.append(path)
        createDir(path, interm_path)
        return os.path.abspath(interm_path[-1])
    else:
        os.makedirs(os.path.abspath(path), exist_ok=True)
        return os.path.abspath(path)


def write_json(json_out, output_config):

    try:
        with open(output_config, "w") as fn:
            json.dump(json_out, fn, indent=4)
    except OSError as error:
        raise error


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
    scheduler = str(Path(p, "utils", "scheduler.py"))

    return scheduler


def get_snakefile(analysis_type, sequencing_type="targeted"):
    """
    Return a string path for variant calling snakefile.
    """

    p = Path(__file__).parents[1]
    snakefile = Path(p, "workflows", "balsamic.smk")
    if analysis_type == "generate_ref":
        snakefile = Path(p, 'workflows', 'reference.smk')
    elif analysis_type == "umi":
        snakefile = Path(p, "workflows", "UMIworkflow.smk")

    return str(snakefile)


def get_config(config_name):
    """
    Return a string path for config file.
    """

    p = Path(__file__).parents[1]
    config_file = str(Path(p, "config", config_name + ".json"))
    if Path(config_file).exists():
        return config_file
    else:
        raise FileNotFoundError(f"Config for {config_name} was not found.")


def recursive_default_dict():
    """
    Recursivly create defaultdict.
    """
    return collections.defaultdict(recursive_default_dict)


def convert_defaultdict_to_regular_dict(inputdict: dict):
    """
    Recursively convert defaultdict to dict.
    """
    if isinstance(inputdict, collections.defaultdict):
        inputdict = {
            key: convert_defaultdict_to_regular_dict(value)
            for key, value in inputdict.items()
        }
    return inputdict


def find_file_index(file_path):
    indexible_files = {
        ".bam": [".bam.bai", ".bai"],
        ".cram": [".cram.crai", ".crai"],
        ".vcf.gz": [".vcf.gz.tbi"],
        ".vcf": [".vcf.tbi"],
    }

    file_path_index = set()
    for file_extension, file_index_extensions in indexible_files.items():
        if file_path.endswith(file_extension):
            for file_index_extension in file_index_extensions:
                new_file_path = file_path.replace(file_extension,
                                                  file_index_extension)
                if os.path.isfile(new_file_path):
                    file_path_index.add(new_file_path)

    return list(file_path_index)


def get_file_extension(file_path):
    known_multi_extensions = [
        ".vcf.gz", ".vcf.gz.tbi", ".vcf.tbi", ".fastq.gz"
    ]
    file_extension = ""
    for known_ext in known_multi_extensions:
        if file_path.endswith(known_ext):
            file_extension = known_ext
            break

    if not file_extension:
        _, file_extension = os.path.splitext(file_path)

    return file_extension[1:]


def get_from_two_key(input_dict, from_key, by_key, by_value, default=None):
    """
    Given two keys with list of values of same length, find matching index of by_value in from_key from by_key.
    
    from_key and by_key should both exist
    """

    matching_value = default
    if from_key in input_dict and by_key in input_dict and by_value in input_dict[
            from_key]:
        idx = input_dict[from_key].index(by_value)
        matching_value = input_dict[by_key][idx]

    return matching_value


def get_file_status_string(file_to_check):
    """
      Checks if file exsits. and returns a string with checkmark or redcorss mark
      if it exists or doesn't exist respectively.
      Always assume file doesn't exist, unless proven otherwise.
    """
    return_str = Color("[{red}\u2717{/red}] File missing: ") + file_to_check

    file_status = os.path.isfile(file_to_check)
    if file_status:
        return_str = Color("[{green}\u2713{/green}] Found: ") + file_to_check

    return return_str, file_status


def singularity(sif_path: str, cmd: str, bind_paths: list) -> str:
    """Run within container

    Excutes input command string via Singularity container image

    Args:
        sif_path: Path to singularity image file (sif)
        cmd: A string for series of commands to run
        bind_path: a path to bind within container

    Returns:
        A sanitized Singularity cmd

    Raises:
        BalsamicError: An error occured while creating cmd
    """

    singularity_cmd = shutil.which("singularity")
    if not singularity_cmd:
        raise BalsamicError("singularity command does not exist")

    if not Path(sif_path).is_file():
        raise BalsamicError("container file does not exist")

    singularity_bind_path = ""
    for bind_path in bind_paths:
        singularity_bind_path += "--bind {} ".format(bind_path)

    shellcmd = "singularity exec {} {}".format(singularity_bind_path, cmd)

    return " ".join(shellcmd.split())


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
    """Finds the correct filename prefix from file path, and returns it.
    An error is raised if sample name has invalid pattern """

    fq_pattern = re.compile(r"R_[12]" + ".fastq.gz$")
    sample_basename = Path(sample).name

    file_str = sample_basename[0:(
        fq_pattern.search(sample_basename).span()[0] + 1)]
    return file_str


def get_panel_chrom(panel_bed) -> list:
    """Returns a set of chromosomes present in PANEL BED"""

    lines = [line.rstrip("\n") for line in open(panel_bed, "r")]
    return {s.split("\t")[0] for s in lines}


def get_bioinfo_tools_version(bioinfo_tools: dict,
                              container_conda_env_path: os.PathLike) -> dict:
    """Parses the names and versions of bioinfo tools 
    used by BALSAMIC from config YAML into a dict """

    bioinfo_tools_version = {}
    for container_conda_env_name in set(bioinfo_tools.values()):
        yaml_file = Path(container_conda_env_path, container_conda_env_name,
                         container_conda_env_name + ".yaml")
        with open(yaml_file, "r") as f:
            packages = yaml.safe_load(f).get("dependencies")
            for p in packages:
                name = p.split("=")[0]
                version = "=".join(p.split("=")[1:])
                if name not in bioinfo_tools:
                    continue
                if name in bioinfo_tools_version:
                    bioinfo_tools_version[name].append(version)
                    bioinfo_tools_version[name] = list(
                        set(bioinfo_tools_version[name]))
                else:
                    bioinfo_tools_version[name] = list([version])
    return bioinfo_tools_version


def get_sample_dict(tumor: str,
                    normal: str,
                    tumor_sample_name: str = None,
                    normal_sample_name: str = None) -> dict:
    """Concatenates sample dicts for all provided files"""
    samples = {}
    if normal:
        for sample in normal:
            key, val = get_sample_names(sample, "normal")
            samples[key] = val
            samples[key]["sample_name"] = normal_sample_name

    for sample in tumor:
        key, val = get_sample_names(sample, "tumor")
        samples[key] = val
        samples[key]["sample_name"] = tumor_sample_name
    return samples


def get_sample_names(filename, sample_type):
    """Creates a dict with sample prefix, sample type, and readpair suffix"""
    file_str = validate_fastq_pattern(filename)
    if file_str:
        return (
            file_str,
            {
                "file_prefix": file_str,
                "type": sample_type,
                "readpair_suffix": ["1", "2"],
            },
        )


def create_fastq_symlink(casefiles, symlink_dir: Path):
    """Creates symlinks for provided files in analysis/fastq directory.
    Identifies file prefix pattern, and also creates symlinks for the 
    second read file, if needed"""

    for filename in casefiles:
        parent_dir = Path(filename).parents[0]
        file_str = validate_fastq_pattern(filename)
        for f in parent_dir.rglob(f"*{file_str}*.fastq.gz"):
            try:
                LOG.info(
                    f"Creating symlink {f} -> {Path(symlink_dir, f.name)}")
                Path(symlink_dir, f.name).symlink_to(f)
            except FileExistsError:
                LOG.info(f"File {symlink_dir / f.name} exists, skipping")


def generate_graph(config_collection_dict, config_path):
    """Generate DAG graph using snakemake stdout output"""

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(
            snakefile=get_snakefile(
                analysis_type=config_collection_dict["analysis"]
                ["analysis_type"],
                sequencing_type=config_collection_dict["analysis"]
                ["sequencing_type"],
            ),
            dryrun=True,
            configfiles=[config_path],
            printrulegraph=True,
        )

    graph_title = "_".join([
        "BALSAMIC",
        BALSAMIC.__version__,
        config_collection_dict["analysis"]["case_id"],
    ])
    graph_dot = "".join(graph_dot).replace(
        "snakemake_dag {",
        'BALSAMIC { label="' + graph_title + '";labelloc="t";')
    graph_obj = graphviz.Source(
        graph_dot,
        filename=".".join(
            config_collection_dict["analysis"]["dag"].split(".")[:-1]),
        format="pdf",
        engine="dot",
    )
    graph_obj.render(cleanup=True)


def get_fastq_bind_path(fastq_path: Path) -> list():
    """Takes a path with symlinked fastq files. 
    Returns unique paths to parent directories for singulatiry bind
    """
    parents = set()
    for fastq_file_path in Path(fastq_path).iterdir():
        parents.add(Path(fastq_file_path).resolve().parent.as_posix())
    return list(parents)


def convert_deliverables_tags(delivery_json: dict,
                              sample_config_dict: dict) -> dict:
    """Replaces values of file_prefix with sample_name in deliverables dict"""

    for file in delivery_json["files"]:
        file_tags = file["tag"].split(",")
        for sample in sample_config_dict["samples"]:
            file_prefix = sample_config_dict["samples"][sample]["file_prefix"]
            sample_name = sample_config_dict["samples"][sample]["sample_name"]
            sample_type = sample_config_dict["samples"][sample]["type"]
            if file_prefix == file["id"]:# or sample_type == file["id"]:
                file["id"] = sample_name
                for tag_index, tag in enumerate(file_tags):
                    if tag == file_prefix or tag == file_prefix.replace(
                            "_", "-"):
                        file_tags[tag_index] = sample_name
                if sample_name not in file_tags:
                    file_tags.append(sample_name)
            if sample_type == file["id"]:
                file["id"] = sample_name
                if sample_name not in file_tags:
                    file_tags.append(sample_name)
        file["tag"] = list(set(file_tags))
    return delivery_json


def check_executable(exe_name: str) -> bool:
    """Checks for executable exe_name in PATH"""
    exe_exist = True

    if find_executable(exe_name) is None:
        exe_exist = False

    return exe_exist


def generate_h5(job_name: str, job_id: str, file_path: str) -> str:
    """Generates H5 file for a finished job. Returns None if it cannot generate H5 file"""
    h5_file_name = Path(file_path, job_name + ".h5")
    sh5util_output = subprocess.check_output(
        ["sh5util", "-o",
         h5_file_name.as_posix(), "-S", "-j", job_id],
        stderr=subprocess.STDOUT)

    if "sh5util: No node-step files found for jobid" in sh5util_output.decode(
            "utf-8"):
        h5_file_name = None

    return h5_file_name


def job_id_dump_to_yaml(job_id_dump: Path, job_id_yaml: Path, case_name: str):
    """Write an input job_id_sacct_file to yaml output"""
    with open(job_id_dump, "r") as jobid_in, open(job_id_yaml,
                                                  "w") as jobid_out:
        jobid_list = jobid_in.read().splitlines()
        yaml.dump({case_name: jobid_list}, jobid_out)
