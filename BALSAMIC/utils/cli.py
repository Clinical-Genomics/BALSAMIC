import logging
import os
import re
import subprocess
import sys
import zlib
from distutils.spawn import find_executable
from io import StringIO
from pathlib import Path
from typing import Dict, Optional, List

import graphviz
import snakemake
import yaml
from colorclass import Color

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.constants.cluster import ClusterConfigType
from BALSAMIC.constants.paths import CONSTANTS_DIR

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
                f"disable_variant_caller={self.disable_variant_caller}"
            )

        if self.dragen:
            snakemake_config_key_value.append("dragen=True")

        if snakemake_config_key_value:
            snakemake_config_key_value.insert(0, "--config")

        if self.use_singularity:
            self.singularity_arg = "--use-singularity --singularity-args ' --cleanenv "
            for bind_path in self.singularity_bind:
                self.singularity_arg += " --bind {}:{}".format(bind_path, bind_path)
            self.singularity_arg += "' "

        if self.run_mode == "cluster":
            sbatch_cmd = (
                " '{} {} "
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
                )
            )

            if self.slurm_profiler:
                sbatch_cmd += " --slurm-profiler {}".format(self.slurm_profiler)

            if self.mail_user:
                sbatch_cmd += " --mail-user {} ".format(self.mail_user)

            if self.mail_type:
                sbatch_cmd += " --mail-type {} ".format(self.mail_type)

            sbatch_cmd += " {dependencies} '"

            cluster_cmd = (
                " --immediate-submit -j 999 "
                "--jobname BALSAMIC.{}.{{rulename}}.{{jobid}}.sh "
                "--cluster-config {} --cluster {} ".format(
                    self.case_name, self.cluster_config, sbatch_cmd
                )
            )

        # Merge snakmake config key value list
        snakemake_config_key_value = " ".join(snakemake_config_key_value)

        sm_cmd = (
            f" snakemake --notemp -p "
            f" --directory {self.working_dir} --snakefile {self.snakefile} --configfiles {self.configfile} "
            f" {self.cluster_config} {self.singularity_arg} {quiet_mode} "
            f" {forceall} {dryrun} {cluster_cmd} "
            f" {report} {snakemake_config_key_value} {sm_opt}"
        )
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


def iterdict(dic):
    """dictionary iteration - returns generator"""
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


def get_snakefile(analysis_type, analysis_workflow="balsamic", reference_genome="hg19"):
    """
    Return a string path for variant calling snakefile.
    """

    p = Path(__file__).parents[1]
    snakefile = Path(p, "workflows", "balsamic.smk")

    if analysis_type == "generate_ref":
        snakefile = Path(p, "workflows", "reference.smk")
        if "canfam3" in reference_genome:
            snakefile = Path(p, "workflows", "reference-canfam3.smk")
            return str(snakefile)

    if analysis_type == "pon":
        snakefile = Path(p, "workflows", "PON.smk")

    if "balsamic-qc" in analysis_workflow:
        snakefile = Path(p, "workflows", "QC.smk")

    return str(snakefile)


def get_config_path(config_type: ClusterConfigType) -> Path:
    """Return a config path given its type."""
    return Path(CONSTANTS_DIR, config_type + ".json")


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
                new_file_path = file_path.replace(file_extension, file_index_extension)
                if os.path.isfile(new_file_path):
                    file_path_index.add(new_file_path)

    return list(file_path_index)


def get_file_extension(file_path):
    known_multi_extensions = [".vcf.gz", ".vcf.gz.tbi", ".vcf.tbi", ".fastq.gz"]
    file_extension = ""
    for known_ext in known_multi_extensions:
        if file_path.endswith(known_ext):
            file_extension = known_ext
            break

    if not file_extension:
        _, file_extension = os.path.splitext(file_path)

    return file_extension[1:]


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


def get_panel_chrom(panel_bed) -> list:
    """Returns a set of chromosomes present in PANEL BED"""

    lines = [line.rstrip("\n") for line in open(panel_bed, "r")]
    return {s.split("\t")[0] for s in lines}


def bioinfo_tool_version_non_conda(packages: dict, bioinfo_tools: dict) -> dict:
    """Parses a non-conda bioinfo tool dictionary

    Args:
        packages: dict. Dictionary of bioinfo tools in tool_name=version format
        bioinfo_tools: dict. Bioinfo tools with key as tool name

    Returns:
        non_conda_bioinfo_version: dict. A dictionary of tools as key and version as list in value
    """
    bioinfo_version = {}
    for p in packages:
        name = p.split("=")[0]
        version = "=".join(p.split("=")[1:])
        if name in bioinfo_tools:
            bioinfo_version[name] = list([version])

    return bioinfo_version


def bioinfo_tool_version_conda(
    packages: list, bioinfo_tools: dict, current_bioinfo_tool_version: dict
) -> dict:
    """Parses conda environment dictionary and extracts dependencies as version

    Args:
        packages: list. List of bioinfo tools in channel::tool_name=version format
        bioinfo_tools: dict. Bioinfo tools with key as tool name
        current_bioinfo_tool_version: dict. Current list of bioinfo tool versions.

    Returns:
        conda_bioinfo_version: dict. A dictionary of tools as key and version as list in value
    """
    conda_bioinfo_version = current_bioinfo_tool_version
    for package in packages:
        if isinstance(package, dict):
            # Extraction of pip specific packages
            bioinfo_tool_version_conda(
                package[list(package.keys())[0]],
                bioinfo_tools,
                current_bioinfo_tool_version,
            )
            continue
        name_version: List[str] = list(
            filter(None, re.split("=|==", package))
        )  # Supporting PIP versioning (double equal sign)
        name, version = name_version[0], name_version[1]
        if name not in bioinfo_tools:
            continue
        if name in conda_bioinfo_version:
            conda_bioinfo_version[name].append(version)
            conda_bioinfo_version[name] = list(set(conda_bioinfo_version[name]))
        else:
            conda_bioinfo_version[name] = list([version])

    return conda_bioinfo_version


def get_bioinfo_tools_version(
    bioinfo_tools: dict, container_conda_env_path: os.PathLike
) -> dict:
    """Parses the names and versions of bioinfo tools
    used by BALSAMIC from config YAML into a dict.

    Args:
        bioinfo_tools: dict. A dictionary with bioinfo tools as keys and container as value
        container_conda_env_path: path. Path to all container and conda yaml

    Returns:
        bioinfo_tools_version: dict. Dictionary of bioinfo tools as key and version as value
    """

    bioinfo_tools_version = {}
    for container_conda_env_name in set(bioinfo_tools.values()):
        yaml_file = Path(
            container_conda_env_path,
            container_conda_env_name,
            container_conda_env_name + ".yaml",
        )
        with open(yaml_file, "r") as f:
            conda_yaml = yaml.safe_load(f)
            if isinstance(conda_yaml, dict):
                packages = conda_yaml.get("dependencies")
                bioinfo_tools_version = {
                    **bioinfo_tools_version,
                    **bioinfo_tool_version_conda(
                        packages, bioinfo_tools, bioinfo_tools_version
                    ),
                }
            else:
                bioinfo_tools_version = {
                    **bioinfo_tools_version,
                    **bioinfo_tool_version_non_conda(conda_yaml, bioinfo_tools),
                }
    return bioinfo_tools_version


def get_sample_dict(
    tumor_sample_name: str, normal_sample_name: Optional[str]
) -> Dict[str, dict]:
    """Returns a sample dictionary given the names of the tumor and/or normal samples."""
    sample_dict: Dict[str, dict] = {tumor_sample_name: {"type": "tumor"}}
    if normal_sample_name:
        sample_dict.update({normal_sample_name: {"type": "normal"}})
    return sample_dict


def get_pon_sample_dict(directory: str) -> Dict[str, dict]:
    """Returns a PON sample dictionary."""
    sample_dict: Dict[str, dict] = {}
    for file in Path(directory).glob("*.fastq.gz"):
        sample_name: str = file.name.split("_")[-4]
        sample_dict.update({sample_name: {"type": "normal"}})
    return sample_dict


def generate_graph(config_collection_dict, config_path):
    """Generate DAG graph using snakemake stdout output."""

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(
            snakefile=get_snakefile(
                analysis_type=config_collection_dict["analysis"]["analysis_type"],
                analysis_workflow=config_collection_dict["analysis"][
                    "analysis_workflow"
                ],
                reference_genome=config_collection_dict["reference"][
                    "reference_genome"
                ],
            ),
            dryrun=True,
            configfiles=[config_path],
            printrulegraph=True,
        )

    graph_title = "_".join(
        [
            "BALSAMIC",
            balsamic_version,
            config_collection_dict["analysis"]["case_id"],
        ]
    )
    graph_dot = "".join(graph_dot).replace(
        "snakemake_dag {", 'BALSAMIC { label="' + graph_title + '";labelloc="t";'
    )
    graph_obj = graphviz.Source(
        graph_dot,
        filename=".".join(config_collection_dict["analysis"]["dag"].split(".")[:-1]),
        format="pdf",
        engine="dot",
    )
    graph_obj.render(cleanup=True)


def convert_deliverables_tags(delivery_json: dict, sample_config_dict: dict) -> dict:
    """Replaces values of sample_type with sample_name in deliverables dict."""

    for delivery_file in delivery_json["files"]:
        file_tags = delivery_file["tag"].split(",")
        for sample in sample_config_dict["samples"]:
            sample_type = sample_config_dict["samples"][sample]["type"]
            if sample == delivery_file["id"]:
                for tag_index, tag in enumerate(file_tags):
                    if tag == sample or tag == sample.replace("_", "-"):
                        file_tags[tag_index] = sample
                if sample not in file_tags:
                    file_tags.append(sample)
            if sample_type == delivery_file["id"]:
                delivery_file["id"] = sample
            if sample_type in file_tags:
                file_tags.remove(sample_type)
        delivery_file["tag"] = list(set(file_tags))
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
        ["sh5util", "-o", h5_file_name.as_posix(), "-S", "-j", job_id],
        stderr=subprocess.STDOUT,
    )

    if "sh5util: No node-step files found for jobid" in sh5util_output.decode("utf-8"):
        h5_file_name = None

    return h5_file_name


def job_id_dump_to_yaml(job_id_dump: Path, job_id_yaml: Path, case_name: str):
    """Write an input job_id_sacct_file to yaml output"""
    with open(job_id_dump, "r") as jobid_in, open(job_id_yaml, "w") as jobid_out:
        jobid_list = jobid_in.read().splitlines()
        yaml.dump({case_name: jobid_list}, jobid_out)


def get_md5(filename):
    with open(filename, "rb") as fh:
        hashed = 0
        while True:
            s = fh.read(65536)
            if not s:
                break
            hashed = zlib.crc32(s, hashed)
    return "%08X" % (hashed & 0xFFFFFFFF)


def create_md5(reference, check_md5):
    """create a md5 file for all reference data"""
    with open(check_md5, "w") as fh:
        for key, value in reference.items():
            if os.path.isfile(value):
                fh.write(get_md5(value) + " " + value + "\n")


def get_resolved_fastq_files_directory(directory: str) -> str:
    """Return the absolute path for the directory containing the input fastq files."""
    input_files: List[Path] = [
        file.absolute() for file in Path(directory).glob("*.fastq.gz")
    ]
    if not input_files or not input_files[0].is_symlink():
        return directory
    return os.path.commonpath([file.resolve().as_posix() for file in input_files])


def get_analysis_fastq_files_directory(case_dir: str, fastq_path: str) -> str:
    """Return analysis fastq directory, linking the fastq files if necessary."""
    analysis_fastq_path: Path = Path(case_dir, "fastq")
    analysis_fastq_path.mkdir(parents=True, exist_ok=True)
    if Path(case_dir) not in Path(fastq_path).parents:
        for fastq in Path(fastq_path).glob("*.fastq.gz"):
            try:
                Path(analysis_fastq_path, fastq.name).symlink_to(fastq)
                LOG.info(f"Created link for {fastq} in {analysis_fastq_path}")
            except FileExistsError:
                LOG.warning(
                    f"File {Path(analysis_fastq_path, fastq.name)} exists. Skipping linking."
                )

        return analysis_fastq_path.as_posix()
    return fastq_path
