import os
import subprocess
import re
import json
import copy
import shutil
import glob
import logging
import click
import snakemake
import graphviz
import pkgutil
import sys
import yaml

from urllib.parse import urlparse
from datetime import datetime
from pathlib import Path
from yapf.yapflib.yapf_api import FormatFile
import BALSAMIC
from BALSAMIC.utils.cli import get_package_split
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_config
from BALSAMIC.utils.cli import get_ref_path
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.rule import get_chrom
from BALSAMIC import __version__ as bv

LOG = logging.getLogger(__name__)


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


def set_panel_bed(json_out, panel_bed):
    """
    Set panel path in config file
    """
    try:
        json_out["panel"] = {"capture_kit": os.path.abspath(panel_bed)}
        json_out["panel"]["chrom"] = get_chrom(panel_bed)

    except OSError as error:
        raise error

    return json_out


def check_exist(path):
    """
    Checks if fastq file readable and accessable.
    """

    try:
        f = open(path, "r")
        f.close()
    except (IOError, FileNotFoundError) as error:
        raise error

    return True


def get_analysis_type(normal, umi):
    """ return analysis type """

    return "paired" if normal else "single"


def get_output_config(config, case_id):
    """ return output config json file"""
    if not config:
        return case_id + "_" + datetime.now().strftime("%Y%m%d") + ".json"
    else:
        return config


def get_sample_config(sample_config, case_id, analysis_dir, analysis_type):
    """
    creating sample config to run the analysis
    """
    with open(sample_config) as sample_json:
        sample_config = json.load(sample_json)

    sample_config["analysis"]["case_id"] = case_id
    sample_config["analysis"]["config_creation_date"] = datetime.now(
    ).strftime("%Y-%m-%d %H:%M")
    sample_config["analysis"]["analysis_dir"] = analysis_dir + "/"
    sample_config["analysis"]["log"] = os.path.join(analysis_dir, case_id,
                                                    "logs/")
    sample_config["analysis"]["script"] = os.path.join(analysis_dir, case_id,
                                                       "scripts/")
    sample_config["analysis"]["result"] = os.path.join(analysis_dir, case_id,
                                                       "analysis")
    sample_config["analysis"]["benchmark"] = os.path.join(
        analysis_dir, case_id, "benchmarks/")
    sample_config["analysis"]["analysis_type"] = analysis_type
    sample_config["samples"] = {}

    return sample_config


def link_fastq(src_files, des_path):
    """
    Creating fastq copy in given destination
    """
    for src_file in src_files:
        basename = os.path.basename(src_file)
        des_file = os.path.join(des_path, basename)
        try:
            shutil.copyfile(Path(src_file).resolve(), des_file)
        except shutil.SameFileError as e:
            LOG.warning(e)
            LOG.warning(
                f"Desitination file {des_file} exists. No copy link was created."
            )


def get_fastq_path(fq_file, fq_pattern):
    # check the fastq if exists
    fq_file = os.path.abspath(fq_file)
    if Path(fq_file).exists():
        file_basename = os.path.basename(fq_file)
        try:
            # extracting file prefix
            file_str = file_basename[0:(
                fq_pattern.search(file_basename).span()[0] + 1)]
        except AttributeError as error:
            LOG.error(
                f"File name is invalid, fastq file should be sample_R_1.fastq.gz"
            )
            raise click.Abort()
    else:
        LOG.error(f"{fq_file} is not found, update correct file path")
        raise click.Abort()

    return file_str, os.path.split(fq_file)[0]


def configure_fastq(fq_path, sample, fastq_prefix):
    """
    Configure the fastq files for analysis
    """

    fq_pattern = re.compile(r"R_[12]" + fastq_prefix + ".fastq.gz$")
    paths = list()

    # get a list of fq files
    sample_str, sample_path = get_fastq_path(sample, fq_pattern)
    paths.append(sample_path)

    fq_files = set()
    for path in paths:
        for fq_file in os.listdir(path):
            if fq_pattern.search(fq_file):
                fq_files.add(os.path.join(path, fq_file))

    # create symlink
    link_fastq(fq_files, fq_path)

    # return file prefix
    return sample_str


def case_config0(
    context,
    umi,
    umi_trim_length,
    quality_trim,
    adapter_trim,
    reference_config,
    panel_bed,
    output_config,
    normal,
    tumor,
    case_id,
    analysis_dir,
    overwrite_config,
    create_dir,
    fastq_prefix,
    singularity,
):
    """
    Prepares a config file for balsamic run_analysis. For now it is just treating json as
    dictionary and merging them as it is. So this is just a placeholder for future.
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")

    config_path = Path(__file__).parents[2] / "config"
    config_path = config_path.absolute()

    balsamic_env = config_path / "balsamic_env.yaml"
    rule_directory = Path(__file__).parents[2]

    install_config = dict()

    install_config["conda_env_yaml"] = balsamic_env.as_posix()
    install_config["rule_directory"] = rule_directory.as_posix() + "/"

    install_config["singularity"] = dict()
    install_config["singularity"]["image"] = Path(
        singularity).absolute().as_posix()

    analysis_type = get_analysis_type(normal, umi)
    sequencing_type = "targeted" if panel_bed else "wgs"
    output_config = get_output_config(output_config, case_id)
    analysis_config = get_config("analysis")

    LOG.info("Reading analysis config file %s" % analysis_config)
    LOG.info("Reading reference config file %s" % reference_config)

    reference_json = get_ref_path(reference_config)

    read_prefix = ["1", "2"]

    sample_config_path = get_config("sample")

    LOG.info("Reading sample config file %s" % sample_config_path)

    analysis_dir = os.path.abspath(analysis_dir)
    sample_config = get_sample_config(sample_config_path, case_id,
                                      analysis_dir, analysis_type)

    output_dir = os.path.join(analysis_dir, case_id)

    if create_dir:
        os.makedirs(output_dir, exist_ok=True)

    # create dir for fastq symlink creation
    fq_path = os.path.join(output_dir, "analysis", "fastq")
    os.makedirs(fq_path, exist_ok=True)

    for t in tumor:
        t = configure_fastq(fq_path, t, fastq_prefix)

        sample_config["samples"][t] = {
            "file_prefix": t,
            "type": "tumor",
            "readpair_suffix": read_prefix,
        }

    if normal:
        normal = configure_fastq(fq_path, normal, fastq_prefix)
        sample_config["samples"][normal] = {
            "file_prefix": normal,
            "type": "normal",
            "readpair_suffix": read_prefix,
        }

    sample_config["analysis"]["fastq_path"] = fq_path + "/"
    sample_config["analysis"]["BALSAMIC_version"] = bv
    sample_config["analysis"]["sequencing_type"] = sequencing_type

    conda_env = glob.glob(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..",
                     "conda/*.yaml"))

    bioinfo_config = dict()
    bioinfo_config["bioinfo_tools"] = get_package_split(conda_env)

    output_config = os.path.join(output_dir, output_config)
    LOG.info("Writing output config file %s" % os.path.abspath(output_config))

    json_out = merge_json(analysis_config, sample_config, reference_json,
                          install_config, bioinfo_config)

    if umi:
        json_out["QC"]["umi_trim"] = umi
        json_out["QC"]["umi_trim_length"] = str(umi_trim_length)

    json_out["QC"]["quality_trim"] = quality_trim
    json_out["QC"]["adapter_trim"] = adapter_trim

    dag_image = os.path.join(output_dir,
                             output_config + "_BALSAMIC_" + bv + "_graph")

    json_out["analysis"]["dag"] = dag_image + ".pdf"

    if panel_bed:
        json_out = set_panel_bed(json_out, panel_bed)

    if overwrite_config:
        write_json(json_out, output_config)

    FormatFile(output_config, in_place=True)

    with CaptureStdout() as graph_dot:
        snakemake.snakemake(
            snakefile=get_snakefile(analysis_type, sequencing_type),
            dryrun=True,
            configfiles=[output_config],
            printrulegraph=True,
        )

    graph_title = "_".join(["BALSAMIC", bv, json_out["analysis"]["case_id"]])
    graph_dot = "".join(graph_dot).replace(
        "snakemake_dag {",
        'BALSAMIC { label="' + graph_title + '";labelloc="t";')
    graph_obj = graphviz.Source(graph_dot,
                                filename=dag_image,
                                format="pdf",
                                engine="dot")

    try:
        graph_pdf = graph_obj.render()
        LOG.info(
            f"BALSAMIC Workflow has been configured successfully - {output_config}"
        )
    except Exception:
        LOG.error(f"BALSAMIC dag graph generation failed - {dag_image}")
        raise click.Abort()



class CaseConfigAttributes:
    """Store and return configuration attributes
    """

    def __init__(self):
        self.QC = {}
        self.vcf = {}
        self.analysis = {}
        self.samples = {}
        self.reference = ""
        self.conda_env_yaml = ""
        self.rule_directory = ""
        self.singularity = ""
        self.bioinfo_tools = {}


class CaseConfigAssembler:
    """Calls functions to construct the JSON config"""

    def __init__(self, context):
        self.context = context
        self.config = CaseConfigAttributes()

    @property
    def package_root_path(self):
        return f'{os.path.dirname(sys.modules["BALSAMIC"].__file__)}/'

    @property
    def conda_env_path(self):
        return f"{self.package_root_path}conda/"

    @property
    def analysis_json_recipe_path(self):
        return f"{self.package_root_path}config/analysis.json"

    @property
    def sample_json_recipe_path(self):
        return f"{self.package_root_path}config/sample.json"

    @property
    def config_output_path(self):
        return f'{self.context["case_id"]}_{datetime.now().strftime("%Y%m%d")}.json'

    @property
    def analysis_output_path(self):
        return os.path.abspath(self.context["analysis_dir"])

    @property
    def case_output_path(self):
        return f'{self.analysis_output_path}/{self.context["case_id"]}'

    @property
    def fastq_path(self):
        return f'{self.case_output_path}/analysis/fastq'

    @property
    def dag_image_path(self):
        return f'{self.case_output_path}/{self.config_output_path}_BALSAMIC_{BALSAMIC.__version__}_graph'

    def assemble(self):
        self.get_bioinfo_tools_list()
        self.get_sample_config()
        self.get_analysis_config()
        self.get_sample_names()
        self.get_read_reference_json()
        self.get_balsamic_env_config_path()
        self.get_singularity_image()
        self.config.rule_directory = self.package_root_path



    def get_balsamic_env_config_path(self):
        self.config.conda_env_yaml = f"{self.package_root_path}config/balsamic_env.yaml"

    def get_singularity_image(self):
        self.config.singularity = {
            "image": os.path.abspath(self.context["singularity"])
        }

    def get_bioinfo_tools_list(self):
        bioinfo_tools = {}
        core_packages = [
            "bwa", "bcftools", "cutadapt", "fastqc", "gatk", "manta", "picard",
            "sambamba", "strelka", "samtools", "tabix", "vardic"
        ]

        for yaml_file in os.listdir(self.conda_env_path):
            if yaml_file.endswith(".yaml"):
                with open(self.conda_env_path + yaml_file, "r") as f:
                    packages = yaml.safe_load(f)["dependencies"]
                    for p in packages:
                        try:
                            name, version = p.split("=")
                        except ValueError:
                            name, version = p, ""
                        finally:
                            if name in core_packages:
                                bioinfo_tools[name] = version

        self.config.bioinfo_tools = bioinfo_tools


    def get_sample_config(self):
        with open(self.sample_json_recipe_path) as sample_json:
            sample_config = json.load(sample_json)
            self.config.analysis = sample_config["analysis"]

        self.config.analysis["case_id"] = self.context["case_id"]
        self.config.analysis["config_creation_date"] = datetime.now().strftime("%Y-%m-%d %H:%M")
        self.config.analysis["analysis_dir"] = f"{self.analysis_output_path}/"
        self.config.analysis["log"] = f"{self.case_output_path}/logs"
        self.config.analysis["script"] = f"{self.case_output_path}/scripts"
        self.config.analysis["result"] = f"{self.case_output_path}/analysis"
        self.config.analysis["benchmark"] = f"{self.case_output_path}/benchmarks"
        self.config.analysis["fastq_path"] = f"{self.fastq_path}/"
        self.config.analysis["BALSAMIC_version"] = BALSAMIC.__version__
        self.config.analysis["dag"] = f"{self.dag_image_path}.pdf"
        self.get_analysis_type()
        self.get_sequencing_type()


    def get_analysis_config(self):
        with open(self.analysis_json_recipe_path) as analysis_json:
            analysis_config = json.load(analysis_json)
            self.config.QC = analysis_config["QC"]
            self.config.vcf = analysis_config["vcf"]

        if self.context["umi"] and self.context["umi_trim_length"] > 0:
            self.config.QC["umi_trim"] = self.context["umi"]
            self.config.QC["umi_trim_length"] = str(self.context["umi_trim_length"])
        self.config.QC["quality_trim"] = self.context["quality_trim"]
        self.config.QC["adapter_trim"] = self.context["adapter_trim"]


    def get_analysis_type(self):
        if self.context["normal"]:
            self.config.analysis["analysis_type"] = "single"
        else:
            self.config.analysis["analysis_type"] = "paired"

    def get_sequencing_type(self):
        if self.context["panel_bed"]:
            self.config.analysis["sequencing_type"] = "targeted"
        else:
            self.config.analysis["sequencing_type"] = "wgs"

    def get_read_reference_json(self):
        with open(self.context["reference_config"]) as fh:
            ref_json = json.load(fh)
            for k, v in ref_json['reference'].items():
                ref_json['reference'][k] = os.path.abspath(v)

            self.config.reference = ref_json

    def create_output_dir(self):
        os.makedirs(self.case_output_path, exist_ok=True)

    def create_fastq_symlink_dir(self):
        os.makedirs(self.fastq_path, exist_ok=True)

    def validate_existing_fastq(self):
        for sample in [self.context["tumor"] + self.context["normal"]]:
            if os.path.exists(sample):
                continue
            else:
                LOG.error(f"{sample} is not found, update correct file path")
                return False
        return True

    def validate_fastq_pattern(self, sample):
        fq_pattern = re.compile(r"R_[12]" + self.context["fastq_prefix"] +
                                ".fastq.gz$")
        sample_basename = os.path.basename(sample)
        try:
            file_str = sample_basename[0:(
                fq_pattern.search(sample_basename).span()[0] + 1)]
            return file_str

        except AttributeError:
            LOG.error(
                f"File name is invalid, fastq file should be sample_R_1.fastq.gz"
            )
            raise click.Abort()

    def get_sample_names(self):
        for sample in self.context["tumor"]:
            file_str = self.validate_fastq_pattern(sample)
            self.config.samples[file_str] = {
                "file_prefix": file_str,
                "type": "tumor",
                "readpair_suffix": ["1", "2"]
            }

        for sample in self.context["normal"]:
            file_str = self.validate_fastq_pattern(sample)
            self.config.samples[file_str] = {
                "file_prefix": file_str,
                "type": "normal",
                "readpair_suffix": ["1", "2"]
            }

    def copy_fastq_symlink_path(self):
        for sample in self.context["tumor"] + self.context["normal"]:
            sample_basename = os.path.basename(sample)
            sample_abspath = os.path.abspath(sample)
            sample_dest = f'{self.fastq_path}/{sample_basename}'
            if sample_abspath == sample_dest:
                LOG.info("File reference already a copy")
                continue
            else:
                try:
                    shutil.copyfile(sample_abspath, sample_dest)
                except shutil.SameFileError as e:
                    LOG.warning(e)
                    LOG.warning(
                        f"Desitination file {sample_dest} exists. No copy link was created."
                    )
        

    def return_dict(self):
        return self.config.__dict__

    def save_config_json(self):
        """Saves json, and returns the path"""
        pass

    def reformat_json(self):
        """Restructures saved json"""
        pass


@click.command("case",
               short_help="Create a sample config file from input sample data")
@click.option(
    "--umi/--no-umi",
    default=True,
    show_default=True,
    is_flag=True,
    help="UMI processing steps for samples with UMI tags",
)
@click.option(
    "--umi-trim-length",
    default=5,
    show_default=True,
    type=int,
    help="Trim N bases from reads in fastq",
)
@click.option(
    "--quality-trim/--no-quality-trim",
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim low quality reads in fastq",
)
@click.option(
    "--adapter-trim/--no-adapter-trim",
    default=False,
    show_default=True,
    is_flag=True,
    help="Trim adapters from reads in fastq",
)
@click.option(
    "-r",
    "--reference-config",
    required=True,
    show_default=True,
    type=click.Path(),
    help="Reference config file.",
)
@click.option("-p",
              "--panel-bed",
              type=click.Path(),
              help="Panel bed file for variant calling.")
@click.option(
    "-t",
    "--tumor",
    required=True,
    multiple=True,
    help="Fastq files for tumor sample. \
              Example: if files are tumor_fqreads_1.fastq.gz tumor_fqreads_2.fastq.gz, \
              the input should be --tumor tumor_fqreads",
)
@click.option(
    "-n",
    "--normal",
    required=False,
    multiple=True,
    help="Fastq files for normal sample. \
              Example: if files are normal_fqreads_1.fastq.gz normal_fqreads_2.fastq.gz, \
              the input should be --normal normal_fqreads",
)
@click.option(
    "--case-id",
    required=True,
    help="Sample id that is used for reporting, \
              naming the analysis jobs, and analysis path",
)
@click.option(
    "--analysis-dir",
    type=click.Path(),
    default=".",
    help="Root analysis path to store \
              analysis logs and results. The final path will be analysis-dir/sample-id",
)
@click.option(
    "--singularity",
    type=click.Path(),
    required=True,
    help="Download singularity image for BALSAMIC",
)
@click.option("--fastq-prefix",
              required=False,
              default="",
              help="Prefix to fastq file. \
              The string that comes after readprefix")
@click.pass_context
def case_config(context, **args):
    assembler = CaseConfigAssembler(context.params)
    assembler.assemble()
    assembler.copy_fastq_symlink_path()
    print(assembler.config.build_config())



if __name__ == "__main__":
    case_config()
