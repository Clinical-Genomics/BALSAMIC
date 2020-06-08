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
        self.assemble_config()


    @property
    def package_root_path(self):
        return f'{os.path.dirname(sys.modules["BALSAMIC"].__file__)}'

    @property
    def conda_env_path(self):
        return f"{self.package_root_path}/conda/"

    @property
    def analysis_json_recipe_path(self):
        return f"{self.package_root_path}/config/analysis.json"

    @property
    def sample_json_recipe_path(self):
        return f"{self.package_root_path}/config/sample.json"

    @property
    def config_output_path(self):
        return f'{self.case_output_path}/{self.context["case_id"]}_{datetime.now().strftime("%Y%m%d")}.json'

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
        return f'{self.config_output_path}_BALSAMIC_{BALSAMIC.__version__}_graph'


    def get_balsamic_env_config_path(self):
        self.config.conda_env_yaml = f"{self.package_root_path}/config/balsamic_env.yaml"

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
        self.config.analysis["config_creation_date"] = datetime.now().strftime(
            "%Y-%m-%d %H:%M")
        self.config.analysis["analysis_dir"] = f"{self.analysis_output_path}/"
        self.config.analysis["log"] = f"{self.case_output_path}/logs"
        self.config.analysis["script"] = f"{self.case_output_path}/scripts"
        self.config.analysis["result"] = f"{self.case_output_path}/analysis"
        self.config.analysis[
            "benchmark"] = f"{self.case_output_path}/benchmarks"
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
            self.config.QC["umi_trim_length"] = str(
                self.context["umi_trim_length"])
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


    def verify_panel_path(self):
        if os.path.exists(self.context["panel_bed"]):
            return True
        else:
            click.Abort()
            return False

    def get_panel_path(self):
        self.config.panel["capture_kit"] = os.path.abspath(
            self.context["panel_bed"])

    def get_panel_chrom(self):
        lines = [
            line.rstrip('\n') for line in open(self.context["panel_bed"], 'r')
        ]
        chrom = list(set([s.split('\t')[0] for s in lines]))
        self.config.panel["chrom"] = chrom

    def get_panel_config(self):
        if self.context["panel_bed"]:
            if self.verify_panel_path():
                self.config.panel = {}
                self.get_panel_path()
                self.get_panel_chrom()

    def assemble_config(self):
        self.get_bioinfo_tools_list()
        self.get_sample_config()
        self.get_analysis_config()
        self.get_sample_names()
        self.get_read_reference_json()
        self.get_balsamic_env_config_path()
        self.get_singularity_image()
        self.config.rule_directory = self.package_root_path
        self.get_panel_config()

        LOG.info("Creating directories")
        self.create_output_dir()
        self.create_fastq_symlink_dir()


    def copy_fastq_symlink_path(self):
        for sample in self.context["tumor"] + self.context["normal"]:
            sample_basename = os.path.basename(sample)
            sample_abspath = os.path.abspath(sample)
            sample_dest = f'{self.fastq_path}/{sample_basename}'
            if sample_abspath == sample_dest:
                continue
            else:
                try:
                    shutil.copyfile(sample_abspath, sample_dest)
                except shutil.SameFileError as e:
                    LOG.warning(e)
                    LOG.warning(
                        f"Desitination file {sample_dest} exists. No copy link was created."
                    )

    def build_config(self):
        return self.config.__dict__

    def save_config_json(self):
        with open(self.config_output_path, "w") as fh:
            json.dump(self.build_config(), fh, indent=4)

    def reformat_config_json(self):
        FormatFile(self.config_output_path)

    def build_graph(self):
        with CaptureStdout() as graph_dot:
            snakemake.snakemake(snakefile=get_snakefile(
                self.config.analysis["analysis_type"],
                self.config.analysis["sequencing_type"]),
                                dryrun=True,
                                configfiles=[self.config_output_path],
                                printrulegraph=True)

        graph_title = "_".join([
            'BALSAMIC', BALSAMIC.__version__, self.config.analysis["case_id"]
        ])
        graph_dot = "".join(graph_dot).replace(
            'snakemake_dag {',
            'BALSAMIC { label="' + graph_title + '";labelloc="t";')
        graph_obj = graphviz.Source(graph_dot,
                                    filename=self.dag_image_path,
                                    format="pdf",
                                    engine="dot")

        try:
            graph_pdf = graph_obj.render()
            LOG.info(
                f'BALSAMIC Workflow has been configured successfully - {self.config_output_path}'
            )
        except Exception:
            LOG.error(
                f'BALSAMIC dag graph generation failed - {self.dag_image_path}'
            )
            raise click.Abort()


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
    assembler.copy_fastq_symlink_path()
    LOG.info("Saving")
    assembler.save_config_json()
    assembler.reformat_config_json()
    assembler.build_graph()
