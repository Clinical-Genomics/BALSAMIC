import os
import re
import json
import shutil
import logging
import click
import snakemake
import graphviz
import sys
import yaml
import BALSAMIC

from pathlib import Path
from datetime import datetime
from yapf.yapflib.yapf_api import FormatFile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.errors import ExtensionNotSupportedError
from BALSAMIC.config.constats import QC, VCF, ANALYSIS, BIOINFO_BASE, CONDA_ENV_PATH, CONDA_ENV_YAML, SUPPORTED_EXTENSIONS, RULE_DIRECTORY

LOG = logging.getLogger(__name__)

CASE_ID_OPTION = click.Option("--case-id",
                              required=True,
                              help="Sample id that is used for reporting, \
                                    naming the analysis jobs, and analysis path"
                              )
UMI_OPTION = click.Option(
    "--umi/--no-umi",
    default=True,
    show_default=True,
    is_flag=True,
    help="UMI processing steps for samples with UMI tags")
UMI_TRIM_LENGTH_OPTION = click.Option("--umi-trim-length",
                                      default=5,
                                      show_default=True,
                                      type=int,
                                      help="Trim N bases from reads in fastq")
QUALITY_TRIM_OPTION = click.Option("--quality-trim/--no-quality-trim",
                                   default=True,
                                   show_default=True,
                                   is_flag=True,
                                   help="Trim low quality reads in fastq")
ADAPTER_TRIM_OPTION = click.Option("--adapter-trim/--no-adapter-trim",
                                   default=False,
                                   show_default=True,
                                   is_flag=True,
                                   help="Trim adapters from reads in fastq")
REFERENCE_CONFIG_OPTION = click.Option("-r",
                                       "--reference-config",
                                       required=True,
                                       type=click.Path(),
                                       exists=True,
                                       resolve_path=True,
                                       help="Reference config file.")
PANEL_BED_OPTION = click.Option("-p",
                                "--panel-bed",
                                type=click.Path(),
                                required=False,
                                exists=True,
                                resolve_path=True,
                                help="Panel bed file for variant calling.")
SINGULARITY_OPTION = click.Option(
    "--singularity",
    type=click.Path(),
    required=True,
    resolve_path=True,
    exists=True,
    help="Download singularity image for BALSAMIC")
FASTQ_PREFIX_OPTION = click.Option(
    "--fastq-prefix",
    required=False,
    default="",
    help="Prefix to fastq file. The string that comes after readprefix")
ANALYSIS_DIR_OPTION = click.Option(
    "--analysis-dir",
    type=click.Path(),
    default=".",
    resolve_path=True,
    exists=True,
    help="Root analysis path to store analysis logs and results. \
                                     The final path will be analysis-dir/sample-id"
)
TUMOR_OPTION = click.Option("-t",
                            "--tumor",
                            required=True,
                            multiple=True,
                            exists=True,
                            resolve_path=True,
                            help="Fastq files for tumor sample.")
NORMAL_OPTION = click.Option("-n",
                             "--normal",
                             required=False,
                             multiple=True,
                             exists=True,
                             resolve_path=True,
                             help="Fastq files for normal sample.")

FILE_FORMAT_OPTION = click.Option("-f", "--format", required=False, type=click.Choice(["fastq", "vcf", "bam"], default="fastq"))




class InputFile:

    """Class for handling input file operations"""
    def __init__(self, path, sampletype, input_format):
        self.path = Path(path)
        self.sampletype = sampletype
        self.input_format = input_format
        self.prefix = str()


    def verify_exists(self, dst) -> bool:
        pass
        

    def create_copy(self, dst) -> None:
        if self.verify_exists(dst):
            continue
        else:
            shutil.copy()

    def extract_prefix(self) -> None:
        pass

 

class CaseConfigIOHandler:
    """Handle arguments relevant to disc use:
        ARGS: case_id, analysis_dir, 

        Check input file format
        Check input files naming convention
        Create directories
        Create copy of input files in working directory

    """

    def __init__(self, case_id, analysis_dir, files, input_format):
        self.case_id = case_id
        self.analysis_dir = analysis_dir
        self.input_files = files

        self.case_output_dir = Path(self.analysis_dir / self.case_id )
        self.working_copy_dir = Path(self.analysis_dir / self.case_id / "analysis" / input_format)


        def create_output_dirs(self):
            Path.mkdir(self.working_copy_dir, parents=True, exist_ok=True)

        def create_working_file_copy(self):
            for sample in self.files:
                sample.create_copy(self.working_copy_dir)


class CaseConfigAttributes:
    def __init__(self):
        self.QC = QC
        self.vcf = VCF
        self.analysis = ANALYSIS
        self.conda_env_yaml = CONDA_ENV_YAML
        self.rule_directory = RULE_DIRECTORY
        self.reference = {}
        self.samples = {}
        self.panel_bed = {}
        self.singularity = {}

    def to_dict(self) -> dict:
        return self.__dict__

    def to_json(self, save_path) -> None:
        pass




class CaseConfigAssembler:
    """Calls functions to construct the JSON config"""

    def __init__(self, ):
        self.config_attributes = CaseConfigAttributes()



    @property
    def config_output_path(self) -> Path:
        return f'{self.case_output_path}/{self.context["case_id"]}_{datetime.now().strftime("%Y%m%d")}.json'


    @property
    def dag_image_path(self):
        return f'{self.config_output_path}_BALSAMIC_{BALSAMIC.__version__}_graph'

    def get_singularity_image(self):
        self.config.singularity = {
            "image": os.path.abspath(self.context["singularity"])
        }

    def get_bioinfo_tools_list(self):
        bioinfo_tools = {}

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

        self.config.analysis[
            "benchmark"] = f"{self.case_output_path}/benchmarks"
        self.config.analysis["fastq_path"] = f"{self.fastq_path}/"

        self.config.analysis["dag"] = f"{self.dag_image_path}.pdf"
        self.get_analysis_type()
        self.get_sequencing_type()



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

    def get_panel_path(self):
        self.config.panel["capture_kit"] = os.path.abspath(
            self.context["panel_bed"])

    def get_panel_chrom(self) -> list:
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

                    )

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
@click.pass_context
def case_config(context, **args):

    assembler = CaseConfigAssembler(context.params)
    assembler.assemble_config()
    assembler.copy_fastq_symlink_path()
    LOG.info("Saving")
    assembler.save_config_json()
    assembler.reformat_config_json()
    assembler.build_graph()
