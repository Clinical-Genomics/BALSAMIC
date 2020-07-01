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
from BALSAMIC.utils.cli import CaptureStdout, get_snakefile
from BALSAMIC.config.constants import BalsamicConfigModel, CONDA_ENV_PATH

LOG = logging.getLogger(__name__)


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
            packages = yaml.safe_load(f)["dependencies"]
            for p in packages:
                try:
                    name, version = p.split("=")
                except ValueError:
                    name, version = p, None
                finally:
                    bioinfo_tools[name] = version
    return bioinfo_tools


def get_sample_dict(tumor, normal):
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


def create_fastq_symlink(filename, symlink_dir: Path):
    parent_dir = Path(filename).parents[0]
    file_str = validate_fastq_pattern(filename)

    for f in parent_dir.rglob(f'*{file_str}*.fastq.gz'):
        try:
            Path(symlink_dir, f.name).symlink_to(f)
        except FileExistsError:
            LOG.info(f"Path {symlink_dir / f.name} exists, skipping")


def create_working_directories(config_collection_dict):
    Path.mkdir(Path(config_collection_dict["analysis"]["fastq_path"]),
               parents=True,
               exist_ok=True)


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


@click.command("case",
               short_help="Create a sample config file from input sample data")
@click.option("--case-id",
              required=True,
              help="Sample id that is used for reporting, \
    naming the analysis jobs, and analysis path")
@click.option("--umi/--no-umi",
              default=True,
              show_default=True,
              is_flag=True,
              help="UMI processing steps for samples with UMI tags")
@click.option("--umi-trim-length",
              default=5,
              show_default=True,
              type=int,
              help="Trim N bases from reads in fastq")
@click.option("--quality-trim/--no-quality-trim",
              default=True,
              show_default=True,
              is_flag=True,
              help="Trim low quality reads in fastq")
@click.option("--adapter-trim/--no-adapter-trim",
              default=False,
              show_default=True,
              is_flag=True,
              help="Trim adapters from reads in fastq")
@click.option("-r",
              "--reference-config",
              required=True,
              type=click.Path(exists=True, resolve_path=True),
              help="Reference config file.")
@click.option("-p",
              "--panel-bed",
              type=click.Path(exists=True, resolve_path=True),
              required=False,
              help="Panel bed file for variant calling.")
@click.option("--singularity",
              type=click.Path(exists=True, resolve_path=True),
              required=True,
              help="Download singularity image for BALSAMIC")
@click.option("--analysis-dir",
              type=click.Path(exists=True, resolve_path=True),
              default=".",
              help="Root analysis path to store analysis logs and results. \
                                     The final path will be analysis-dir/sample-id"
              )
@click.option("-t",
              "--tumor",
              type=click.Path(exists=True, resolve_path=True),
              required=True,
              multiple=True,
              help="Fastq files for tumor sample.")
@click.option("-n",
              "--normal",
              type=click.Path(exists=True, resolve_path=True),
              required=False,
              multiple=True,
              help="Fastq files for normal sample.")
@click.option(
    "-f",
    "--format",
    required=False,
    type=click.Choice(["fastq", "vcf", "bam"]),
    default="fastq",
    show_default=True,
)
@click.pass_context
def case_config(context, case_id, umi, umi_trim_length, adapter_trim,
                quality_trim, reference_config, panel_bed, singularity,
                analysis_dir, tumor, normal, format):

    try:
        samples = get_sample_dict(tumor, normal)
    except AttributeError:
        LOG.error(
            f"File name is invalid, use convention [SAMPLE_ID]_R_[1,2].fastq.gz")
        click.Abort()

    try:
        reference_dict = json.load(open(reference_config))["reference"]
    except Exception as e:
        LOG.error(f"Reference config {reference_config} does not follow correct format: {e}")
        click.Abort()

    try:
        bioinfo_tools = get_bioinfo_tools_list(CONDA_ENV_PATH)
    except Exception as e:
        LOG.error(f"Error generating a list of bioinfo tools: {e}")
        click.Abort()


    try:
        config_collection = BalsamicConfigModel(
            QC={
                "quality_trim": quality_trim,
                "adapter_trim": adapter_trim,
                "umi_trim": umi,
                "umi_trim_length": umi_trim_length,
            },
            analysis={
                "case_id": case_id,
                "analysis_dir": analysis_dir,
                "analysis_type": "paired" if normal else "single",
                "sequencing_type": "targeted" if panel_bed else "wgs",
            },
            panel={
                "capture_kit": panel_bed,
                "chrom": get_panel_chrom(panel_bed),
            } if panel_bed else None,
            bioinfo_tools=bioinfo_tools,
            reference=reference_dict,
            singularity=singularity,
            samples=samples,
            vcf={},
        )
    except Exception as e:
        LOG.error(f"Failed to generate config file: {e}")
        click.Abort()

    config_collection_dict = config_collection.dict(by_alias=True)

    try:
        create_working_directories(config_collection_dict)
    except Exception as e:
        LOG.error(f"Could not create directories: {e}")
        click.Abort()

    for filename in tumor + normal:
        create_fastq_symlink(
            filename,
            Path(config_collection_dict["analysis"]["fastq_path"]))

    config_path = Path(analysis_dir) / case_id / (case_id + ".json")
    with open(config_path, "w+") as fh:
        fh.write(json.dumps(config_collection_dict, indent=4))

    try:
        generate_graph(config_collection_dict, config_path)
        LOG.info(
            f'BALSAMIC Workflow has been configured successfully - {config_path}')
    except ValueError as e:
        LOG.error(
            f'BALSAMIC dag graph generation failed - {config_collection_dict["analysis"]["dag"]}',
        )
        click.Abort()
            

