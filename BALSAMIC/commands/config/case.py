import os
import json
import logging
from pathlib import Path

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.cli import (
    get_sample_dict,
    get_panel_chrom,
    get_bioinfo_tools_version,
    create_fastq_symlink,
    generate_graph,
)
from BALSAMIC.constants.common import (
    CONTAINERS_CONDA_ENV_PATH,
    BIOINFO_TOOL_ENV,
    GENDER_OPTIONS,
)
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.utils.models import BalsamicConfigModel

LOG = logging.getLogger(__name__)


@click.command("case", short_help="Create a sample config file from input sample data")
@click.option(
    "--case-id",
    required=True,
    help="Sample id that is used for reporting, \
              naming the analysis jobs, and analysis path",
)
@click.option(
    "--gender",
    required=False,
    default="female",
    show_default=True,
    type=click.Choice(GENDER_OPTIONS),
    help="Case associated gender",
)
@click.option(
    "--umi/--no-umi",
    default=True,
    show_default=True,
    is_flag=True,
    help=(
        "UMI processing steps for samples with UMI tags."
        "For WGS cases, UMI is always disabled."
    ),
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
    default=True,
    show_default=True,
    is_flag=True,
    help="Trim adapters from reads in fastq",
)
@click.option(
    "-p",
    "--panel-bed",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel bed file for variant calling.",
)
@click.option(
    "-b",
    "--background-variants",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Background set of valid variants for UMI",
)
@click.option(
    "--pon-cnn",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="Panel of normal reference (.cnn) for cnvkit",
)
@click.option(
    "--balsamic-cache",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Path to BALSAMIC cache",
)
@click.option(
    "--container-version",
    show_default=True,
    default=balsamic_version,
    type=click.Choice(["develop", "master", balsamic_version]),
    help="Container for BALSAMIC version to download",
)
@click.option(
    "--analysis-dir",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    help="Root analysis path to store analysis logs and results. \
                                     The final path will be analysis-dir/sample-id",
)
@click.option(
    "-t",
    "--tumor",
    type=click.Path(exists=True, resolve_path=True),
    required=True,
    multiple=True,
    help="Fastq files for tumor sample.",
)
@click.option(
    "-n",
    "--normal",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    multiple=True,
    help="Fastq files for normal sample.",
)
@click.option("--tumor-sample-name", help="Tumor sample name")
@click.option("--normal-sample-name", help="Normal sample name")
@click.option(
    "--clinical-snv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of clinical SNV observations (WGS analysis workflow)",
)
@click.option(
    "--clinical-sv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of clinical SV observations (WGS analysis workflow)",
)
@click.option(
    "--cancer-all-snv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of cancer SNV normal observations (WGS analysis workflow)",
)
@click.option(
    "--cancer-somatic-snv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of cancer SNV tumor observations (WGS analysis workflow)",
)
@click.option(
    "--cancer-somatic-sv-observations",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of cancer SV observations (WGS analysis workflow)",
)
@click.option(
    "--swegen-snv",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of Swegen SNV frequency database (WGS analysis workflow)",
)
@click.option(
    "--swegen-sv",
    type=click.Path(exists=True, resolve_path=True),
    required=False,
    help="VCF path of Swegen SV frequency database (WGS analysis workflow)",
)
@click.option(
    "-g",
    "--genome-version",
    default="hg19",
    type=click.Choice(["hg19", "hg38", "canfam3"]),
    help=(
        "Genome version to prepare reference. Path to genome"
        "will be <outdir>/genome_version"
    ),
)
@click.option(
    "-w",
    "--analysis-workflow",
    default="balsamic",
    show_default=True,
    type=click.Choice(["balsamic", "balsamic-umi", "balsamic-qc"]),
    help=(
        'Analysis workflow to run. By default: "balsamic" only '
        "workflow will be running. If you want to run both "
        "balsamic and UMI workflow together for panel data; "
        'choose "balsamic-umi" option '
    ),
)
@click.pass_context
def case_config(
    context,
    case_id,
    gender,
    umi,
    umi_trim_length,
    adapter_trim,
    quality_trim,
    panel_bed,
    background_variants,
    pon_cnn,
    analysis_dir,
    tumor,
    normal,
    tumor_sample_name,
    normal_sample_name,
    clinical_snv_observations,
    clinical_sv_observations,
    cancer_all_snv_observations,
    cancer_somatic_snv_observations,
    cancer_somatic_sv_observations,
    swegen_snv,
    swegen_sv,
    genome_version,
    balsamic_cache,
    container_version,
    analysis_workflow,
):

    try:
        samples = get_sample_dict(
            tumor=tumor,
            normal=normal,
            tumor_sample_name=tumor_sample_name,
            normal_sample_name=normal_sample_name,
        )
    except AttributeError:
        LOG.error(f"File name is invalid, use convention [SAMPLE_ID]_R_[1,2].fastq.gz")
        raise click.Abort()

    if container_version:
        balsamic_version = container_version

    reference_config = os.path.join(
        balsamic_cache, balsamic_version, genome_version, "reference.json"
    )
    with open(reference_config, "r") as f:
        reference_dict = json.load(f)["reference"]

    variants_observations = {
        "clinical_snv_observations": clinical_snv_observations,
        "clinical_sv_observations": clinical_sv_observations,
        "cancer_all_snv_observations": cancer_all_snv_observations,
        "cancer_somatic_snv_observations": cancer_somatic_snv_observations,
        "cancer_somatic_sv_observations": cancer_somatic_sv_observations,
        "swegen_snv_frequency": swegen_snv,
        "swegen_sv_frequency": swegen_sv,
    }
    reference_dict.update(
        {
            observations: path
            for observations, path in variants_observations.items()
            if path is not None
        }
    )

    config_collection_dict = BalsamicConfigModel(
        QC={
            "quality_trim": quality_trim,
            "adapter_trim": adapter_trim,
            "umi_trim": umi if panel_bed else False,
            "umi_trim_length": umi_trim_length,
        },
        analysis={
            "case_id": case_id,
            "gender": gender,
            "analysis_dir": analysis_dir,
            "analysis_type": "paired" if normal else "single",
            "sequencing_type": "targeted" if panel_bed else "wgs",
            "analysis_workflow": analysis_workflow,
        },
        reference=reference_dict,
        singularity=os.path.join(balsamic_cache, balsamic_version, "containers"),
        background_variants=background_variants,
        samples=samples,
        vcf=VCF_DICT,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            BIOINFO_TOOL_ENV, CONTAINERS_CONDA_ENV_PATH
        ),
        panel={
            "capture_kit": panel_bed,
            "chrom": get_panel_chrom(panel_bed),
            "pon_cnn": pon_cnn,
        }
        if panel_bed
        else None,
    ).dict(by_alias=True, exclude_none=True)
    LOG.info("Config file generated successfully")

    Path.mkdir(
        Path(config_collection_dict["analysis"]["fastq_path"]),
        parents=True,
        exist_ok=True,
    )
    LOG.info("Directories created successfully")

    create_fastq_symlink(
        casefiles=(tumor + normal),
        symlink_dir=Path(config_collection_dict["analysis"]["fastq_path"]),
    )
    LOG.info(f"Symlinks generated successfully")

    config_path = Path(analysis_dir) / case_id / (case_id + ".json")
    with open(config_path, "w+") as fh:
        fh.write(json.dumps(config_collection_dict, indent=4))
    LOG.info(f"Config file saved successfully - {config_path}")

    try:
        generate_graph(config_collection_dict, config_path)
        LOG.info(f"BALSAMIC Workflow has been configured successfully!")
    except ValueError as e:
        LOG.error(
            f'BALSAMIC dag graph generation failed - {config_collection_dict["analysis"]["dag"]}',
        )
        raise click.Abort()
