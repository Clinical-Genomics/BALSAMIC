"""Balsamic config case CLI."""
import json
import logging
import os
from pathlib import Path

import click

from BALSAMIC.commands.config.options import (
    OPTION_CASE_ID,
    OPTION_GENDER,
    OPTION_UMI,
    OPTION_UMI_TRIM_LENGTH,
    OPTION_QUALITY_TRIM,
    OPTION_ADAPTER_TRIM,
    OPTION_FASTQ_PATH,
    OPTION_PANEL_BED,
    OPTION_BACKGROUND_VARIANTS,
    OPTION_PON_CNN,
    OPTION_BALSAMIC_CACHE,
    OPTION_CONTAINER_VERSION,
    OPTION_ANALYSIS_DIR,
    OPTION_TUMOR_SAMPLE_NAME,
    OPTION_NORMAL_SAMPLE_NAME,
    OPTION_CADD_ANNOTATIONS,
    OPTION_CLINICAL_SNV_OBSERVATIONS,
    OPTION_CLINICAL_SV_OBSERVATIONS,
    OPTION_CANCER_GERMLINE_SNV_OBSERVATIONS,
    OPTION_CANCER_SOMATIC_SNV_OBSERVATIONS,
    OPTION_CANCER_SOMATIC_SV_OBSERVATIONS,
    OPTION_SWEGEN_SNV,
    OPTION_SWEGEN_SV,
    OPTION_ANALYSIS_WORKFLOW,
)
from BALSAMIC.commands.options import OPTION_GENOME_VERSION
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV, Gender, AnalysisWorkflow
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.paths import CONTAINERS_DIR
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.models.analysis import ConfigModel
from BALSAMIC.utils.cli import (
    get_sample_list,
    get_panel_chrom,
    get_bioinfo_tools_version,
    generate_graph,
    get_analysis_fastq_files_directory,
)
from BALSAMIC.utils.io import write_json

LOG = logging.getLogger(__name__)


@click.command("case", short_help="Create a sample config file from input sample data")
@OPTION_ADAPTER_TRIM
@OPTION_ANALYSIS_DIR
@OPTION_ANALYSIS_WORKFLOW
@OPTION_BACKGROUND_VARIANTS
@OPTION_BALSAMIC_CACHE
@OPTION_CADD_ANNOTATIONS
@OPTION_CANCER_GERMLINE_SNV_OBSERVATIONS
@OPTION_CANCER_SOMATIC_SNV_OBSERVATIONS
@OPTION_CANCER_SOMATIC_SV_OBSERVATIONS
@OPTION_CASE_ID
@OPTION_CLINICAL_SNV_OBSERVATIONS
@OPTION_CLINICAL_SV_OBSERVATIONS
@OPTION_CONTAINER_VERSION
@OPTION_FASTQ_PATH
@OPTION_GENDER
@OPTION_GENOME_VERSION
@OPTION_NORMAL_SAMPLE_NAME
@OPTION_PANEL_BED
@OPTION_PON_CNN
@OPTION_QUALITY_TRIM
@OPTION_SWEGEN_SNV
@OPTION_SWEGEN_SV
@OPTION_TUMOR_SAMPLE_NAME
@OPTION_UMI
@OPTION_UMI_TRIM_LENGTH
@click.pass_context
def case_config(
    context: click.Context,
    adapter_trim: bool,
    analysis_dir: Path,
    analysis_workflow: AnalysisWorkflow,
    background_variants: Path,
    balsamic_cache: Path,
    cadd_annotations: Path,
    cancer_germline_snv_observations: Path,
    cancer_somatic_snv_observations: Path,
    cancer_somatic_sv_observations: Path,
    case_id: str,
    clinical_snv_observations: Path,
    clinical_sv_observations: Path,
    container_version: str,
    fastq_path: Path,
    gender: Gender,
    genome_version: GenomeVersion,
    normal_sample_name: str,
    panel_bed: Path,
    pon_cnn: Path,
    quality_trim: bool,
    swegen_snv: Path,
    swegen_sv: Path,
    tumor_sample_name: str,
    umi: bool,
    umi_trim_length: int,
):
    if container_version:
        balsamic_version = container_version

    reference_config = os.path.join(
        balsamic_cache, balsamic_version, genome_version, "reference.json"
    )
    with open(reference_config, "r") as config_file:
        reference_dict = json.load(config_file)

    cadd_annotations_path = {
        "cadd_annotations": cadd_annotations,
    }
    if cadd_annotations:
        reference_dict.update(cadd_annotations_path)

    variants_observations = {
        "clinical_snv_observations": clinical_snv_observations,
        "clinical_sv_observations": clinical_sv_observations,
        "cancer_germline_snv_observations": cancer_germline_snv_observations,
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

    analysis_fastq_dir: str = get_analysis_fastq_files_directory(
        case_dir=Path(analysis_dir, case_id).as_posix(), fastq_path=fastq_path
    )

    config_collection_dict = ConfigModel(
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
            "fastq_path": analysis_fastq_dir,
            "analysis_type": "paired" if normal_sample_name else "single",
            "sequencing_type": "targeted" if panel_bed else "wgs",
            "analysis_workflow": analysis_workflow,
        },
        reference=reference_dict,
        singularity={
            "image": os.path.join(balsamic_cache, balsamic_version, "containers")
        },
        background_variants=background_variants,
        samples=get_sample_list(
            tumor_sample_name=tumor_sample_name,
            normal_sample_name=normal_sample_name,
            fastq_path=analysis_fastq_dir,
        ),
        vcf=VCF_DICT,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            bioinfo_tools=BIOINFO_TOOL_ENV,
            container_conda_env_path=CONTAINERS_DIR,
        ),
        panel={
            "capture_kit": panel_bed,
            "chrom": get_panel_chrom(panel_bed),
            "pon_cnn": pon_cnn,
        }
        if panel_bed
        else None,
    ).dict(by_alias=True, exclude_none=True)
    LOG.info("Balsamic config model instantiated successfully")

    result_path: Path = Path(config_collection_dict["analysis"]["result"])
    log_path: Path = Path(config_collection_dict["analysis"]["log"])
    script_path: Path = Path(config_collection_dict["analysis"]["script"])
    benchmark_path: Path = Path(config_collection_dict["analysis"]["benchmark"])

    # Create directories for results, logs, scripts and benchmark files
    analysis_directories_list = [result_path, log_path, script_path, benchmark_path]

    for analysis_sub_dir in analysis_directories_list:
        analysis_sub_dir.mkdir(exist_ok=True)

    config_path = Path(analysis_dir, case_id, case_id + ".json").as_posix()
    write_json(json_obj=config_collection_dict, path=config_path)
    LOG.info(f"Config file saved successfully - {config_path}")

    generate_graph(config_collection_dict, config_path)
    LOG.info(f"BALSAMIC Workflow has been configured successfully!")
