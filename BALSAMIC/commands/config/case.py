"""Balsamic config case CLI."""

import logging
from datetime import datetime
from pathlib import Path
from typing import Dict

import click

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.commands.options import (
    OPTION_ANALYSIS_DIR,
    OPTION_ANALYSIS_WORKFLOW,
    OPTION_ARTEFACT_SNV_OBSERVATIONS,
    OPTION_BACKGROUND_VARIANTS,
    OPTION_BALSAMIC_CACHE,
    OPTION_CACHE_VERSION,
    OPTION_CADD_ANNOTATIONS,
    OPTION_CANCER_GERMLINE_SNV_OBSERVATIONS,
    OPTION_CANCER_SOMATIC_SNV_OBSERVATIONS,
    OPTION_CANCER_SOMATIC_SV_OBSERVATIONS,
    OPTION_CASE_ID,
    OPTION_CLINICAL_SNV_OBSERVATIONS,
    OPTION_CLINICAL_SV_OBSERVATIONS,
    OPTION_SOFT_FILTER_NORMAL,
    OPTION_EXOME,
    OPTION_TUMOR_FASTQ_PATH,
    OPTION_NORMAL_FASTQ_PATH,
    OPTION_GENDER,
    OPTION_GENOME_INTERVAL,
    OPTION_GENOME_VERSION,
    OPTION_GENS_COV_PON,
    OPTION_GNOMAD_AF5,
    OPTION_NORMAL_SAMPLE_NAME,
    OPTION_PANEL_BED,
    OPTION_PON_CNN,
    OPTION_SENTIEON_INSTALL_DIR,
    OPTION_SENTIEON_LICENSE,
    OPTION_SWEGEN_SNV,
    OPTION_SWEGEN_SV,
    OPTION_TUMOR_SAMPLE_NAME,
    OPTION_UMI_MIN_READS,
)
from BALSAMIC.constants.analysis import (
    BIOINFO_TOOL_ENV,
    AnalysisWorkflow,
    Gender,
    LogFile,
)
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.constants import FileType
from BALSAMIC.constants.paths import (
    CONTAINERS_DIR,
    SENTIEON_DNASCOPE_MODEL,
    SENTIEON_TNSCOPE_MODEL,
)
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.utils.cli import (
    generate_graph,
    get_analysis_fastq_files_directory,
    get_bioinfo_tools_version,
    get_panel_chrom,
    get_sample_list,
    get_gens_references,
)
from BALSAMIC.utils.io import read_json, write_json
from BALSAMIC.utils.utils import get_absolute_paths_dict
from BALSAMIC.utils.logging import add_file_logging

LOG = logging.getLogger(__name__)


@click.command("case", short_help="Create a sample config file from input sample data")
@OPTION_ANALYSIS_DIR
@OPTION_ANALYSIS_WORKFLOW
@OPTION_ARTEFACT_SNV_OBSERVATIONS
@OPTION_BACKGROUND_VARIANTS
@OPTION_BALSAMIC_CACHE
@OPTION_CACHE_VERSION
@OPTION_CADD_ANNOTATIONS
@OPTION_CANCER_GERMLINE_SNV_OBSERVATIONS
@OPTION_CANCER_SOMATIC_SNV_OBSERVATIONS
@OPTION_CANCER_SOMATIC_SV_OBSERVATIONS
@OPTION_CASE_ID
@OPTION_CLINICAL_SNV_OBSERVATIONS
@OPTION_CLINICAL_SV_OBSERVATIONS
@OPTION_SOFT_FILTER_NORMAL
@OPTION_EXOME
@OPTION_TUMOR_FASTQ_PATH
@OPTION_NORMAL_FASTQ_PATH
@OPTION_GENDER
@OPTION_GENOME_VERSION
@OPTION_GENOME_INTERVAL
@OPTION_GENS_COV_PON
@OPTION_GNOMAD_AF5
@OPTION_NORMAL_SAMPLE_NAME
@OPTION_PANEL_BED
@OPTION_PON_CNN
@OPTION_SENTIEON_INSTALL_DIR
@OPTION_SENTIEON_LICENSE
@OPTION_SWEGEN_SNV
@OPTION_SWEGEN_SV
@OPTION_TUMOR_SAMPLE_NAME
@OPTION_UMI_MIN_READS
@click.pass_context
def case_config(
    context: click.Context,
    analysis_dir: Path,
    analysis_workflow: AnalysisWorkflow,
    artefact_snv_observations: Path,
    background_variants: Path,
    balsamic_cache: Path,
    cache_version: str,
    cadd_annotations: Path,
    cancer_germline_snv_observations: Path,
    cancer_somatic_snv_observations: Path,
    cancer_somatic_sv_observations: Path,
    case_id: str,
    clinical_snv_observations: Path,
    clinical_sv_observations: Path,
    exome: bool,
    tumor_fastq_path: Path,
    normal_fastq_path: Path,
    gender: Gender,
    genome_version: GenomeVersion,
    genome_interval: Path,
    gens_coverage_pon: Path,
    gnomad_min_af5: Path,
    normal_sample_name: str,
    panel_bed: Path,
    pon_cnn: Path,
    sentieon_install_dir: Path,
    sentieon_license: str,
    soft_filter_normal: bool,
    swegen_snv: Path,
    swegen_sv: Path,
    tumor_sample_name: str,
    umi_min_reads: str | None,
):
    """Configure BALSAMIC workflow based on input arguments."""

    LOG.info(f"Starting configuring analysis case: {case_id}.")

    LOG.info(f"Creating case analysis directory: {analysis_dir}/{case_id}.")
    Path(analysis_dir, case_id).mkdir(exist_ok=True)

    log_file = Path(analysis_dir, case_id, LogFile.LOGNAME).as_posix()
    LOG.info(f"Setting BALSAMIC logfile path to: {log_file}.")
    add_file_logging(log_file, logger_name=__name__)

    LOG.info(f"Running BALSAMIC version {balsamic_version} -- CONFIG CASE")
    LOG.info(f"BALSAMIC started with log level {context.obj['log_level']}.")

    LOG.info("Collecting reference and annotation file paths.")
    references_path: Path = Path(balsamic_cache, cache_version, genome_version)
    references: Dict[str, Path] = get_absolute_paths_dict(
        base_path=references_path,
        data=read_json(Path(references_path, f"reference.{FileType.JSON}").as_posix()),
    )
    cadd_annotations_path = {"cadd_annotations": cadd_annotations}
    if cadd_annotations:
        references.update(cadd_annotations_path)

    if analysis_workflow is not AnalysisWorkflow.BALSAMIC_QC:
        gens_references: dict[str, str] = get_gens_references(
            genome_interval=genome_interval,
            gens_coverage_pon=gens_coverage_pon,
            gnomad_min_af5=gnomad_min_af5,
            panel_bed=panel_bed,
        )
        if gens_references:
            # Update references dictionary with GENS reference-files
            references.update(
                {
                    gens_file: path
                    for gens_file, path in gens_references.items()
                    if path is not None
                }
            )
    variants_observations = {
        "artefact_snv_observations": artefact_snv_observations,
        "clinical_snv_observations": clinical_snv_observations,
        "clinical_sv_observations": clinical_sv_observations,
        "cancer_germline_snv_observations": cancer_germline_snv_observations,
        "cancer_somatic_snv_observations": cancer_somatic_snv_observations,
        "cancer_somatic_sv_observations": cancer_somatic_sv_observations,
        "swegen_snv_frequency": swegen_snv,
        "swegen_sv_frequency": swegen_sv,
    }
    references.update(
        {
            observations: path
            for observations, path in variants_observations.items()
            if path is not None
        }
    )
    LOG.info(f"Collected references: {references}")

    result_dir: Path = Path(analysis_dir, case_id, "analysis")
    log_dir: Path = Path(analysis_dir, case_id, "logs")
    script_dir: Path = Path(analysis_dir, case_id, "scripts")
    benchmark_dir: Path = Path(analysis_dir, case_id, "benchmarks")
    dag_path: Path = Path(
        analysis_dir, case_id, f"{case_id}_BALSAMIC_{balsamic_version}_graph.pdf"
    )
    for directory in [result_dir, log_dir, script_dir, benchmark_dir]:
        directory.mkdir(exist_ok=True)

    LOG.info("Created analysis and log directories.")
    LOG.info("Validating configuration data in pydantic model.")
    config_collection_dict = ConfigModel(
        sentieon={
            "sentieon_install_dir": sentieon_install_dir,
            "sentieon_license": sentieon_license,
            "sentieon_exec": Path(sentieon_install_dir, "bin", "sentieon").as_posix(),
            "dnascope_model": SENTIEON_DNASCOPE_MODEL.as_posix(),
            "tnscope_model": SENTIEON_TNSCOPE_MODEL.as_posix(),
        },
        analysis={
            "case_id": case_id,
            "soft_filter_normal": soft_filter_normal,
            "gender": gender,
            "analysis_dir": analysis_dir,
            "analysis_type": "paired" if normal_sample_name else "single",
            "tumor_fastq_path": tumor_fastq_path,
            "normal_fastq_path": normal_fastq_path,
            "log": log_dir.as_posix(),
            "script": script_dir.as_posix(),
            "result": result_dir.as_posix(),
            "benchmark": benchmark_dir.as_posix(),
            "dag": dag_path.as_posix(),
            "sequencing_type": "targeted" if panel_bed else "wgs",
            "analysis_workflow": analysis_workflow,
            "config_creation_date": datetime.now().strftime("%Y-%m-%d %H:%M"),
        },
        custom_filters={"umi_min_reads": umi_min_reads if umi_min_reads else None},
        reference=references,
        singularity={
            "image": Path(balsamic_cache, cache_version, "containers").as_posix()
        },
        background_variants=background_variants,
        samples=get_sample_list(
            tumor_sample_name=tumor_sample_name,
            normal_sample_name=normal_sample_name,
            tumor_fastq_path=tumor_fastq_path,
            normal_fastq_path=normal_fastq_path,
        ),
        vcf=VCF_DICT,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        bioinfo_tools_version=get_bioinfo_tools_version(
            bioinfo_tools=BIOINFO_TOOL_ENV,
            container_conda_env_path=CONTAINERS_DIR,
        ),
        panel=(
            {
                "exome": exome,
                "capture_kit": panel_bed,
                "chrom": get_panel_chrom(panel_bed),
                "pon_cnn": pon_cnn,
            }
            if panel_bed
            else None
        ),
    ).model_dump(by_alias=True, exclude_none=True)
    LOG.info("Balsamic config model instantiated successfully")

    config_path = Path(analysis_dir, case_id, case_id + ".json").as_posix()
    write_json(json_obj=config_collection_dict, path=config_path)
    LOG.info(f"Config file saved successfully - {config_path}")

    generate_graph(config_collection_dict, config_path)
    LOG.info(f"BALSAMIC Workflow has been configured successfully!")
