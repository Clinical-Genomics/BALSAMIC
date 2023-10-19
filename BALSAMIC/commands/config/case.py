"""Balsamic config case CLI."""
import logging
from pathlib import Path
from typing import Dict

import click

from BALSAMIC.commands.options import (
    OPTION_GENOME_VERSION,
    OPTION_ADAPTER_TRIM,
    OPTION_ANALYSIS_DIR,
    OPTION_ANALYSIS_WORKFLOW,
    OPTION_BACKGROUND_VARIANTS,
    OPTION_BALSAMIC_CACHE,
    OPTION_CADD_ANNOTATIONS,
    OPTION_CANCER_GERMLINE_SNV_OBSERVATIONS,
    OPTION_CANCER_SOMATIC_SNV_OBSERVATIONS,
    OPTION_CANCER_SOMATIC_SV_OBSERVATIONS,
    OPTION_CASE_ID,
    OPTION_CLINICAL_SNV_OBSERVATIONS,
    OPTION_CLINICAL_SV_OBSERVATIONS,
    OPTION_FASTQ_PATH,
    OPTION_GENDER,
    OPTION_GENOME_INTERVAL,
    OPTION_GENS_COV_PON,
    OPTION_GNOMAD_AF5,
    OPTION_NORMAL_SAMPLE_NAME,
    OPTION_PANEL_BED,
    OPTION_PON_CNN,
    OPTION_QUALITY_TRIM,
    OPTION_SWEGEN_SNV,
    OPTION_SWEGEN_SV,
    OPTION_TUMOR_SAMPLE_NAME,
    OPTION_UMI,
    OPTION_UMI_TRIM_LENGTH,
    OPTION_CACHE_VERSION,
)
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV, Gender, AnalysisWorkflow
from BALSAMIC.constants.cache import GenomeVersion
from BALSAMIC.constants.constants import FileType
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
from BALSAMIC.utils.io import write_json, read_json
from BALSAMIC.utils.utils import get_absolute_paths_dict

LOG = logging.getLogger(__name__)


@click.command("case", short_help="Create a sample config file from input sample data")
@OPTION_ADAPTER_TRIM
@OPTION_ANALYSIS_DIR
@OPTION_ANALYSIS_WORKFLOW
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
@OPTION_FASTQ_PATH
@OPTION_GENDER
@OPTION_GENOME_VERSION
@OPTION_GENOME_INTERVAL
@OPTION_GENS_COV_PON
@OPTION_GNOMAD_AF5
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
    cache_version: str,
    cadd_annotations: Path,
    cancer_germline_snv_observations: Path,
    cancer_somatic_snv_observations: Path,
    cancer_somatic_sv_observations: Path,
    case_id: str,
    clinical_snv_observations: Path,
    clinical_sv_observations: Path,
    fastq_path: Path,
    gender: Gender,
    genome_version: GenomeVersion,
    genome_interval: Path,
    gens_coverage_pon: Path,
    gnomad_min_af5: Path,
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
    references_path: Path = Path(balsamic_cache, cache_version, genome_version)
    references: Dict[str, Path] = get_absolute_paths_dict(
        base_path=references_path,
        data=read_json(Path(references_path, f"reference.{FileType.JSON}").as_posix()),
    )
    cadd_annotations_path = {"cadd_annotations": cadd_annotations}
    if cadd_annotations:
        references.update(cadd_annotations_path)

    if any([genome_interval, gens_coverage_pon, gnomad_min_af5]) and not all([genome_interval, gens_coverage_pon, gnomad_min_af5]):
        raise click.BadParameter(
            "All three arguments (genome_interval gens_coverage_pon, gnomad_min_af5) are required for GENS."
        )

        gens_ref_files = {
            "genome_interval": genome_interval,
            "gens_coverage_pon": gens_coverage_pon,
            "gnomad_min_af5": gnomad_min_af5,
        }

        references.update(
            {
                gens_file: path
                for gens_file, path in gens_ref_files.items()
                if path is not None
            }
        )

    variants_observations = {
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
        reference=references,
        singularity={
            "image": Path(balsamic_cache, cache_version, "containers").as_posix()
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
