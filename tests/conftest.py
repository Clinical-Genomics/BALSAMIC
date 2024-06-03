import copy
import os
import shutil
from datetime import datetime
from functools import partial
from pathlib import Path
from typing import Any, Dict, List
from unittest import mock

import pytest
from _pytest.tmpdir import TempPathFactory
from click.testing import CliRunner
from pydantic_core import Url

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.assets.scripts.preprocess_gens import cli as gens_preprocessing_cli
from BALSAMIC.commands.base import cli
from BALSAMIC.constants.analysis import (
    BIOINFO_TOOL_ENV,
    AnalysisWorkflow,
    PONWorkflow,
    RunMode,
)
from BALSAMIC.constants.cache import REFERENCE_FILES, DockerContainers, GenomeVersion
from BALSAMIC.constants.cluster import (
    QOS,
    ClusterAccount,
    ClusterConfigType,
    ClusterProfile,
    ClusterMailType,
)
from BALSAMIC.constants.constants import FileType
from BALSAMIC.constants.paths import CONSTANTS_DIR, FASTQ_TEST_INFO, TEST_DATA_DIR
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.models.cache import (
    AnalysisReferencesHg,
    CacheAnalysis,
    CacheConfig,
    References,
    ReferencesHg,
)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.models.scheduler import Scheduler
from BALSAMIC.models.snakemake import SingularityBindPath, SnakemakeExecutable
from BALSAMIC.utils.io import read_json, read_yaml, write_json
from .helpers import ConfigHelper, Map

MOCKED_OS_ENVIRON = "os.environ"


def fastq_patterns() -> list:
    """
    Returns a list of dicts containing different formatted fastq-file names to be used in parameterized tests.
    """
    fastq_test_info_path = Path(FASTQ_TEST_INFO).as_posix()
    fastq_test_info_dict = read_json(fastq_test_info_path)
    return fastq_test_info_dict["fastq_pattern_types"]


def fastq_pattern_ids() -> list:
    """
    Returns a list of IDs for the parameterized testing of different fastq file name formats.
    """
    fastq_test_info_path = Path(FASTQ_TEST_INFO).as_posix()
    fastq_test_info_dict = read_json(fastq_test_info_path)
    fastq_pattern_types = fastq_test_info_dict["fastq_pattern_types"]
    fastq_pattern_ids = ["FastqPattern{}".format(p["id"]) for p in fastq_pattern_types]
    return fastq_pattern_ids


@pytest.fixture(scope="session")
def empty_dir(tmp_path_factory: TempPathFactory) -> Path:
    """Return an empty directory path."""
    return tmp_path_factory.mktemp("empty_dir")


@pytest.fixture(scope="session")
def empty_file(empty_dir: Path) -> Path:
    """Return an empty directory path."""
    empty_file: Path = Path(empty_dir, "file.empty")
    empty_file.touch()
    return empty_file


@pytest.fixture(scope="session")
def sacct_file(case_id_tumor_only: str) -> Path:
    """Return sacct file path."""
    return Path(TEST_DATA_DIR, "logs", f"{case_id_tumor_only}.sacct")


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """
    Creates path for test data directory.
    """
    return TEST_DATA_DIR


@pytest.fixture(scope="session")
def load_test_fastq_data(test_data_dir) -> Dict:
    """Returns dict from loaded json containing strings of fastq-names."""
    fastq_test_info_path = Path(FASTQ_TEST_INFO).as_posix()
    return read_json(fastq_test_info_path)


@pytest.fixture(scope="session")
def pon_fastq_list(load_test_fastq_data) -> list:
    """Returns list of fastq names to be used in PON creation testing."""
    return load_test_fastq_data["pon_fastq_list"]


@pytest.fixture(scope="session")
def standard_samples_list(
    load_test_fastq_data,
) -> Dict[str, List]:
    """Returns a list of standard tumor normal sample dicts."""
    return load_test_fastq_data["samples_standard_fastq_names"]


@pytest.fixture(scope="session")
def standard_samples_list_pon(
    load_test_fastq_data,
) -> Dict[str, List]:
    """Returns a list of standard tumor normal sample dicts for PON."""
    return load_test_fastq_data["pon_samples_standard_fastq_names"]


@pytest.fixture(scope="session")
def tumor_fastq_names(
    load_test_fastq_data,
) -> Dict[str, List]:
    """Returns a list of standard tumor fastq-names."""
    return load_test_fastq_data["standard_fastq_names"]["tumor"]


@pytest.fixture(scope="session")
def normal_fastq_names(
    load_test_fastq_data,
) -> Dict[str, List]:
    """Returns a list of standard normal fastq-names."""
    return load_test_fastq_data["standard_fastq_names"]["normal"]


@pytest.fixture(scope="session")
def fastq_names_duplicate_assigned_fastq_patterns(
    load_test_fastq_data,
) -> Dict[str, List]:
    """Returns dict with list of fastq file names for testing of duplicate assigned fastq patterns."""
    return load_test_fastq_data["fastq_fails"]["duplicate_fastq_patterns"]


@pytest.fixture(scope="session")
def sample_list_duplicate_assigned_fastq_patterns_model(
    load_test_fastq_data,
) -> Dict[str, List]:
    """Returns List of sample-dicts with fastq-dicts with duplicate assigned fastq patterns."""
    return load_test_fastq_data["fastq_fails"]["duplicate_fastq_patterns_model"]


@pytest.fixture(scope="session")
def tumor_normal_fastq_info_correct(load_test_fastq_data) -> Dict[str, Dict]:
    """Mock tumor normal fastq info in sample_dict."""
    return load_test_fastq_data["test_fastq_info"]


@pytest.fixture(scope="session")
def valid_dnascope_variant() -> str:
    """Mock valid DNAscope variant."""
    return (
        "1\t100\trs1\tT\tC\t389.77\t.\tINFO\tGT:AD:DP:GQ:PL\t0/1:9,14:23:99:418,0,257"
    )


@pytest.fixture(scope="session")
def invalid_dnascope_variant_no_ad() -> str:
    """Mock invalid DNAscope variant without any read support."""
    return "1\t200\t.\tCAAA\tCAAAA,C\t0.00\tLowQual\tINFO\tGT:AD:DP:GQ:PL\t0/0:0,0,0:0:0:0,0,0,3,3,19"


@pytest.fixture(scope="session")
def invalid_dnascope_variant_illegal_chrom() -> str:
    """Mock invalid DNAscope variant with non-standard chromosome."""
    return (
        "25\t100\trs1\tT\tC\t389.77\t.\tINFO\tGT:AD:DP:GQ:PL\t0/1:9,14:23:99:418,0,257"
    )


@pytest.fixture(scope="session", name="session_tmp_path")
def fixture_session_tmp_path(tmp_path_factory: TempPathFactory) -> Path:
    """Return a non-existent files directory path."""
    return tmp_path_factory.mktemp("session_tests")


@pytest.fixture(scope="session")
def tumor_sample_name() -> str:
    """Create mock name for tumor sample."""
    return "ACC1"


@pytest.fixture(scope="session")
def normal_sample_name() -> str:
    """Create mock name for normal sample."""
    return "ACC2"


@pytest.fixture(scope="session")
def case_id_tumor_only() -> str:
    """Create mock case-id for TGA tumor-only."""
    return "sample_tumor_only"


@pytest.fixture(scope="session")
def case_id_tumor_only_dummy_vep() -> str:
    """Mock TGA tumor-only case ID for dummy vep file testing."""
    return "sample_tumor_only_dummy_vep"


@pytest.fixture(scope="session")
def case_id_tumor_only_qc() -> str:
    """Mock TGA tumor-only case ID for QC workflow."""
    return "sample_tumor_only_qc"


@pytest.fixture(scope="session")
def case_id_tumor_only_pon_cnn() -> str:
    """Mock TGA tumor-only case ID for testing with PON CNN file."""
    return "sample_tumor_only_pon_cnn"


@pytest.fixture(scope="session")
def case_id_pon() -> str:
    """
    Creates mock case-id for PON creation workflow
    """
    return "sample_pon_creation"


@pytest.fixture(scope="session")
def case_id_gens_pon() -> str:
    """
    Creates mock case-id for PON creation workflow
    """
    return "genscreation"


def case_id_tumor_only_pon() -> str:
    """Create mock case-id for TGA PON tumor-only."""
    return "sample_tumor_only_pon"


@pytest.fixture(scope="session")
def case_id_tumor_only_umi() -> str:
    """Creates mock case-id for TGA tumor-only UMI workflow."""
    return "sample_tumor_only_umi"


@pytest.fixture(scope="session")
def case_id_tumor_normal_umi() -> str:
    """Creates mock case-id for TGA tumor-only UMI workflow."""
    return "sample_tumor_normal_umi"


@pytest.fixture(scope="session")
def case_id_tumor_normal_fastqdir() -> str:
    """Mock case ID for dummy fastqdir."""
    return "sample_tumor_normal_fastqdir"


@pytest.fixture(scope="session")
def case_id_tumor_normal() -> str:
    """Create mock case-id for TGA tumor-normal."""
    return "sample_tumor_normal"


@pytest.fixture(scope="session")
def case_id_tumor_normal_qc() -> str:
    """Mock TGA tumor-normal case ID for QC TGA test."""
    return "sample_tumor_normal_qc"


@pytest.fixture(scope="session")
def case_id_tumor_normal_qc_wgs() -> str:
    """Mock TGA tumor-normal case ID for QC WGS test."""
    return "sample_tumor_normal_qc_wgs"


@pytest.fixture(scope="session")
def case_id_tumor_only_wgs() -> str:
    """Create mock case-id for WGS tumor-only."""
    return "sample_tumor_only_wgs"


@pytest.fixture(scope="session")
def case_id_tumor_normal_wgs() -> str:
    """Create mock case-id for WGS tumor-normal."""
    return "sample_tumor_normal_wgs"


@pytest.fixture
def cli_runner():
    """Run click for command line interface testing."""
    runner = CliRunner()
    return runner


@pytest.fixture
def invoke_cli(cli_runner):
    """Invoke cli commands with options."""
    return partial(cli_runner.invoke, cli)


@pytest.fixture
def invoke_gens_cli(cli_runner):
    """Invoke cli commands with options."""
    return partial(cli_runner.invoke, gens_preprocessing_cli)


@pytest.fixture(scope="session")
def environ():
    """Create operating system's environment object."""
    return "os.environ"


@pytest.fixture(scope="session")
def cluster_analysis_config_path() -> str:
    """Return cluster analysis configuration file."""
    return Path(
        CONSTANTS_DIR, f"{ClusterConfigType.ANALYSIS}.{FileType.JSON}"
    ).as_posix()


@pytest.fixture(scope="session")
def reference():
    """Return a dictionary for reference json model."""
    return {
        "reference_genome": "genome/human_g1k_v37_decoy.fasta",
        "dbsnp": "variants/dbsnp_grch37_b138.vcf.gz",
        "vcf_1kg": "variants/1k_genome_wgs_p1_v3_all_sites.vcf.gz",
        "hc_vcf_1kg": "variants/1kg_phase1_snps_high_confidence_b37.vcf.gz",
        "known_indel_1kg": "variants/1kg_known_indels_b37.vcf.gz",
        "mills_1kg": "variants/mills_1kg_index.vcf.gz",
        "gnomad_variant": "variants/gnomad.genomes.r2.1.1.sites.vcf.bgz",
        "cosmic": "variants/cosmic_coding_muts_v89.vcf.gz",
        "vep_dir": "vep/",
        "refgene_flat": "genome/refseq.flat",
        "refgene_txt": "genome/refGene.txt",
        "wgs_calling_regions": "genome/wgs_calling_regions.v1",
        "genome_chrom_size": "genome/hg19.chrom.sizes",
        "refgene_bed": "genome/refseq.flat.bed",
        "rank_score": "genome/cancer_rank_model_-v0.1-.ini",
        "access_regions": "genome/access-5k-mappable.hg19.bed",
        "delly_exclusion": "genome/delly_exclusion.tsv",
        "delly_exclusion_converted": "genome/delly_exclusion_converted.tsv",
        "delly_mappability": "genome/delly_mappability.gz",
        "delly_mappability_gindex": "genome/delly_mappability.gz.gzi",
        "delly_mappability_findex": "genome/delly_mappability.fai",
        "ascat_gc_correction": "genome/GRCh37_SnpGcCorrections.tsv",
        "ascat_chr_y_loci": "genome/GRCh37_Y.loci",
        "clinvar": "genome/clinvar.vcf.gz",
        "somalier_sites": "variants/GRCh37.somalier.sites.vcf.gz",
        "cadd_snv": "variants/hg19.cadd_snv.tsv.gz",
        "cadd_annotations": "cadd/",
        "simple_repeat": "genome/simpleRepeat.txt.gz",
    }


@pytest.fixture(scope="session")
def reference_panel_dir_path(test_data_dir: str) -> str:
    """Return path for reference panel directory."""
    return Path(test_data_dir, "references", "panel").as_posix()


@pytest.fixture(scope="session")
def reference_variants_dir_path(test_data_dir: str) -> str:
    """Return path for reference variants directory."""
    return Path(test_data_dir, "references", "variants").as_posix()


@pytest.fixture(scope="session")
def config_path(test_data_dir: str) -> str:
    """Created path for config json file."""
    return Path(test_data_dir, f"config.{FileType.JSON}").as_posix()


@pytest.fixture(scope="session")
def config_dict(config_path: str) -> str:
    """Read and return config from json."""
    return read_json(config_path)


@pytest.fixture(scope="session")
def config_dict_w_singularity(config_dict: str, balsamic_cache: str) -> str:
    """Read and return config from json with singularity image path."""
    modify_dict = copy.deepcopy(config_dict)
    modify_dict["singularity"] = {
        "image": f"{balsamic_cache}/{balsamic_version}/containers"
    }
    return modify_dict


@pytest.fixture(scope="session")
def pon_config_path(test_data_dir: str) -> str:
    """Created path for PON config json file."""
    return Path(test_data_dir, f"config_pon.{FileType.JSON}").as_posix()


@pytest.fixture(scope="session")
def pon_config_dict(pon_config_path: str) -> str:
    """Read and return PON config from json."""
    return read_json(pon_config_path)


@pytest.fixture(scope="session")
def pon_config_dict_w_singularity(pon_config_dict: str, balsamic_cache: str) -> str:
    """Read and return PON config from json with singularity image path."""
    modify_pon_config_dict = copy.deepcopy(pon_config_dict)
    modify_pon_config_dict["singularity"] = {
        "image": f"{balsamic_cache}/{balsamic_version}/containers"
    }
    return modify_pon_config_dict


@pytest.fixture(scope="session")
def cadd_annotations(test_data_dir: str) -> str:
    """Return path for CADD annotations."""
    return Path(test_data_dir, "references", "cadd").as_posix()


@pytest.fixture(scope="session")
def vcf_file_path(test_data_dir: str) -> str:
    """Return path for minimal VCF."""
    return Path(test_data_dir, "vcfs", "SNV.germline.sample.dnascope.vcf").as_posix()


@pytest.fixture(scope="session")
def vcf_file_gz_path(test_data_dir: str) -> str:
    """Return path for minimal gzipped VCF."""
    return Path(test_data_dir, "vcfs", "SNV.germline.sample.dnascope.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def gens_cov_pon_file(test_data_dir: str) -> str:
    """Return path for dummy GENS male PON file."""
    return Path(
        test_data_dir, "references", "gens", "grch37_gens_male_pon_100bp.hdf5"
    ).as_posix()


@pytest.fixture(scope="session")
def gens_min_5_af_gnomad_file(test_data_dir: str) -> str:
    """Return path for dummy GENS minimum af 5 gnomad file."""
    return Path(
        test_data_dir, "references", "gens", "gnomad.genomes.r2.1.1.sites_0.05AF.vcf.gz"
    ).as_posix()


@pytest.fixture(scope="session")
def gens_hg19_interval_list(test_data_dir: str) -> str:
    """Return path for dummy hg19 genome 100bp interval list used in GENS."""
    return Path(
        test_data_dir,
        "references",
        "gens",
        "grch37_gens_targets_preprocessed_100bp.interval_list",
    ).as_posix()


@pytest.fixture(scope="session")
def gens_dummy_gnomad_baf_bed(test_data_dir: str) -> str:
    """Return path expected dummy result-file created from GENS pre-processing test."""
    return Path(
        test_data_dir,
        "gens_files",
        "dummy.baf.bed",
    ).as_posix()


@pytest.fixture(scope="session")
def gens_dummy_gnomad_vcf(test_data_dir: str) -> str:
    """Return path dummy vcf called in given gnomad for GENS pre-processing test."""
    return Path(
        test_data_dir,
        "gens_files",
        "SNV.germline.dummy.dnascope_gnomad_af5.vcf",
    ).as_posix()


@pytest.fixture(scope="session")
def gens_dummy_cov_bed(test_data_dir: str) -> str:
    """Return path expected dummy result-file created from GENS pre-processing test."""
    return Path(
        test_data_dir,
        "gens_files",
        "dummy.cov.bed",
    ).as_posix()


@pytest.fixture(scope="session")
def gens_dummy_denoised_cov(test_data_dir: str) -> str:
    """Return path dummy coverage file for GENS pre-processing test."""
    return Path(
        test_data_dir,
        "gens_files",
        "dummy.denoisedCR.tsv",
    ).as_posix()


@pytest.fixture(scope="session")
def panel_bed_file(reference_panel_dir_path: str) -> str:
    """Return path for panel bed file."""
    return Path(reference_panel_dir_path, "panel.bed").as_posix()


@pytest.fixture(scope="session")
def background_variant_file(reference_panel_dir_path: str) -> str:
    """Return path for background variants for TGA."""
    return Path(reference_panel_dir_path, "background_variants.txt").as_posix()


@pytest.fixture(scope="session")
def pon_cnn_path(reference_panel_dir_path: str) -> str:
    """Creates path for Panel Of Normal (PON), cnn file for cnvkit."""
    return Path(reference_panel_dir_path, "test_panel_ponn.cnn").as_posix()


@pytest.fixture(scope="session")
def clinical_snv_observations_path(reference_variants_dir_path: str) -> str:
    """Return path for clinical SNVs from loqusDB."""
    return Path(reference_variants_dir_path, "clinical_snv_variants.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def cancer_germline_snv_observations_path(reference_variants_dir_path: str) -> str:
    """Return path of cancer germline SNVs from loqusDB."""
    return Path(
        reference_variants_dir_path, "cancer_germline_snv_variants.vcf.gz"
    ).as_posix()


@pytest.fixture(scope="session")
def cancer_somatic_snv_observations_path(reference_variants_dir_path: str) -> str:
    """Return path for somatic SNVs from loqusDB."""
    return Path(
        reference_variants_dir_path, "cancer_somatic_snv_variants.vcf.gz"
    ).as_posix()


@pytest.fixture(scope="session")
def clinical_sv_observations_path(reference_variants_dir_path: str) -> str:
    """Return path for clinical SV observations from loqusDB."""
    return Path(reference_variants_dir_path, "clinical_sv_variants.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def somatic_sv_observations_path(reference_variants_dir_path: str) -> str:
    """Return path for somatic SV observations from loqusDB."""
    return Path(reference_variants_dir_path, "somatic_sv_variants.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def swegen_snv_frequency_path(reference_variants_dir_path: str) -> str:
    """Return path for Swegen SNVs."""
    return Path(reference_variants_dir_path, "swegen_snv.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def swegen_sv_frequency_path(reference_variants_dir_path: str) -> str:
    """Create path for Swegen SVs."""
    return Path(reference_variants_dir_path, "swegen_sv.vcf.gz").as_posix()


@pytest.fixture(scope="session", name="invalid_json_file")
def fixture_invalid_json_file(session_tmp_path: Path) -> Path:
    """Return a non-existent json file path."""
    return Path(session_tmp_path, f"invalid_file.{FileType.JSON}")


@pytest.fixture(scope="session", name="json_file")
def fixture_json_file(session_tmp_path: Path) -> Path:
    """Return a mocked json file path."""
    return Path(session_tmp_path, f"write_json.{FileType.JSON}")


@pytest.fixture(scope="session", name="config_json")
def fixture_config_json() -> str:
    """Return Balsamic analysis config json file name."""
    return f"config.{FileType.JSON}"


@pytest.fixture(scope="session", name="reference_graph")
def fixture_reference_graph() -> str:
    """Return Balsamic reference graph pdf file name."""
    return "reference_graph.pdf"


@pytest.fixture(scope="session")
def sentieon_license(tmp_path_factory):
    """Create Sentieon's license path"""
    sentieon_license_dir = tmp_path_factory.mktemp("sentieon_licence")
    sentieon_license_path = sentieon_license_dir / "license_file.lic"
    sentieon_license_path.touch()

    return sentieon_license_path.as_posix()


@pytest.fixture(scope="session")
def sentieon_install_dir(tmp_path_factory):
    """Create install directory for Sentieon tools"""
    sentieon_install_dir = tmp_path_factory.mktemp("sentieon_install_dir")
    Path(sentieon_install_dir / "bin").mkdir(exist_ok=True)
    sentieon_executable = sentieon_install_dir / "bin" / "sentieon"
    sentieon_executable.touch()

    return sentieon_install_dir.as_posix()


@pytest.fixture()
def no_write_perm_path(tmp_path_factory) -> str:
    """Return path with no write permissions."""
    bad_perm_path: Path = tmp_path_factory.mktemp("bad_perm_path")
    bad_perm_path.chmod(0o444)
    return bad_perm_path.as_posix()


@pytest.fixture(scope="session")
def references_dir(test_data_dir) -> Path:
    """Return a references directory path."""
    return Path(test_data_dir, "references")


@pytest.fixture(scope="session")
def purity_csv_path(test_data_dir) -> Path:
    """Return pureCN purity CSV path."""
    return Path(test_data_dir, "cnv_report", "CNV.somatic.case_id.purecn.purity.csv")


@pytest.fixture(scope="session")
def cnv_plot_path(test_data_dir) -> Path:
    """Return AscatNgs CNV plot path."""
    return Path(
        test_data_dir,
        "cnv_report",
        "CNV.somatic.sample_tumor_normal_wgs.ascat.ASPCF.png",
    )


@pytest.fixture(scope="session")
def cnv_statistics_path(test_data_dir) -> Path:
    """Return CNV sample statistics path."""
    return Path(
        test_data_dir,
        "cnv_report",
        "CNV.somatic.sample_tumor_normal_wgs.ascat.samplestatistics.txt",
    )


@pytest.fixture(scope="session")
def balsamic_cache(
    tmp_path_factory: TempPathFactory, reference: Dict[str, Path], references_dir: Path
) -> str:
    """Create and return the path for balsamic-cache."""
    balsamic_cache: Path = tmp_path_factory.mktemp("balsamic_cache")
    # Mocked containers directory
    container: Path = Path(balsamic_cache, balsamic_version, "containers")
    container.mkdir(parents=True, exist_ok=True)
    # Mocked cache directory
    hg19_reference: Path = Path(balsamic_cache, balsamic_version, "hg19")
    shutil.copytree(references_dir, hg19_reference)
    reference_json: Path = Path(hg19_reference, f"reference.{FileType.JSON}")
    reference_json.touch()
    write_json(json_obj=reference, path=reference_json.as_posix())
    return balsamic_cache.as_posix()


@pytest.fixture(scope="session")
def analysis_dir(tmp_path_factory: TempPathFactory) -> str:
    """Create and return the directory where the case analysis will be saved."""
    analysis_dir = tmp_path_factory.mktemp("analysis", numbered=False)
    return analysis_dir.as_posix()


@pytest.fixture(scope="session", params=fastq_patterns(), ids=fastq_pattern_ids())
def fastq_dir(case_id_tumor_normal_fastqdir: str, analysis_dir: str, request):
    """Mock directory with tumor and normal FASTQs."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal_fastqdir, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    fastq_test_dict = request.param

    for fastq in fastq_test_dict["tumor"]:
        Path(fastq_dir, fastq).touch()

    for fastq in fastq_test_dict["normal"]:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()

    for fastq in fastq_test_dict["tumor"]:
        Path.unlink(fastq_dir / fastq)

    for fastq in fastq_test_dict["normal"]:
        Path.unlink(fastq_dir / fastq)


@pytest.fixture(scope="session")
def fastq_dir_tumor_only(
    analysis_dir: str, case_id_tumor_only: str, tumor_fastq_names: List[str]
) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-only.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_symlinked(
    tumor_fastq_names: List[str], session_tmp_path: Path, reference_file: Path
) -> Path:
    """Return directory containing symlinked FASTQs for tumor-only."""
    fastq_dir: Path = Path(session_tmp_path, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    for fastq in tumor_fastq_names:
        Path(session_tmp_path, fastq).touch()
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).symlink_to(reference_file)
    return fastq_dir


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_pon_cnn(
    analysis_dir: str, case_id_tumor_only_pon_cnn: str, tumor_fastq_names: List[str]
) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-only w pon cnn.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_pon_cnn, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_wgs(
    analysis_dir: str, case_id_tumor_only_wgs: str, tumor_fastq_names: List[str]
) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-only WGS.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_wgs, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_normal(
    analysis_dir: str,
    case_id_tumor_normal: str,
    tumor_fastq_names: List[str],
    normal_fastq_names: List[str],
) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-normal.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    for fastq in normal_fastq_names:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_normal_wgs(
    analysis_dir: str,
    case_id_tumor_normal_wgs: str,
    tumor_fastq_names: List[str],
    normal_fastq_names: List[str],
) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-normal WGS.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal_wgs, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    for fastq in normal_fastq_names:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_normal_qc(
    analysis_dir: str,
    case_id_tumor_normal_qc: str,
    tumor_fastq_names: List[str],
    normal_fastq_names: List[str],
) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-normal QC workflow.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal_qc, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    for fastq in normal_fastq_names:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_qc(
    analysis_dir: str, case_id_tumor_only_qc: str, tumor_fastq_names: List[str]
) -> str:
    """Creates and returns the directory containing the FASTQs."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_qc, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_dummy_vep(
    analysis_dir: str, case_id_tumor_only_dummy_vep: str, tumor_fastq_names: List[str]
) -> str:
    """Creates and returns the directory containing the FASTQs."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_dummy_vep, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in tumor_fastq_names:
        Path(fastq_dir, fastq).touch()

    vep_dir: Path = Path(analysis_dir, case_id_tumor_only_dummy_vep, "analysis", "vep")
    vep_dir.mkdir(parents=True, exist_ok=True)
    vep_test_file = (
        "SNV.somatic.sample_tumor_only.vardict.research.filtered.pass.vcf.gz"
    )
    Path(vep_dir, vep_test_file).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_pon(analysis_dir: str, case_id_pon: str, pon_fastq_list: list) -> str:
    """
    Creates and returns the directory containing the FASTQs for PON creation workflow.
    """

    fastq_dir: Path = Path(analysis_dir, case_id_pon, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    for fastq in pon_fastq_list:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_gens_pon(
    analysis_dir: str, case_id_gens_pon: str, pon_fastq_list: list
) -> str:
    """
    Creates and returns the directory containing the FASTQs for PON creation workflow.
    """

    fastq_dir: Path = Path(analysis_dir, case_id_gens_pon, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    for fastq in pon_fastq_list:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def empty_fastq_dir(analysis_dir: str, case_id_tumor_normal: str) -> str:
    """
    Creates and returns an empty FASTQ directory.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal, "fastq_empty")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session", params=fastq_patterns(), ids=fastq_pattern_ids())
def fastq_dir_tumor_normal_parameterize(
    analysis_dir: str, case_id_tumor_normal: str, request
) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-normal once for each fastq-name structure.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal, "fastq_parameterize")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    fastq_test_dict = request.param
    for fastq in fastq_test_dict["tumor"]:
        Path(fastq_dir, fastq).touch()

    for fastq in fastq_test_dict["normal"]:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()

    for fastq in fastq_test_dict["tumor"]:
        Path.unlink(fastq_dir / fastq)

    for fastq in fastq_test_dict["normal"]:
        Path.unlink(fastq_dir / fastq)


@pytest.fixture(scope="session")
def fastq_dir_tumor_duplicate_fastq_patterns(
    analysis_dir: str,
    case_id_tumor_normal: str,
    fastq_names_duplicate_assigned_fastq_patterns: Dict,
) -> str:
    """Creates and returns the directory containing the FASTQs to test duplicate fastq-patterns."""
    fastq_dir: Path = Path(
        analysis_dir, case_id_tumor_normal, "fastq_duplicate_assigned_fastq_patterns"
    )
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    for fastq in fastq_names_duplicate_assigned_fastq_patterns["tumor"]:
        Path(fastq_dir, fastq).touch()

    yield fastq_dir.as_posix()


@pytest.fixture(scope="session")
def config_dict_w_fastqs(
    analysis_dir: str,
    case_id_tumor_normal: str,
    config_dict_w_singularity: str,
    standard_samples_list: List[Dict],
) -> str:
    """Change samples-list in config and create test fastq-files."""

    fastq_dir: Path = Path(
        analysis_dir,
        case_id_tumor_normal,
        "fastq_standard_names",
    )
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Change fastq_path to be the newly created test fastq dir
    modified_config = copy.deepcopy(config_dict_w_singularity)
    modified_config["analysis"]["fastq_path"] = fastq_dir.as_posix()

    # Create analysis_dirs
    analysis_sub_dirs = ["script", "log", "result", "benchmark"]
    for analysis_sub_dir in analysis_sub_dirs:
        analysis_sub_dir_path: Path = Path(
            analysis_dir, case_id_tumor_normal, analysis_sub_dir
        )
        analysis_sub_dir_path.mkdir(parents=True, exist_ok=True)
        modified_config["analysis"][analysis_sub_dir] = analysis_sub_dir_path.as_posix()

    # Fill the fastq path folder with the test fastq-files
    samples = standard_samples_list
    for sample_dict in samples:
        for fastq_pattern, values in sample_dict["fastq_info"].items():
            values["fwd"] = fastq_dir.joinpath(values["fwd"]).as_posix()
            values["rev"] = fastq_dir.joinpath(values["rev"]).as_posix()
            # Create dummy fastq files
            Path(values["fwd"]).touch()
            Path(values["rev"]).touch()

    # Modify input config sample list to correspond to current test sample list
    modified_config["samples"] = samples

    return modified_config


@pytest.fixture(scope="session")
def pon_config_dict_w_fastq(
    analysis_dir: str,
    case_id_pon,
    pon_config_dict_w_singularity: str,
    balsamic_cache: str,
    standard_samples_list_pon: List[Dict],
) -> str:
    """Create fastqs and modify pon config to contain created fastq paths."""
    fastq_dir: Path = Path(
        analysis_dir,
        case_id_pon,
        "fastq_standard_names_pon",
    )
    fastq_dir.mkdir(parents=True, exist_ok=True)

    pon_config_w_fastq = copy.deepcopy(pon_config_dict_w_singularity)
    pon_config_w_fastq["analysis"]["fastq_path"] = fastq_dir.as_posix()

    # Create analysis_dirs and modify config
    analysis_sub_dirs = ["script", "log", "result", "benchmark"]
    for analysis_sub_dir in analysis_sub_dirs:
        analysis_sub_dir_path: Path = Path(analysis_dir, case_id_pon, analysis_sub_dir)
        analysis_sub_dir_path.mkdir(parents=True, exist_ok=True)
        pon_config_w_fastq["analysis"][
            analysis_sub_dir
        ] = analysis_sub_dir_path.as_posix()

    # Fill the fastq path folder with the test fastq-files
    samples_list = standard_samples_list_pon
    for sample_dict in samples_list:
        for fastq_pattern, values in sample_dict["fastq_info"].items():
            fwd_fastq_path = f"{fastq_dir}/{os.path.basename(values['fwd'])}"
            rev_fastq_path = f"{fastq_dir}/{os.path.basename(values['rev'])}"
            values["fwd"] = fwd_fastq_path
            values["rev"] = rev_fastq_path
            Path(fwd_fastq_path).touch()
            Path(rev_fastq_path).touch()

    # Modify input config sample list to correspond to current test sample list
    pon_config_w_fastq["samples"] = samples_list

    return pon_config_w_fastq


@pytest.fixture(scope="session")
def config_w_fastq_dir_for_duplicate_fastq_patterns_model(
    analysis_dir: str,
    case_id_tumor_normal: str,
    sample_list_duplicate_assigned_fastq_patterns_model: List[Dict],
    config_dict_w_singularity: Dict,
) -> str:
    """Creates and returns the directory containing the FASTQs to test duplicate fastq-patterns to test model."""

    fastq_dir: Path = Path(
        analysis_dir,
        case_id_tumor_normal,
        "fastq_duplicate_assigned_fastq_patterns_model",
    )
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Change fastq_path to be the newly created test fastq dir
    modified_config = copy.deepcopy(config_dict_w_singularity)
    modified_config["analysis"]["fastq_path"] = fastq_dir.as_posix()

    # Create analysis_dirs and modify config
    analysis_sub_dirs = ["script", "log", "result", "benchmark"]
    for analysis_sub_dir in analysis_sub_dirs:
        analysis_sub_dir_path: Path = Path(
            analysis_dir, case_id_tumor_normal, analysis_sub_dir
        )
        analysis_sub_dir_path.mkdir(parents=True, exist_ok=True)
        modified_config["analysis"][analysis_sub_dir] = analysis_sub_dir_path.as_posix()

    # Fill the fastq path folder with the test fastq-files
    samples = sample_list_duplicate_assigned_fastq_patterns_model
    for sample_dict in samples:
        for fastq_pattern, values in sample_dict["fastq_info"].items():
            values["fwd"] = fastq_dir.joinpath(values["fwd"]).as_posix()
            values["rev"] = fastq_dir.joinpath(values["rev"]).as_posix()
            # Create dummy fastq files
            Path(values["fwd"]).touch()
            Path(values["rev"]).touch()

    # Modify input config sample list to correspond to current test sample list
    modified_config["samples"] = samples

    return modified_config


@pytest.fixture(scope="session")
def config_tumor_normal_extrafile(
    analysis_dir: str,
    case_id_tumor_normal: str,
    config_dict_w_singularity: Dict,
) -> str:
    """Creates and returns the directory containing the FASTQs to test detection of unassigned fastq-files."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal, "fastq_extrafile")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the fastq path folder with the test fastq-files
    modified_config = copy.deepcopy(config_dict_w_singularity)

    # Change fastq_path to be the newly created test fastq dir
    modified_config["analysis"]["fastq_path"] = fastq_dir.as_posix()

    # Create analysis_dirs and modify config
    analysis_sub_dirs = ["script", "log", "result", "benchmark"]
    for analysis_sub_dir in analysis_sub_dirs:
        analysis_sub_dir_path: Path = Path(
            analysis_dir, case_id_tumor_normal, analysis_sub_dir
        )
        analysis_sub_dir_path.mkdir(parents=True, exist_ok=True)
        modified_config["analysis"][analysis_sub_dir] = analysis_sub_dir_path.as_posix()

    samples_list = modified_config["samples"]
    for sample_dict in samples_list:
        for fastq_pattern, values in sample_dict["fastq_info"].items():
            fwd_fastq_path = f"{fastq_dir}/{os.path.basename(values['fwd'])}"
            rev_fastq_path = f"{fastq_dir}/{os.path.basename(values['rev'])}"
            values["fwd"] = fwd_fastq_path
            values["rev"] = rev_fastq_path
            Path(fwd_fastq_path).touch()
            Path(rev_fastq_path).touch()
    modified_config["samples"] = samples_list

    # Add extra files not assigned to dict
    extra_file1 = "ACC3fail_S1_L001_R1_001.fastq.gz"
    extra_file2 = "ACC3fail_S1_L001_R2_001.fastq.gz"
    Path(fastq_dir, extra_file1).touch()
    Path(fastq_dir, extra_file2).touch()

    # Returned modified dict
    return modified_config


@pytest.fixture(scope="session")
def snakemake_job_script(tmp_path_factory, tumor_normal_config):
    """Create a dummy snakemake jobscript"""
    script_dir = tmp_path_factory.mktemp("snakemake_script")
    snakemake_script_file = script_dir / "example_script.sh"
    snakemake_script = """#!/bin/sh
# properties = {"type": "single", "rule": "all", "local": false, "input": ["dummy_path"], "output": ["dummy_path"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 0, "cluster": {"name": "BALSAMIC.all.", "time": "00:15:00", "n": 1, "mail_type": "END", "partition": "core"}}
ls -l # dummy command
"""
    snakemake_script_file.touch()
    with open(snakemake_script_file, "w") as fn:
        fn.write(snakemake_script)

    return {"snakescript": str(snakemake_script_file)}


@pytest.fixture(name="helpers")
def fixture_config_helpers():
    """Return helper for case config files"""
    return ConfigHelper()


@pytest.fixture(scope="session")
def balsamic_model(
    config_dict_w_fastqs: Dict,
) -> ConfigModel:
    """Return ConfigModel parsed from static tumor normal config dict."""
    # Initialize balsamic model
    balsamic_config = ConfigModel.model_validate(config_dict_w_fastqs)
    return balsamic_config


@pytest.fixture(scope="session")
def balsamic_pon_model(
    pon_config_dict_w_fastq: Dict,
) -> ConfigModel:
    """Return ConfigModel parsed from static PON config dict."""
    # Initialize ConfigModel
    balsamic_pon_config = ConfigModel.model_validate(pon_config_dict_w_fastq)
    return balsamic_pon_config


@pytest.fixture(scope="session")
def config_case_cli(
    balsamic_cache: str,
    background_variant_file: str,
    cadd_annotations: str,
    swegen_snv_frequency_path: str,
    swegen_sv_frequency_path: str,
    clinical_snv_observations_path: str,
    clinical_sv_observations_path: str,
    somatic_sv_observations_path: str,
    cancer_germline_snv_observations_path: str,
    cancer_somatic_snv_observations_path: str,
) -> List[str]:
    """Return common config case CLI."""
    return [
        "--balsamic-cache",
        balsamic_cache,
        "--background-variants",
        background_variant_file,
        "--cadd-annotations",
        cadd_annotations,
        "--swegen-snv",
        swegen_snv_frequency_path,
        "--swegen-sv",
        swegen_sv_frequency_path,
        "--clinical-snv-observations",
        clinical_snv_observations_path,
        "--clinical-sv-observations",
        clinical_sv_observations_path,
        "--cancer-somatic-sv-observations",
        somatic_sv_observations_path,
        "--cancer-germline-snv-observations",
        cancer_germline_snv_observations_path,
        "--cancer-somatic-snv-observations",
        cancer_somatic_snv_observations_path,
    ]


@pytest.fixture(scope="session")
def tumor_only_config_qc(
    case_id_tumor_only_qc: str,
    analysis_dir: str,
    fastq_dir_tumor_only_qc: str,
    tumor_sample_name: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: list[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-only TGA."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only_qc,
                "--analysis-workflow",
                AnalysisWorkflow.BALSAMIC_QC,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only_qc,
                "--tumor-sample-name",
                tumor_sample_name,
                "-p",
                panel_bed_file,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir, case_id_tumor_only_qc, f"{case_id_tumor_only_qc}.{FileType.JSON}"
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_normal_config_qc(
    case_id_tumor_normal_qc: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal_qc: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: list[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-normal TGA QC workflow."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_normal_qc,
                "--analysis-workflow",
                AnalysisWorkflow.BALSAMIC_QC,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_normal_qc,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir,
        case_id_tumor_normal_qc,
        f"{case_id_tumor_normal_qc}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_normal_config_qc_wgs(
    case_id_tumor_normal_qc_wgs: str,
    analysis_dir: str,
    fastq_dir_tumor_normal_qc_wgs: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: List[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-normal WGS QC workflow."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_normal_qc_wgs,
                "--analysis-workflow",
                AnalysisWorkflow.BALSAMIC_QC,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_normal_qc_wgs,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir,
        case_id_tumor_normal_qc_wgs,
        f"{case_id_tumor_normal_qc_wgs}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_config(
    case_id_tumor_only: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: list[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-only TGA."""
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only,
                "--tumor-sample-name",
                tumor_sample_name,
                "-p",
                panel_bed_file,
            ]
            + config_case_cli,
        )
    return Path(
        analysis_dir,
        case_id_tumor_only,
        f"{case_id_tumor_only}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_normal_config(
    case_id_tumor_normal: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: list[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-normal TGA."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_normal,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_normal,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
                "-p",
                panel_bed_file,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir,
        case_id_tumor_normal,
        f"{case_id_tumor_normal}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_umi_config(
    case_id_tumor_only_umi: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: list[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-only TGA."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only_umi,
                "--analysis-workflow",
                AnalysisWorkflow.BALSAMIC_UMI,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only,
                "--tumor-sample-name",
                tumor_sample_name,
                "-p",
                panel_bed_file,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir,
        case_id_tumor_only_umi,
        f"{case_id_tumor_only_umi}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_normal_umi_config(
    case_id_tumor_normal_umi: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: list[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-normal TGA."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "-p",
                panel_bed_file,
                "--case-id",
                case_id_tumor_normal_umi,
                "--analysis-workflow",
                AnalysisWorkflow.BALSAMIC_UMI,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_normal,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
                "-p",
                panel_bed_file,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir,
        case_id_tumor_normal_umi,
        f"{case_id_tumor_normal_umi}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_wgs_config(
    case_id_tumor_only_wgs: str,
    analysis_dir: str,
    fastq_dir_tumor_only_wgs: str,
    tumor_sample_name: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: List[str],
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-only WGS."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only_wgs,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only_wgs,
                "--tumor-sample-name",
                tumor_sample_name,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir,
        case_id_tumor_only_wgs,
        f"{case_id_tumor_only_wgs}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_normal_wgs_config(
    case_id_tumor_normal_wgs: str,
    analysis_dir: str,
    fastq_dir_tumor_normal_wgs: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    config_case_cli: str,
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-normal WGS."""
    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_normal_wgs,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_normal_wgs,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
            ]
            + config_case_cli,
        )

    return Path(
        analysis_dir,
        case_id_tumor_normal_wgs,
        f"{case_id_tumor_normal_wgs}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_config_dummy_vep(
    case_id_tumor_only_dummy_vep: str,
    tumor_sample_name: str,
    balsamic_cache: str,
    analysis_dir: str,
    fastq_dir_tumor_only_dummy_vep: str,
    panel_bed_file: str,
    background_variant_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
) -> str:
    """Invoke balsamic config sample to create sample configuration file for tumor-only TGA with dummy VEP file."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only_dummy_vep,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only_dummy_vep,
                "-p",
                panel_bed_file,
                "--balsamic-cache",
                balsamic_cache,
                "--background-variants",
                background_variant_file,
                "--tumor-sample-name",
                tumor_sample_name,
            ],
        )
    return Path(
        analysis_dir,
        case_id_tumor_only_dummy_vep,
        f"{case_id_tumor_only_dummy_vep}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_pon_config(
    case_id_tumor_only_pon_cnn: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only_pon_cnn: str,
    panel_bed_file: str,
    pon_cnn_path: str,
    balsamic_cache: str,
    sentieon_license: str,
    sentieon_install_dir: str,
) -> str:
    """Invoke balsamic PON config sample to create sample configuration file for tumor-only TGA."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "case",
                "--case-id",
                case_id_tumor_only_pon_cnn,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only_pon_cnn,
                "-p",
                panel_bed_file,
                "--pon-cnn",
                pon_cnn_path,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
            ],
        )

    return Path(
        analysis_dir,
        case_id_tumor_only_pon_cnn,
        f"{case_id_tumor_only_pon_cnn}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def cnvkit_pon_creation_config(
    case_id_pon: str,
    analysis_dir: str,
    fastq_dir_pon: str,
    panel_bed_file: str,
    balsamic_cache: str,
    sentieon_license: str,
    sentieon_install_dir: str,
) -> str:
    """Invoke PON creation config configuration file for CNVkit PON workflow."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "pon",
                "--case-id",
                case_id_pon,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_pon,
                "-p",
                panel_bed_file,
                "--version",
                "v5",
                "--balsamic-cache",
                balsamic_cache,
                "--pon-workflow",
                PONWorkflow.CNVKIT,
            ],
        )

    return Path(
        analysis_dir, case_id_pon, f"{case_id_pon}_PON.{FileType.JSON}"
    ).as_posix()


@pytest.fixture(scope="session")
def gens_pon_creation_config(
    case_id_gens_pon: str,
    analysis_dir: str,
    fastq_dir_gens_pon: str,
    balsamic_cache: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    gens_hg19_interval_list: str,
) -> str:
    """Invoke PON creation config configuration file for GENS PON workflow."""

    with mock.patch.dict(
        MOCKED_OS_ENVIRON,
        {
            "SENTIEON_LICENSE": sentieon_license,
            "SENTIEON_INSTALL_DIR": sentieon_install_dir,
        },
    ):
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "config",
                "pon",
                "--case-id",
                case_id_gens_pon,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_gens_pon,
                "--version",
                "v5",
                "--balsamic-cache",
                balsamic_cache,
                "--pon-workflow",
                PONWorkflow.GENS_MALE,
                "--genome-interval",
                gens_hg19_interval_list,
            ],
        )

    return Path(
        analysis_dir, case_id_gens_pon, f"{case_id_gens_pon}_PON.{FileType.JSON}"
    ).as_posix()


@pytest.fixture(scope="session")
def sample_config(
    tumor_sample_name: str,
    normal_sample_name: str,
    tumor_normal_fastq_info_correct: dict,
):
    """Create and return sample config dict to test workflow utils"""
    sample_config = {
        "QC": {
            "picard_rmdup": "False",
            "adapter": "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
            "min_seq_length": "25",
            "quality_trim": "True",
            "adapter_trim": "False",
            "umi_trim": "True",
            "umi_trim_length": "5",
        },
        "analysis": {
            "case_id": "id1",
            "analysis_type": "paired",
            "analysis_dir": "tests/test_data/",
            "fastq_path": "tests/test_data/id1/fastq/",
            "script": "tests/test_data/id1/scripts/",
            "log": "tests/test_data/id1/logs/",
            "result": "tests/test_data/id1/analysis/",
            "config_creation_date": "yyyy-mm-dd xx",
            "BALSAMIC_version": "2.9.8",
            "dag": "tests/test_data/id1/id1_analysis.json_BALSAMIC_2.9.8_graph.pdf",
        },
        "reference": {
            "reference_genome": "tests/test_data/references/genome/human_g1k_v37_decoy.fasta",
            "genome_chrom_size": "tests/test_data/references/genome/hg19.chrom.sizes",
        },
        "vcf": VCF_DICT,
        "samples": tumor_normal_fastq_info_correct,
        "panel": {
            "capture_kit": "tests/test_data/references/panel/panel.bed",
            "pon_cnn": "tests/test_data/references/panel/test_panel_ponn.cnn",
        },
        "umiworkflow": "true",
    }

    return sample_config


@pytest.fixture(scope="session")
def analysis_path():
    """Return path for test analysis."""
    return "tests/test_data/qc_files/analysis"


@pytest.fixture(scope="session")
def multiqc_data_dir(analysis_path: str) -> Path:
    """Return path of tje MultiQC test data directory."""
    return Path(analysis_path, "qc", "multiqc_data")


@pytest.fixture(scope="session")
def multiqc_data_path(multiqc_data_dir: Path) -> str:
    """Return path of JSON for MultiQC test data."""
    return Path(multiqc_data_dir, "multiqc_data.json").as_posix()


@pytest.fixture(scope="session")
def multiqc_data_dict(multiqc_data_path: str) -> dict:
    """Read and Return test data from JASON of MultiQC test data."""
    return read_json(multiqc_data_path)


@pytest.fixture(scope="session")
def metrics_yaml_path(analysis_path: str) -> str:
    """Return path for Tumor-Only deliverable metrics from YAML."""
    return Path(
        analysis_path, "qc", "sample_tumor_only_metrics_deliverables.yaml"
    ).as_posix()


@pytest.fixture(scope="session")
def bcftools_counts_path(analysis_path: str) -> str:
    """Return path for svdb.clinical.filtered.pass.stats."""
    return Path(
        analysis_path, "vep", "SNV.somatic.case.svdb.clinical.filtered.pass.stats"
    ).as_posix()


@pytest.fixture(scope="session")
def qc_requested_metrics():
    """Return raw requested metrics."""
    return {
        "targeted": {
            "default": {
                "METRIC_1": {"condition": None},
                "METRIC_2": {"condition": {"norm": "gt", "threshold": 2}},
            },
            "panel_1": {
                "METRIC_3": {"condition": {"norm": "gt", "threshold": 3}},
            },
            "panel_2": {
                "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
                "METRIC_2": {"condition": {"norm": "gt", "threshold": 22}},
                "METRIC_4": {"condition": {"norm": "gt", "threshold": 4}},
            },
        },
        "wgs": {
            "METRIC_1": {"condition": {"norm": "gt", "threshold": 1}},
        },
    }


@pytest.fixture(scope="session")
def qc_extracted_metrics(metrics_yaml_path: str) -> dict:
    """Return extracted and formatted QC metrics."""
    return read_yaml(metrics_yaml_path)


@pytest.fixture(scope="function")
def snakemake_bcftools_filter_vardict_research_tumor_only(
    tumor_only_config_dummy_vep, helpers
):
    """bcftools_filter_vardict_research_tumor_only snakemake mock rule"""

    helpers.read_config(tumor_only_config_dummy_vep)
    vep_path = os.path.join(
        helpers.analysis_dir,
        helpers.case_id,
        "analysis",
        "vep",
        "{var_type}.somatic.{case_name}.vardict.research.filtered.pass.vcf.gz",
    )
    return Map(
        {
            "bcftools_filter_vardict_research_tumor_only": Map(
                {
                    "params": Map(
                        {
                            "housekeeper_id": {
                                "id": "sample_tumor_only_single",
                                "tags": "research",
                            }
                        }
                    ),
                    "output": Map(
                        {
                            "_names": Map({"vcf_pass_vardict": vep_path}),
                            "vcf_pass_vardict": vep_path,
                        }
                    ),
                    "rule": Map(
                        {
                            "name": "bcftools_filter_vardict_research_tumor_only",
                            "output": [
                                vep_path,
                            ],
                            "temp_output": set(),
                        }
                    ),
                }
            )
        }
    )


@pytest.fixture(scope="session", name="timestamp_now")
def fixture_timestamp_now() -> datetime:
    """Return a time stamp of current date in date time format."""
    return datetime.now()


@pytest.fixture(scope="session", name="cosmic_key")
def fixture_cosmic_key() -> str:
    """Return a mocked COSMIC key."""
    return "ZW1haWxAZXhhbXBsZS5jb206bXljb3NtaWNwYXNzd29yZAo="


@pytest.fixture(scope="session", name="develop_containers")
def fixture_develop_containers() -> Dict[str, str]:
    """Return a dictionary of docker hub containers for develop branch."""
    return {
        DockerContainers.ASCAT: "docker://clinicalgenomics/balsamic:develop-ascatNgs",
        DockerContainers.MULTIQC: "docker://clinicalgenomics/balsamic:develop-multiqc",
        DockerContainers.VCF2CYTOSURE: "docker://clinicalgenomics/balsamic:develop-vcf2cytosure",
        DockerContainers.PYTHON_3: "docker://clinicalgenomics/balsamic:develop-varcall_py3",
        DockerContainers.SOMALIER: "docker://clinicalgenomics/balsamic:develop-somalier",
        DockerContainers.CNVPYTOR: "docker://clinicalgenomics/balsamic:develop-cnvpytor",
        DockerContainers.ALIGN_QC: "docker://clinicalgenomics/balsamic:develop-align_qc",
        DockerContainers.ANNOTATE: "docker://clinicalgenomics/balsamic:develop-annotate",
        DockerContainers.PYTHON_27: "docker://clinicalgenomics/balsamic:develop-varcall_py27",
        DockerContainers.CNVKIT: "docker://clinicalgenomics/balsamic:develop-cnvkit",
        DockerContainers.COVERAGE_QC: "docker://clinicalgenomics/balsamic:develop-coverage_qc",
        DockerContainers.DELLY: "docker://clinicalgenomics/balsamic:develop-delly",
        DockerContainers.CADD: "docker://clinicalgenomics/balsamic:develop-cadd",
        DockerContainers.HTSLIB: "docker://clinicalgenomics/balsamic:develop-htslib",
        DockerContainers.PURECN: "docker://clinicalgenomics/balsamic:develop-purecn",
        DockerContainers.GATK: "docker://clinicalgenomics/balsamic:develop-gatk",
    }


@pytest.fixture(scope="session", name="cache_config_data")
def fixture_cache_config_data(
    cache_analysis: CacheAnalysis,
    develop_containers: Dict[str, str],
    cosmic_key: str,
    timestamp_now: datetime,
    session_tmp_path: Path,
) -> Dict[str, Any]:
    """Return mocked cache config data."""

    return {
        "analysis": cache_analysis,
        "references_dir": session_tmp_path,
        "genome_dir": session_tmp_path,
        "variants_dir": session_tmp_path,
        "vep_dir": session_tmp_path,
        "containers_dir": session_tmp_path,
        "genome_version": GenomeVersion.HG19,
        "cosmic_key": cosmic_key,
        "bioinfo_tools": BIOINFO_TOOL_ENV,
        "containers": develop_containers,
        "references": REFERENCE_FILES[GenomeVersion.HG19],
        "references_date": timestamp_now.strftime("%Y-%m-%d %H:%M"),
    }


@pytest.fixture(scope="session", name="cache_config")
def fixture_cache_config(cache_config_data: Dict[str, dict]) -> CacheConfig:
    """Return mocked cache config model."""
    cache_config: CacheConfig = CacheConfig(**cache_config_data)
    for reference in cache_config.references:
        reference_file: Path = Path(reference[1].file_path)
        reference_file.parent.mkdir(parents=True, exist_ok=True)
        reference_file.touch()
    return cache_config


@pytest.fixture(scope="session", name="cache_analysis_data")
def fixture_cache_analysis_data(case_id_tumor_only: str) -> Dict[str, str]:
    """Return mocked cache analysis data."""
    return {"case_id": case_id_tumor_only}


@pytest.fixture(scope="session", name="cache_analysis")
def fixture_cache_analysis(cache_analysis_data: Dict[str, str]) -> CacheAnalysis:
    """Return mocked cache analysis model."""
    return CacheAnalysis(**cache_analysis_data)


@pytest.fixture(scope="session", name="refgene_bed_file")
def fixture_refgene_bed_file(session_tmp_path: Path) -> Path:
    """Return dummy refseq's gene bed file."""
    refgene_bed_file: Path = Path(session_tmp_path, "genome", "refGene.flat.bed")
    refgene_bed_file.touch()
    return refgene_bed_file


@pytest.fixture(scope="session", name="refgene_flat_file")
def fixture_refgene_flat_file(session_tmp_path: Path) -> Path:
    """Return dummy refseq's gene flat file."""
    refgene_flat_file: Path = Path(session_tmp_path, "genome", "refGene.flat")
    refgene_flat_file.touch()
    return refgene_flat_file


@pytest.fixture(scope="session", name="analysis_references_data")
def fixture_analysis_references_data(
    cache_config: CacheConfig,
    refgene_bed_file: Path,
    refgene_flat_file: Path,
) -> Dict[str, Path]:
    """Return analysis references model data."""
    return {
        "genome_chrom_size": Path(cache_config.references.genome_chrom_size.file_path),
        "reference_genome": Path(cache_config.references.reference_genome.file_path),
        "refgene_bed": refgene_bed_file,
        "refgene_flat": refgene_flat_file,
        "refgene_txt": Path(cache_config.references.refgene_txt.file_path),
    }


@pytest.fixture(scope="session", name="cadd_snv_indexed_file")
def fixture_cadd_snv_indexed_file(session_tmp_path: Path) -> Path:
    """Return dummy cadd snv indexed file."""
    reference_file: Path = Path(
        session_tmp_path, "variants", "hg19.cadd_snv.tsv.gz.tbi"
    )
    reference_file.touch()
    return reference_file


@pytest.fixture(scope="session", name="delly_exclusion_converted_file")
def fixture_delly_exclusion_converted_file(session_tmp_path: Path) -> Path:
    """Return dummy delly exclusion converted file."""
    reference_file: Path = Path(
        session_tmp_path, "genome", "delly_exclusion_converted.tsv"
    )
    reference_file.touch()
    return reference_file


@pytest.fixture(scope="session", name="clinvar_file")
def fixture_clinvar_file(session_tmp_path: Path) -> Path:
    """Return dummy clinvar file."""
    clinvar_file: Path = Path(session_tmp_path, "variants", "clinvar.vcf.gz")
    clinvar_file.touch()
    return clinvar_file


@pytest.fixture(scope="session", name="cosmic_file")
def fixture_cosmic_file(session_tmp_path: Path) -> Path:
    """Return dummy cosmic file."""
    cosmic_file: Path = Path(
        session_tmp_path, "variants", "cosmic_coding_muts_v97.vcf.gz"
    )
    cosmic_file.touch()
    return cosmic_file


@pytest.fixture(scope="session", name="dbsnp_file")
def fixture_dbsnp_file(session_tmp_path: Path) -> Path:
    """Return dummy dbsnp file."""
    dbsnp_file: Path = Path(session_tmp_path, "variants", "dbsnp_grch37_b138.vcf.gz")
    dbsnp_file.touch()
    return dbsnp_file


@pytest.fixture(scope="session", name="hc_vcf_1kg_file")
def fixture_hc_vcf_1kg(session_tmp_path: Path) -> Path:
    """Return dummy high confidence 1000 genome vcf file."""
    hc_vcf_1kg_file: Path = Path(
        session_tmp_path, "variants", "1kg_phase1_snps_high_confidence_b37.vcf.gz"
    )
    hc_vcf_1kg_file.touch()
    return hc_vcf_1kg_file


@pytest.fixture(scope="session", name="known_indel_1kg_file")
def fixture_known_indel_1kg_file(session_tmp_path: Path) -> Path:
    """Return dummy 1000 genome known indels vcf file."""
    known_indel_1kg_file: Path = Path(
        session_tmp_path, "variants", "1kg_known_indels_b37.vcf.gz"
    )
    known_indel_1kg_file.touch()
    return known_indel_1kg_file


@pytest.fixture(scope="session", name="mills_1kg_file")
def fixture_mills_1kg_file(session_tmp_path: Path) -> Path:
    """Return dummy Mills' high confidence indels vcf file."""
    mills_1kg_file: Path = Path(session_tmp_path, "variants", "mills_1kg_index.vcf.gz")
    mills_1kg_file.touch()
    return mills_1kg_file


@pytest.fixture(scope="session", name="somalier_sites_file")
def fixture_somalier_sites_file(session_tmp_path: Path) -> Path:
    """Return dummy somalier sites vcf file."""
    somalier_sites_file: Path = Path(
        session_tmp_path, "variants", "GRCh37.somalier.sites.vcf.gz"
    )
    somalier_sites_file.touch()
    return somalier_sites_file


@pytest.fixture(scope="session", name="vcf_1kg_file")
def fixture_vcf_1kg_file(session_tmp_path: Path) -> Path:
    """Return dummy 1000 genome all snps file."""
    vcf_1kg_file: Path = Path(
        session_tmp_path, "variants", "1k_genome_wgs_p1_v3_all_sites.vcf.gz"
    )
    vcf_1kg_file.touch()
    return vcf_1kg_file


@pytest.fixture(scope="session", name="analysis_references_hg_data")
def fixture_analysis_references_hg_data(
    cache_config: CacheConfig,
    analysis_references_data: Dict[str, Path],
    delly_exclusion_converted_file: Path,
    clinvar_file: Path,
    cosmic_file: Path,
    dbsnp_file: Path,
    hc_vcf_1kg_file: Path,
    known_indel_1kg_file: Path,
    mills_1kg_file: Path,
    somalier_sites_file: Path,
    vcf_1kg_file: Path,
) -> Dict[str, Path]:
    """Return human genome analysis references model data."""
    analysis_references_hg_data: Dict[str, Path] = {
        "access_regions": Path(cache_config.references.access_regions.file_path),
        "ascat_chr_y_loci": Path(cache_config.references.ascat_chr_y_loci.file_path),
        "ascat_gc_correction": Path(
            cache_config.references.ascat_gc_correction.file_path
        ),
        "cadd_snv": Path(cache_config.references.cadd_snv.file_path),
        "simple_repeat": Path(cache_config.references.simple_repeat.file_path),
        "clinvar": clinvar_file,
        "cosmic": cosmic_file,
        "dbsnp": dbsnp_file,
        "delly_exclusion": Path(cache_config.references.delly_exclusion.file_path),
        "delly_exclusion_converted": delly_exclusion_converted_file,
        "delly_mappability": Path(cache_config.references.delly_mappability.file_path),
        "gnomad_variant": Path(cache_config.references.gnomad_variant.file_path),
        "hc_vcf_1kg": hc_vcf_1kg_file,
        "known_indel_1kg": known_indel_1kg_file,
        "mills_1kg": mills_1kg_file,
        "rank_score": Path(cache_config.references.rank_score.file_path),
        "somalier_sites": somalier_sites_file,
        "vcf_1kg": vcf_1kg_file,
        "vep_dir": cache_config.references_dir,
        "wgs_calling_regions": Path(
            cache_config.references.wgs_calling_regions.file_path
        ),
    }
    analysis_references_hg_data.update(analysis_references_data)
    return analysis_references_hg_data


@pytest.fixture(scope="session", name="analysis_references_hg")
def fixture_analysis_references_hg(
    analysis_references_hg_data: Dict[str, Path]
) -> AnalysisReferencesHg:
    """Return mocked human genome analysis references model."""
    return AnalysisReferencesHg(**analysis_references_hg_data)


@pytest.fixture(scope="session", name="reference_url")
def fixture_reference_url() -> Url:
    """Return dummy reference url."""
    return Url("gs://gatk-legacy-bundles/b37/reference.vcf.gz")


@pytest.fixture(scope="session", name="reference_file")
def fixture_reference_file(session_tmp_path: Path) -> Path:
    """Return dummy reference file."""
    reference_file: Path = Path(session_tmp_path, "reference.vcf")
    reference_file.touch()
    return reference_file


@pytest.fixture(scope="session", name="reference_url_data")
def fixture_reference_url_data(
    reference_url: Url, reference_file: Path, cosmic_key: str
) -> Dict[str, Any]:
    """return reference url model data."""
    return {
        "url": reference_url,
        "file_type": FileType.VCF,
        "gzip": False,
        "file_name": "reference.vcf",
        "dir_name": "variants",
        "file_path": reference_file.as_posix(),
        "secret": cosmic_key,
    }


@pytest.fixture(scope="session", name="references_data")
def fixture_references_data(
    cache_config: CacheConfig,
) -> Dict[str, dict]:
    """Return references model data."""
    return {
        "genome_chrom_size": cache_config.references.genome_chrom_size,
        "reference_genome": cache_config.references.reference_genome,
        "refgene_sql": cache_config.references.refgene_sql,
        "refgene_txt": cache_config.references.refgene_txt,
    }


@pytest.fixture(scope="session", name="references")
def fixture_references(references_data: Dict[str, dict]) -> References:
    """Return mocked references model."""
    return References(**references_data)


@pytest.fixture(scope="session", name="references_hg_data")
def fixture_references_hg_data(
    cache_config: CacheConfig,
) -> Dict[str, dict]:
    """Return human genome references model data."""
    return dict(cache_config.references)


@pytest.fixture(scope="session", name="references_hg")
def fixture_references_hg(references_hg_data: Dict[str, dict]) -> ReferencesHg:
    """Return mocked human genome references model."""
    return ReferencesHg(**references_hg_data)


@pytest.fixture(scope="session", name="singularity_bind_path_data")
def fixture_singularity_bind_path_data(session_tmp_path: Path) -> Dict[str, Path]:
    """Return singularity bind path data."""
    return {"source": session_tmp_path, "destination": Path("/")}


@pytest.fixture(scope="session", name="singularity_bind_path")
def fixture_singularity_bind_path(
    singularity_bind_path_data: Dict[str, Path]
) -> SingularityBindPath:
    """Return mocked singularity bind path model."""
    return SingularityBindPath(**singularity_bind_path_data)


@pytest.fixture(scope="session", name="snakemake_options_command")
def fixture_snakemake_options_command() -> List[str]:
    """Return mocked singularity bind path model."""
    return ["--cores", "36"]


@pytest.fixture(scope="session", name="mail_user_option")
def fixture_mail_user_option() -> str:
    """Return mail user option."""
    return "balsamic@scilifelab.se"


@pytest.fixture(scope="session", name="snakemake_executable_data")
def fixture_snakemake_executable_data(
    case_id_tumor_only: str,
    reference_file: Path,
    session_tmp_path: Path,
    mail_user_option: str,
    singularity_bind_path: SingularityBindPath,
    snakemake_options_command: List[str],
) -> Dict[str, Any]:
    """Return snakemake executable model data."""
    return {
        "account": ClusterAccount.DEVELOPMENT.value,
        "case_id": case_id_tumor_only,
        "cluster_config_path": reference_file,
        "config_path": reference_file,
        "disable_variant_caller": "tnscope,vardict",
        "log_dir": session_tmp_path,
        "mail_user": mail_user_option,
        "profile": ClusterProfile.SLURM,
        "qos": QOS.HIGH,
        "quiet": True,
        "run_analysis": True,
        "run_mode": RunMode.CLUSTER,
        "script_dir": session_tmp_path,
        "singularity_bind_paths": [singularity_bind_path],
        "snakefile": reference_file,
        "snakemake_options": snakemake_options_command,
        "working_dir": session_tmp_path,
    }


@pytest.fixture(scope="session", name="snakemake_executable")
def fixture_snakemake_executable(
    snakemake_executable_data: Dict[str, Any]
) -> SnakemakeExecutable:
    """Return mocked snakemake executable model."""
    return SnakemakeExecutable(**snakemake_executable_data)


@pytest.fixture(scope="session", name="snakemake_executable_validated_data")
def fixture_snakemake_executable_validated_data(
    case_id_tumor_only: str,
    reference_file: Path,
    session_tmp_path: Path,
    mail_user_option: str,
    singularity_bind_path: SingularityBindPath,
    snakemake_options_command: List[str],
) -> Dict[str, Any]:
    """Return snakemake model expected data."""
    return {
        "account": ClusterAccount.DEVELOPMENT.value,
        "benchmark": False,
        "case_id": case_id_tumor_only,
        "cluster_config_path": reference_file,
        "config_path": reference_file,
        "disable_variant_caller": "disable_variant_caller=tnscope,vardict",
        "dragen": False,
        "force": False,
        "log_dir": session_tmp_path,
        "mail_type": None,
        "mail_user": mail_user_option,
        "profile": ClusterProfile.SLURM,
        "qos": QOS.HIGH,
        "quiet": True,
        "report_path": None,
        "run_analysis": True,
        "run_mode": RunMode.CLUSTER,
        "script_dir": session_tmp_path,
        "singularity_bind_paths": [singularity_bind_path],
        "snakefile": reference_file,
        "snakemake_options": snakemake_options_command,
        "working_dir": session_tmp_path,
    }


@pytest.fixture(scope="session")
def job_id() -> str:
    """Return cluster job identifier."""
    return "12345"


@pytest.fixture(scope="session")
def job_properties() -> Dict[str, Any]:
    """Cluster job properties."""
    return {
        "cluster": {
            "partition": "core",
            "n": "1",
            "time": "10:00:00",
            "mail_type": ClusterMailType.ALL.value,
        }
    }


@pytest.fixture(scope="session")
def scheduler_data(
    case_id_tumor_only: str,
    job_properties: Dict[str, Any],
    empty_file: Path,
    empty_dir: Path,
    mail_user_option: str,
) -> Dict[str, Any]:
    """Return raw scheduler model data."""
    return {
        "account": ClusterAccount.DEVELOPMENT.value,
        "case_id": case_id_tumor_only,
        "dependencies": ["1", "2", "3"],
        "job_properties": job_properties,
        "job_script": empty_file.as_posix(),
        "log_dir": empty_dir.as_posix(),
        "mail_user": mail_user_option,
        "mail_type": ClusterMailType.FAIL.value,
        "profile": ClusterProfile.SLURM.value,
        "qos": QOS.HIGH,
    }


@pytest.fixture(scope="session")
def scheduler_validated_data(
    case_id_tumor_only: str,
    job_properties: Dict[str, Any],
    empty_file: Path,
    empty_dir: Path,
    mail_user_option: str,
) -> Dict[str, Any]:
    """Return scheduler model validated data."""
    return {
        "account": f"--account {ClusterAccount.DEVELOPMENT}",
        "benchmark": False,
        "case_id": case_id_tumor_only,
        "dependencies": ["1", "2", "3"],
        "job_properties": job_properties,
        "job_script": empty_file,
        "log_dir": empty_dir,
        "mail_type": "--mail-type FAIL",
        "mail_user": f"--mail-user {mail_user_option}",
        "profile": ClusterProfile.SLURM,
        "profiling_interval": 10,
        "profiling_type": "task",
        "qos": "--qos high",
    }


@pytest.fixture(scope="session")
def scheduler_model(scheduler_data: Dict[str, Any]) -> Scheduler:
    """Return scheduler pydantic model."""
    return Scheduler(**scheduler_data)
