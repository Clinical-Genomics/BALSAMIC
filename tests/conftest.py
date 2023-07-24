import json
import os
from datetime import datetime
from distutils.dir_util import copy_tree
from functools import partial
from pathlib import Path
from typing import Dict, Any
from unittest import mock

import pytest
from _pytest.tmpdir import TempPathFactory
from click.testing import CliRunner

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.commands.base import cli
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV, RunMode
from BALSAMIC.constants.cache import (
    DockerContainers,
    GenomeVersion,
    REFERENCE_FILES,
)
from BALSAMIC.constants.cluster import (
    ClusterConfigType,
    ClusterAccount,
    QOS,
    ClusterMailType,
    ClusterProfile,
)
from BALSAMIC.constants.constants import FileType
from BALSAMIC.constants.paths import CONSTANTS_DIR
from BALSAMIC.constants.workflow_params import VCF_DICT
from BALSAMIC.models.cache import (
    CacheAnalysis,
    CacheConfig,
    References,
    ReferencesHg,
    AnalysisReferencesHg,
)
from BALSAMIC.models.snakemake import SingularityBindPath, Snakemake
from BALSAMIC.utils.io import read_json, read_yaml
from .helpers import ConfigHelper, Map

MOCKED_OS_ENVIRON = "os.environ"


@pytest.fixture(scope="session", name="session_tmp_path")
def fixture_session_tmp_path(tmp_path_factory: TempPathFactory) -> Path:
    """Return a non-existent files directory path."""
    return tmp_path_factory.mktemp("session_tests")


@pytest.fixture(scope="session")
def tumor_sample_name() -> str:
    """
    Creates mock name for tumor sample.
    """
    return "ACC1"


@pytest.fixture(scope="session")
def normal_sample_name() -> str:
    """
    Creates mock name for normal sample.
    """
    return "ACC2"


@pytest.fixture(scope="session")
def case_id_tumor_only() -> str:
    """
    creates mock case-id for TGA tumor-only.
    """
    return "sample_tumor_only"


@pytest.fixture(scope="session")
def case_id_tumor_only_pon() -> str:
    """
    Creates mock case-id for TGA PON tumor-only.
    """
    return "sample_tumor_only_pon"


@pytest.fixture(scope="session")
def case_id_tumor_only_umi() -> str:
    """
    Creates mock case-id for TGA PON tumor-only.
    """
    return "sample_tumor_only_umi"


@pytest.fixture(scope="session")
def case_id_tumor_normal() -> str:
    """
    Creates mock case-id for TGA tumor-normal.
    """
    return "sample_tumor_normal"


@pytest.fixture(scope="session")
def case_id_tumor_only_wgs() -> str:
    """
    Creates mock case-id for WGS tumor-only.
    """
    return "sample_tumor_only_wgs"


@pytest.fixture(scope="session")
def case_id_tumor_normal_wgs() -> str:
    """
    Creates mock case-id for WGS tumor-normal.
    """
    return "sample_tumor_normal_wgs"


@pytest.fixture(scope="session")
def fastq_dir(case_id_tumor_only: str, analysis_dir: str):
    """
    Creates path for mock FastQ directory.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    Path(fastq_dir, "ACC1_XXXXX_R_1.fastq.gz").touch()
    Path(fastq_dir, "ACC1_XXXXX_R_2.fastq.gz").touch()
    Path(fastq_dir, "ACC2_XXXXX_R_1.fastq.gz").touch()
    Path(fastq_dir, "ACC2_XXXXX_R_2.fastq.gz").touch()
    return fastq_dir.as_posix()


@pytest.fixture
def cli_runner():
    """
    Runs click for command line interface testing.
    """
    runner = CliRunner()
    return runner


@pytest.fixture
def invoke_cli(cli_runner):
    """
    Invokes cli commands with options.
    """
    return partial(cli_runner.invoke, cli)


@pytest.fixture(scope="session")
def environ():
    """
    Creates operating system's environment object.
    """
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
        "reference_genome": "tests/test_data/references/genome/human_g1k_v37_decoy.fasta",
        "dbsnp": "tests/test_data/references/variants/dbsnp_grch37_b138.vcf.gz",
        "vcf_1kg": "tests/test_data/references/variants/1k_genome_wgs_p1_v3_all_sites.vcf.gz",
        "hc_vcf_1kg": "tests/test_data/references/variants/1kg_phase1_snps_high_confidence_b37.vcf.gz",
        "known_indel_1kg": "tests/test_data/references/variants/1kg_known_indels_b37.vcf.gz",
        "mills_1kg": "tests/test_data/references/variants/mills_1kg_index.vcf.gz",
        "gnomad_variant": "tests/test_data/reference/variants/gnomad.genomes.r2.1.1.sites.vcf.bgz",
        "cosmic": "tests/test_data/references/variants/cosmic_coding_muts_v89.vcf.gz",
        "vep_dir": "tests/test_data/references/vep/",
        "refgene_flat": "tests/test_data/references/genome/refseq.flat",
        "refgene_txt": "tests/test_data/references/genome/refGene.txt",
        "wgs_calling_regions": "tests/test_data/references/genome/wgs_calling_regions.v1",
        "genome_chrom_size": "tests/test_data/references/genome/hg19.chrom.sizes",
        "refgene_bed": "tests/test_data/references/genome/refseq.flat.bed",
        "rank_score": "tests/test_data/references/genome/cancer_rank_model_-v0.1-.ini",
        "access_regions": "tests/test_data/references/genome/access-5k-mappable.hg19.bed",
        "delly_exclusion": "tests/test_data/references/genome/delly_exclusion.tsv",
        "delly_exclusion_converted": "tests/test_data/references/genome/delly_exclusion_converted.tsv",
        "delly_mappability": "tests/test_data/references/genome/delly_mappability.gz",
        "delly_mappability_gindex": "tests/test_data/references/genome/delly_mappability.gz.gzi",
        "delly_mappability_findex": "tests/test_data/references/genome/delly_mappability.fai",
        "ascat_gc_correction": "tests/test_data/references/genome/GRCh37_SnpGcCorrections.tsv",
        "ascat_chr_y_loci": "tests/test_data/references/genome/GRCh37_Y.loci",
        "clinvar": "tests/test_data/references/genome/clinvar.vcf.gz",
        "somalier_sites": "tests/test_data/references/variants/GRCh37.somalier.sites.vcf.gz",
        "cadd_snv": "tests/test_data/references/variants/hg19.cadd_snv.tsv.gz",
    }


@pytest.fixture(scope="session")
def test_data_dir() -> str:
    """
    Creates path for test data directory.
    """
    return "tests/test_data"


@pytest.fixture(scope="session")
def pon_fastq_path(test_data_dir: str) -> str:
    """
    Creates path for FASTQ directory for Panel Of Normal (PON).
    """
    return Path(test_data_dir, "fastq").as_posix()


@pytest.fixture(scope="session")
def reference_panel_dir_path(test_data_dir: str) -> str:
    """
    Creates path for reference panel directory.
    """
    return Path(test_data_dir, "references", "panel").as_posix()


@pytest.fixture(scope="session")
def reference_variants_dir_path(test_data_dir: str) -> str:
    """
    Created path for reference variants directory.
    """
    return Path(test_data_dir, "references", "variants").as_posix()


@pytest.fixture(scope="session")
def config_path(test_data_dir: str) -> str:
    """
    Created path for config json file.
    """
    return Path(test_data_dir, "config.json").as_posix()


@pytest.fixture(scope="session")
def config_dict(config_path: str) -> str:
    """
    Reads and returns config from json.
    """
    return read_json(config_path)


@pytest.fixture(scope="session")
def panel_bed_file(reference_panel_dir_path: str) -> str:
    """
    Created path for panel bed file.
    """
    return Path(reference_panel_dir_path, "panel.bed").as_posix()


@pytest.fixture(scope="session")
def background_variant_file(reference_panel_dir_path: str) -> str:
    """
    Created path for background variants for TGA.
    """
    return Path(reference_panel_dir_path, "background_variants.txt").as_posix()


@pytest.fixture(scope="session")
def pon_cnn(reference_panel_dir_path: str) -> str:
    """
    Creates path for Panel Of Normal (PON), cnn file for cnvkit.
    """
    return Path(reference_panel_dir_path, "test_panel_ponn.cnn").as_posix()


@pytest.fixture(scope="session")
def clinical_snv_observations(reference_variants_dir_path: str) -> str:
    """
    Creates path for clinical SNVs from loqusDB.
    """
    return Path(reference_variants_dir_path, "clinical_snv_variants.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def cancer_germline_snv_observations(reference_variants_dir_path: str) -> str:
    """
    Creates path of cancer germline SNVs from loqusDB.
    """
    return Path(
        reference_variants_dir_path, "cancer_germline_snv_variants.vcf.gz"
    ).as_posix()


@pytest.fixture(scope="session")
def cancer_somatic_snv_observations(reference_variants_dir_path: str) -> str:
    """
    Creates path for somatic SNVs from loqusDB.
    """
    return Path(
        reference_variants_dir_path, "cancer_somatic_snv_variants.vcf.gz"
    ).as_posix()


@pytest.fixture(scope="session")
def clinical_sv_observations(reference_variants_dir_path: str) -> str:
    """
    Creates path for clinical SVs from loqusDB.
    """
    return Path(reference_variants_dir_path, "clinical_sv_variants.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def swegen_snv_frequency(reference_variants_dir_path: str) -> str:
    """
    Creates path for Swegen SNVs.
    """
    return Path(reference_variants_dir_path, "swegen_snv.vcf.gz").as_posix()


@pytest.fixture(scope="session")
def swegen_sv_frequency(reference_variants_dir_path: str) -> str:
    """
    Creates path for Swegen SVs.
    """
    return Path(reference_variants_dir_path, "swegen_sv.vcf.gz").as_posix()


@pytest.fixture(scope="session", name="invalid_json_file")
def fixture_invalid_json_file(session_tmp_path: Path) -> Path:
    """Return a non-existent json file path."""
    return Path(session_tmp_path, "invalid_file.json")


@pytest.fixture(scope="session", name="json_file")
def fixture_json_file(session_tmp_path: Path) -> Path:
    """Return a mocked json file path."""
    return Path(session_tmp_path, "write_json.json")


@pytest.fixture(scope="session", name="config_json")
def fixture_config_json() -> str:
    """Return Balsamic analysis config json file name."""
    return "config.json"


@pytest.fixture(scope="session", name="reference_graph")
def fixture_reference_graph() -> str:
    """Return Balsamic reference graph pdf file name."""
    return "reference_graph.pdf"


@pytest.fixture(scope="session")
def sentieon_license(tmp_path_factory):
    """
    Creates Sentieon's license path
    """
    sentieon_license_dir = tmp_path_factory.mktemp("sentieon_licence")
    sentieon_license_path = sentieon_license_dir / "license_file.lic"
    sentieon_license_path.touch()

    return sentieon_license_path.as_posix()


@pytest.fixture(scope="session")
def sentieon_install_dir(tmp_path_factory):
    """
    Creates install directory for Sentieon tools
    """
    sentieon_install_dir = tmp_path_factory.mktemp("sentieon_install_dir")
    Path(sentieon_install_dir / "bin").mkdir(exist_ok=True)
    sentieon_executable = sentieon_install_dir / "bin" / "sentieon"
    sentieon_executable.touch()

    return sentieon_install_dir.as_posix()


@pytest.fixture()
def no_write_perm_path(tmp_path_factory) -> str:
    """
    A path with no write permissions.
    """
    bad_perm_path: Path = tmp_path_factory.mktemp("bad_perm_path")
    bad_perm_path.chmod(0o444)
    return bad_perm_path.as_posix()


@pytest.fixture(scope="session")
def balsamic_cache(tmp_path_factory, reference):
    """
    Creates balsamic-cache and returns the path for balsamic-cache
    """

    cache_dir = tmp_path_factory.mktemp("balsmic_cache")

    cache_container = cache_dir / balsamic_version / "containers" / "align_qc"
    cache_container.mkdir(parents=True, exist_ok=True)
    cache_container_example = cache_container / "example.sif"
    cache_container_example.touch()

    cache_reference = cache_dir / balsamic_version / "hg19"
    cache_reference.mkdir(parents=True, exist_ok=True)

    cache_reference_json = cache_reference / "reference.json"
    cache_reference_json.touch()
    with open(cache_reference_json, "w") as fp:
        json.dump(reference, fp)

    return cache_dir.as_posix()


@pytest.fixture(scope="session")
def analysis_dir(tmp_path_factory: TempPathFactory) -> str:
    """
    Creates and returns the directory where the case analysis will be saved.
    """
    analysis_dir = tmp_path_factory.mktemp("analysis", numbered=False)
    return analysis_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only(analysis_dir: str, case_id_tumor_only: str) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-only.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Fill the concat fastq path folder with the concatenated fastq files
    concat_dir = Path(analysis_dir, case_id_tumor_only, "analysis", "concat")
    concat_dir.mkdir(parents=True, exist_ok=True)
    for read in [1, 2]:
        Path(concat_dir, f"ACC1_R_{read}.fastq.gz").touch()

    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_pon(analysis_dir: str, case_id_tumor_only_pon: str) -> str:
    """
    Creates and returns the directory containing the FASTQs for PON, tumor-only.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_pon, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_umi(analysis_dir: str, case_id_tumor_only_umi: str) -> str:
    """
    Creates and returns the directory containing the FASTQs for UMI, tumor only.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_umi, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_normal(analysis_dir: str, case_id_tumor_normal: str) -> str:
    """
    Creates and returns the directory containing the FASTQs for tumor-normal.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_wgs(analysis_dir: str, case_id_tumor_only_wgs: str) -> str:
    """
    Creates and returns the directory containing the FASTQs for WGS tumor-only.
    """
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_wgs, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_normal_wgs(analysis_dir: str, case_id_tumor_normal_wgs: str) -> str:
    """
    Creates and returns the directory containing the FASTQs for WGS tumor-normal.
    """

    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal_wgs, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def snakemake_job_script(tmp_path_factory, tumor_normal_config):
    """
    Creates a dummy snakemake jobscript
    """

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


@pytest.fixture(scope="session")
def tumor_normal_config(
    case_id_tumor_normal: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal: str,
    balsamic_cache: str,
    background_variant_file: str,
    panel_bed_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    swegen_snv_frequency: str,
    swegen_sv_frequency: str,
    clinical_snv_observations: str,
    clinical_sv_observations: str,
    cancer_germline_snv_observations: str,
    cancer_somatic_snv_observations: str,
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
                case_id_tumor_normal,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_normal,
                "--background-variants",
                background_variant_file,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
                "--swegen-snv",
                swegen_snv_frequency,
                "--swegen-sv",
                swegen_sv_frequency,
                "--clinical-snv-observations",
                clinical_snv_observations,
                "--clinical-sv-observations",
                clinical_sv_observations,
                "--cancer-germline-snv-observations",
                cancer_germline_snv_observations,
                "--cancer-somatic-snv-observations",
                cancer_somatic_snv_observations,
            ],
        )

    qc_dir = Path(analysis_dir, case_id_tumor_normal, "analysis", "qc")
    qc_dir.mkdir(parents=True, exist_ok=False)
    copy_tree("tests/test_data/qc_files/analysis/qc/", qc_dir.as_posix())

    return Path(
        analysis_dir, case_id_tumor_normal, f"{case_id_tumor_normal}.{FileType.JSON}"
    ).as_posix()


@pytest.fixture(name="helpers")
def fixture_config_helpers():
    """Helper fixture for case config files"""
    return ConfigHelper()


@pytest.fixture(scope="session")
def tumor_normal_wgs_config(
    case_id_tumor_normal_wgs: str,
    tumor_sample_name: str,
    normal_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_normal_wgs: str,
    balsamic_cache: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    swegen_snv_frequency: str,
    swegen_sv_frequency: str,
    clinical_snv_observations: str,
    clinical_sv_observations: str,
    cancer_germline_snv_observations: str,
    cancer_somatic_snv_observations: str,
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
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
                "--normal-sample-name",
                normal_sample_name,
                "--swegen-snv",
                swegen_snv_frequency,
                "--swegen-sv",
                swegen_sv_frequency,
                "--clinical-snv-observations",
                clinical_snv_observations,
                "--clinical-sv-observations",
                clinical_sv_observations,
                "--cancer-germline-snv-observations",
                cancer_germline_snv_observations,
                "--cancer-somatic-snv-observations",
                cancer_somatic_snv_observations,
            ],
        )

    return Path(
        analysis_dir,
        case_id_tumor_normal_wgs,
        f"{case_id_tumor_normal_wgs}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_config(
    case_id_tumor_only: str,
    tumor_sample_name: str,
    balsamic_cache: str,
    analysis_dir: str,
    fastq_dir_tumor_only: str,
    panel_bed_file: str,
    background_variant_file: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    swegen_snv_frequency: str,
    swegen_sv_frequency: str,
    clinical_snv_observations: str,
    clinical_sv_observations: str,
    cancer_germline_snv_observations: str,
    cancer_somatic_snv_observations: str,
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
                "-p",
                panel_bed_file,
                "--balsamic-cache",
                balsamic_cache,
                "--background-variants",
                background_variant_file,
                "--tumor-sample-name",
                tumor_sample_name,
                "--swegen-snv",
                swegen_snv_frequency,
                "--swegen-sv",
                swegen_sv_frequency,
                "--clinical-snv-observations",
                clinical_snv_observations,
                "--clinical-sv-observations",
                clinical_sv_observations,
                "--cancer-germline-snv-observations",
                cancer_germline_snv_observations,
                "--cancer-somatic-snv-observations",
                cancer_somatic_snv_observations,
            ],
        )

    qc_dir = Path(analysis_dir, case_id_tumor_only, "analysis", "qc")
    qc_dir.mkdir(parents=True, exist_ok=False)
    copy_tree("tests/test_data/qc_files/analysis/qc/", qc_dir.as_posix())

    return Path(
        analysis_dir, case_id_tumor_only, f"{case_id_tumor_only}.{FileType.JSON}"
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_wgs_config(
    case_id_tumor_only_wgs: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only_wgs: str,
    balsamic_cache: str,
    sentieon_license: str,
    sentieon_install_dir: str,
    swegen_snv_frequency: str,
    swegen_sv_frequency: str,
    clinical_snv_observations: str,
    clinical_sv_observations: str,
    cancer_germline_snv_observations: str,
    cancer_somatic_snv_observations: str,
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
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
                "--swegen-snv",
                swegen_snv_frequency,
                "--swegen-sv",
                swegen_sv_frequency,
                "--clinical-snv-observations",
                clinical_snv_observations,
                "--clinical-sv-observations",
                clinical_sv_observations,
                "--cancer-germline-snv-observations",
                cancer_germline_snv_observations,
                "--cancer-somatic-snv-observations",
                cancer_somatic_snv_observations,
            ],
        )

    return Path(
        analysis_dir,
        case_id_tumor_only_wgs,
        f"{case_id_tumor_only_wgs}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def tumor_only_pon_config(
    case_id_tumor_only_pon: str,
    tumor_sample_name: str,
    analysis_dir: str,
    fastq_dir_tumor_only_pon: str,
    panel_bed_file: str,
    pon_cnn: str,
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
                case_id_tumor_only_pon,
                "--analysis-dir",
                analysis_dir,
                "--fastq-path",
                fastq_dir_tumor_only_pon,
                "-p",
                panel_bed_file,
                "--pon-cnn",
                pon_cnn,
                "--balsamic-cache",
                balsamic_cache,
                "--tumor-sample-name",
                tumor_sample_name,
            ],
        )

    return Path(
        analysis_dir,
        case_id_tumor_only_pon,
        f"{case_id_tumor_only_pon}.{FileType.JSON}",
    ).as_posix()


@pytest.fixture(scope="session")
def sample_config(tumor_sample_name: str, normal_sample_name: str):
    """
    sample config dict to test workflow utils
    """
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
        "vcf": VCF_DICT,
        "samples": {
            tumor_sample_name: {"type": "tumor"},
            normal_sample_name: {"type": "normal"},
        },
        "umiworkflow": "true",
    }

    return sample_config


@pytest.fixture(scope="session")
def analysis_path():
    """Analysis test path"""
    return "tests/test_data/qc_files/analysis"


@pytest.fixture(scope="session")
def multiqc_data_path(analysis_path):
    """multiqc_data.json test path"""
    return os.path.join(analysis_path, "qc", "multiqc_data", "multiqc_data.json")


@pytest.fixture(scope="session")
def multiqc_data_dict(multiqc_data_path):
    """multiqc_data.json test path"""
    return read_json(multiqc_data_path)


@pytest.fixture(scope="session")
def metrics_yaml_path(analysis_path):
    """sample_tumor_only_metrics_deliverables.yaml test path"""
    return os.path.join(
        analysis_path, "qc", "sample_tumor_only_metrics_deliverables.yaml"
    )


@pytest.fixture(scope="session")
def bcftools_counts_path(analysis_path):
    """svdb.all.filtered.pass.stats test path"""
    return os.path.join(
        analysis_path, "vep", "SNV.somatic.case.svdb.all.filtered.pass.stats"
    )


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
def snakemake_fastqc_rule(tumor_only_config, helpers):
    """FastQC snakemake mock rule"""

    helpers.read_config(tumor_only_config)
    fastq_path = os.path.join(
        helpers.analysis_dir,
        helpers.case_id,
        "analysis",
        "concat",
        "ACC1_R_{read}.fastq.gz",
    )

    return Map(
        {
            "fastqc": Map(
                {
                    "params": Map(
                        {
                            "housekeeper_id": {
                                "id": "sample_tumor_only",
                                "tags": "quality-trimmed-seq",
                            }
                        }
                    ),
                    "output": Map(
                        {
                            "_names": Map({"fastqc": fastq_path}),
                            "fastqc": fastq_path,
                        }
                    ),
                    "rule": Map(
                        {
                            "name": "fastq",
                            "output": [
                                fastq_path,
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
    """Return a time stamp of today's date in date time format."""
    return datetime.now()


@pytest.fixture(scope="session", name="cosmic_key")
def fixture_cosmic_key() -> str:
    """Return mocked COSMIC key."""
    return "ZW1haWxAZXhhbXBsZS5jb206bXljb3NtaWNwYXNzd29yZAo="


@pytest.fixture(scope="session", name="develop_containers")
def fixture_develop_containers() -> Dict[str, str]:
    """Return develop docker hub containers dictionary."""
    return {
        DockerContainers.ASCAT.value: "docker://clinicalgenomics/balsamic:develop-ascatNgs",
        DockerContainers.VCF2CYTOSURE.value: "docker://clinicalgenomics/balsamic:develop-vcf2cytosure",
        DockerContainers.PYTHON_3.value: "docker://clinicalgenomics/balsamic:develop-varcall_py3",
        DockerContainers.BALSAMIC.value: "docker://clinicalgenomics/balsamic:develop-balsamic",
        DockerContainers.SOMALIER.value: "docker://clinicalgenomics/balsamic:develop-somalier",
        DockerContainers.CNVPYTOR.value: "docker://clinicalgenomics/balsamic:develop-cnvpytor",
        DockerContainers.ALIGN_QC.value: "docker://clinicalgenomics/balsamic:develop-align_qc",
        DockerContainers.ANNOTATE.value: "docker://clinicalgenomics/balsamic:develop-annotate",
        DockerContainers.PYTHON_27.value: "docker://clinicalgenomics/balsamic:develop-varcall_py27",
        DockerContainers.CNVKIT.value: "docker://clinicalgenomics/balsamic:develop-varcall_cnvkit",
        DockerContainers.COVERAGE_QC.value: "docker://clinicalgenomics/balsamic:develop-coverage_qc",
        DockerContainers.DELLY.value: "docker://clinicalgenomics/balsamic:develop-delly",
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
) -> CacheAnalysis:
    """Return mocked human genome analysis references model."""
    return AnalysisReferencesHg(**analysis_references_hg_data)


@pytest.fixture(scope="session", name="reference_url")
def fixture_reference_url() -> str:
    """Return dummy reference url."""
    return "gs://gatk-legacy-bundles/b37/reference.vcf.gz"


@pytest.fixture(scope="session", name="reference_file")
def fixture_reference_file(session_tmp_path: Path) -> Path:
    """Return dummy reference file."""
    reference_file: Path = Path(session_tmp_path, "reference.vcf")
    reference_file.touch()
    return reference_file


@pytest.fixture(scope="session", name="reference_url_data")
def fixture_reference_url_data(
    reference_url: str, reference_file: Path, cosmic_key: str
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
def fixture_singularity_bind_path_data(reference_file: Path) -> Dict[str, Path]:
    """Return singularity bind path data."""
    return {"source": reference_file, "destination": Path("/")}


@pytest.fixture(scope="session", name="singularity_bind_path")
def fixture_singularity_bind_path(
    singularity_bind_path_data: Dict[str, Path]
) -> SingularityBindPath:
    """Return mocked singularity bind path model."""
    return SingularityBindPath(**singularity_bind_path_data)


@pytest.fixture(scope="session", name="snakemake_data")
def fixture_snakemake_data(
    case_id_tumor_only: str,
    reference_file: Path,
    session_tmp_path: Path,
    singularity_bind_path: SingularityBindPath,
) -> Dict[str, Any]:
    """Return snakemake model data."""
    return {
        "account": ClusterAccount.DEVELOPMENT.value,
        "case_id": case_id_tumor_only,
        "cluster_config_path": reference_file,
        "config_path": reference_file,
        "disable_variant_caller": "tnscope,vardict",
        "log_dir": session_tmp_path,
        "profile": ClusterProfile.SLURM,
        "qos": QOS.HIGH,
        "quiet": True,
        "result_dir": session_tmp_path,
        "run_analysis": True,
        "run_mode": RunMode.CLUSTER,
        "script_dir": session_tmp_path,
        "singularity_bind_paths": [singularity_bind_path],
        "snakefile": reference_file,
        "working_dir": session_tmp_path,
    }


@pytest.fixture(scope="session", name="snakemake")
def fixture_snakemake(snakemake_data: Dict[str, Any]) -> SingularityBindPath:
    """Return mocked singularity model."""
    return Snakemake(**snakemake_data)


@pytest.fixture(scope="session", name="snakemake_validated_data")
def fixture_snakemake_validated_data(
    case_id_tumor_only: str,
    reference_file: Path,
    session_tmp_path: Path,
    singularity_bind_path: SingularityBindPath,
) -> Dict[str, Any]:
    """Return snakemake model expected data."""
    return {
        "account": ClusterAccount.DEVELOPMENT.value,
        "case_id": case_id_tumor_only,
        "cluster_config_path": reference_file,
        "config_path": reference_file,
        "disable_variant_caller": "disable_variant_caller=tnscope,vardict",
        "dragen": False,
        "force": False,
        "log_dir": session_tmp_path,
        "mail_type": None,
        "mail_user": "",
        "profile": ClusterProfile.SLURM,
        "qos": QOS.HIGH,
        "quiet": True,
        "report_path": None,
        "result_dir": session_tmp_path,
        "run_analysis": True,
        "run_mode": RunMode.CLUSTER,
        "script_dir": session_tmp_path,
        "singularity_bind_paths": [singularity_bind_path],
        "slurm_profiler": "",
        "snakefile": reference_file,
        "snakemake_options": None,
        "working_dir": session_tmp_path,
    }
