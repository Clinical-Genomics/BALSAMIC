from datetime import datetime
from typing import Dict, Any

import pytest
import json
import os

from unittest import mock
from distutils.dir_util import copy_tree
from pathlib import Path
from functools import partial


from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV

from BALSAMIC.constants.cache import (
    DockerContainers,
    GenomeVersion,
    REFERENCE_FILES,
    FileType,
)
from _pytest.tmpdir import TempPathFactory

from BALSAMIC.constants.cluster import ClusterConfigType
from BALSAMIC.constants.paths import CONSTANTS_DIR
from BALSAMIC.constants.workflow_params import VCF_DICT
from click.testing import CliRunner

from BALSAMIC.models.cache import (
    CacheAnalysisModel,
    CacheConfigModel,
    ReferencesModel,
    HgReferencesModel,
)
from BALSAMIC.utils.io import read_json, read_yaml
from .helpers import ConfigHelper, Map
from BALSAMIC.commands.base import cli
from BALSAMIC import __version__ as balsamic_version

MOCKED_OS_ENVIRON = "os.environ"


@pytest.fixture(scope="session")
def tumor_sample_name() -> str:
    """Mock tumor sample name."""
    return "ACC1"


@pytest.fixture(scope="session")
def normal_sample_name() -> str:
    """Mock normal sample name."""
    return "ACC2"


@pytest.fixture(scope="session")
def case_id_tumor_only() -> str:
    """Mock TGA tumor-only case ID."""
    return "sample_tumor_only"


@pytest.fixture(scope="session")
def case_id_tumor_only_pon() -> str:
    """Mock TGA PON tumor-only case ID."""
    return "sample_tumor_only_pon"


@pytest.fixture(scope="session")
def case_id_tumor_only_umi() -> str:
    """Mock TGA PON tumor-only case ID."""
    return "sample_tumor_only_umi"


@pytest.fixture(scope="session")
def case_id_tumor_normal() -> str:
    """Mock TGA tumor-normal case ID."""
    return "sample_tumor_normal"


@pytest.fixture(scope="session")
def case_id_tumor_only_wgs() -> str:
    """Mock WGS tumor-only case ID."""
    return "sample_tumor_only_wgs"


@pytest.fixture(scope="session")
def case_id_tumor_normal_wgs() -> str:
    """Mock WGS tumor-normal case ID."""
    return "sample_tumor_normal_wgs"


@pytest.fixture(scope="session")
def fastq_dir(case_id_tumor_only: str, analysis_dir: str):
    """Mock FastQ directory."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    Path(fastq_dir, "ACC1_XXXXX_R_1.fastq.gz").touch()
    Path(fastq_dir, "ACC1_XXXXX_R_2.fastq.gz").touch()
    Path(fastq_dir, "ACC2_XXXXX_R_1.fastq.gz").touch()
    Path(fastq_dir, "ACC2_XXXXX_R_2.fastq.gz").touch()
    return fastq_dir.as_posix()


@pytest.fixture
def cli_runner():
    """click - cli testing"""
    runner = CliRunner()
    return runner


@pytest.fixture
def invoke_cli(cli_runner):
    """invoking cli commands with options"""
    return partial(cli_runner.invoke, cli)


@pytest.fixture(scope="session")
def environ():
    """environment process"""
    return "os.environ"


@pytest.fixture(scope="session")
def cluster_analysis_config_path() -> str:
    """Return cluster analysis configuration file."""
    return Path(CONSTANTS_DIR, ClusterConfigType.ANALYSIS + ".json").as_posix()


@pytest.fixture(scope="session")
def reference():
    """reference json model"""
    return {
        "reference": {
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
            "clinical_snv_observations": "tests/test_data/references/variants/clinical_snv_variants.vcf.gz",
            "cancer_germline_snv_observations": "tests/test_data/references/variants/cancer_germline_snv_variants.vcf.gz",
            "cancer_somatic_snv_observations": "tests/test_data/references/variants/cancer_somatic_snv_variants.vcf.gz",
            "clinical_sv_observations": "tests/test_data/references/variants/clinical_sv_variants.vcf.gz",
            "swegen_snv_frequency": "tests/test_data/references/variants/swegen_snv.vcf.gz",
            "swegen_sv_frequency": "tests/test_data/references/variants/swegen_sv.vcf.gz",
            "somalier_sites": "tests/test_data/references/variants/GRCh37.somalier.sites.vcf.gz",
            "cadd_snv": "tests/test_data/references/variants/hg19.cadd_snv.tsv.gz",
        }
    }


@pytest.fixture(scope="session")
def test_data_dir():
    return "tests/test_data"


@pytest.fixture(scope="session")
def config_path():
    return "tests/test_data/config.json"


@pytest.fixture(scope="session")
def config_dict(config_path):
    return read_json(config_path)


@pytest.fixture(scope="session")
def pon_fastq_path():
    return "tests/test_data/fastq/"


@pytest.fixture(scope="session")
def panel_bed_file():
    return "tests/test_data/references/panel/panel.bed"


@pytest.fixture(scope="session")
def background_variant_file():
    return "tests/test_data/references/panel/background_variants.txt"


@pytest.fixture(scope="session")
def pon_cnn():
    return "tests/test_data/references/panel/test_panel_ponn.cnn"


@pytest.fixture(scope="session")
def sentieon_license(tmp_path_factory):
    """
    Sentieon's license path fixture
    """
    sentieon_license_dir = tmp_path_factory.mktemp("sentieon_licence")
    sentieon_license_path = sentieon_license_dir / "license_file.lic"
    sentieon_license_path.touch()

    return sentieon_license_path.as_posix()


@pytest.fixture(scope="session")
def sentieon_install_dir(tmp_path_factory):
    """
    Sentieon's license path fixture
    """
    sentieon_install_dir = tmp_path_factory.mktemp("sentieon_install_dir")
    Path(sentieon_install_dir / "bin").mkdir(exist_ok=True)
    sentieon_executable = sentieon_install_dir / "bin" / "sentieon"
    sentieon_executable.touch()

    return sentieon_install_dir.as_posix()


@pytest.fixture()
def no_write_perm_path(tmp_path_factory) -> str:
    """A path with no write permissions."""
    bad_perm_path: Path = tmp_path_factory.mktemp("bad_perm_path")
    bad_perm_path.chmod(0o444)
    return bad_perm_path.as_posix()


@pytest.fixture(scope="session")
def balsamic_cache(tmp_path_factory, reference):
    """
    Create singularity container
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
    """Creates and returns the directory where the case analysis will be saved."""
    analysis_dir = tmp_path_factory.mktemp("analysis", numbered=False)
    return analysis_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only(analysis_dir: str, case_id_tumor_only: str) -> str:
    """Creates and returns the directory containing the FASTQs."""
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
    """Creates and returns the directory containing the FASTQs."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_pon, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_umi(analysis_dir: str, case_id_tumor_only_umi: str) -> str:
    """Creates and returns the directory containing the FASTQs."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_umi, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_normal(analysis_dir: str, case_id_tumor_normal: str) -> str:
    """Creates and returns the directory containing the FASTQs."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_normal, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_only_wgs(analysis_dir: str, case_id_tumor_only_wgs: str) -> str:
    """Creates and returns the directory containing the FASTQs."""
    fastq_dir: Path = Path(analysis_dir, case_id_tumor_only_wgs, "fastq")
    fastq_dir.mkdir(parents=True, exist_ok=True)
    return fastq_dir.as_posix()


@pytest.fixture(scope="session")
def fastq_dir_tumor_normal_wgs(analysis_dir: str, case_id_tumor_normal_wgs: str) -> str:
    """Creates and returns the directory containing the FASTQs."""

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
            ],
        )

    qc_dir = Path(analysis_dir, case_id_tumor_normal, "analysis", "qc")
    qc_dir.mkdir(parents=True, exist_ok=False)
    copy_tree("tests/test_data/qc_files/analysis/qc/", qc_dir.as_posix())

    return Path(
        analysis_dir, case_id_tumor_normal, case_id_tumor_normal + ".json"
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
            ],
        )

    return Path(
        analysis_dir, case_id_tumor_normal_wgs, case_id_tumor_normal_wgs + ".json"
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
            ],
        )

    qc_dir = Path(analysis_dir, case_id_tumor_only, "analysis", "qc")
    qc_dir.mkdir(parents=True, exist_ok=False)
    copy_tree("tests/test_data/qc_files/analysis/qc/", qc_dir.as_posix())

    return Path(
        analysis_dir, case_id_tumor_only, case_id_tumor_only + ".json"
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
            ],
        )

    return Path(
        analysis_dir, case_id_tumor_only_wgs, case_id_tumor_only_wgs + ".json"
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
        analysis_dir, case_id_tumor_only_pon, case_id_tumor_only_pon + ".json"
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
    """Raw requested metrics"""
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
def qc_extracted_metrics(metrics_yaml_path):
    """Extracted and formatted QC metrics"""
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
    """Mocked COSMIC key."""
    return "ZW1haWxAZXhhbXBsZS5jb206bXljb3NtaWNwYXNzd29yZAo="


@pytest.fixture(scope="session", name="cluster_account")
def fixture_cluster_account() -> str:
    """Mocked cluster account for job submission."""
    return "development"


@pytest.fixture(scope="session", name="develop_containers")
def fixture_develop_containers() -> Dict[str, str]:
    """Develop containers fixture."""
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


@pytest.fixture(scope="function", name="cache_config_model_data")
def fixture_cache_config_model_data(
    cache_analysis_model: CacheAnalysisModel,
    develop_containers: Dict[str, str],
    cosmic_key: str,
    timestamp_now: datetime,
    tmp_path: Path,
) -> Dict[str, Any]:
    """Mocked cache config data."""

    return {
        "analysis": cache_analysis_model,
        "references_dir": tmp_path,
        "genome_dir": tmp_path,
        "variants_dir": tmp_path,
        "vep_dir": tmp_path,
        "containers_dir": tmp_path,
        "genome_version": GenomeVersion.HG19,
        "cosmic_key": cosmic_key,
        "bioinfo_tools": BIOINFO_TOOL_ENV,
        "containers": develop_containers,
        "references": REFERENCE_FILES[GenomeVersion.HG19],
        "references_date": timestamp_now.strftime("%Y-%m-%d %H:%M"),
    }


@pytest.fixture(scope="function", name="cache_config_model")
def fixture_cache_config_model(
    cache_config_model_data: Dict[str, dict]
) -> CacheConfigModel:
    """Mocked cache config model."""
    cache_config_model: CacheConfigModel = CacheConfigModel(**cache_config_model_data)
    for reference in cache_config_model.references:
        reference_file: Path = Path(reference[1].file_path)
        reference_file.parent.mkdir(parents=True, exist_ok=True)
        reference_file.touch()
    return cache_config_model


@pytest.fixture(scope="function", name="cache_analysis_model_data")
def fixture_cache_analysis_model_data(case_id_tumor_only: str) -> Dict[str, str]:
    """Mocked cache analysis data."""
    return {"case_id": case_id_tumor_only}


@pytest.fixture(scope="function", name="cache_analysis_model")
def fixture_cache_analysis_model(
    cache_analysis_model_data: Dict[str, str]
) -> CacheAnalysisModel:
    """Mocked cache analysis model."""
    return CacheAnalysisModel(**cache_analysis_model_data)


[
    "/private/var/folders/j0/swnl57794_n2pv9yll88hmnm437kpc/T/pytest-of-vadym.ivanchuk/pytest-320/test_get_refgene_files0/genome/refGene.txt",
    "/private/var/folders/j0/swnl57794_n2pv9yll88hmnm437kpc/T/pytest-of-vadym.ivanchuk/pytest-320/test_get_refgene_files0/genome/refGene.flat",
    "/private/var/folders/j0/swnl57794_n2pv9yll88hmnm437kpc/T/pytest-of-vadym.ivanchuk/pytest-320/test_get_refgene_files0/genome/refGene.flat.bed",
]


@pytest.fixture(scope="function", name="refgene_bed_file")
def fixture_refgene_bed_file(tmp_path: Path) -> Path:
    """Dummy RefSeq's gene BED file."""
    refgene_bed_file: Path = Path(tmp_path, "genome", "refGene.flat.bed")
    refgene_bed_file.touch()
    return refgene_bed_file


@pytest.fixture(scope="function", name="refgene_flat_file")
def fixture_refgene_flat_file(tmp_path: Path) -> Path:
    """Dummy RefSeq's gene flat file."""
    refgene_flat_file: Path = Path(tmp_path, "genome", "refGene.flat")
    refgene_flat_file.touch()
    return refgene_flat_file


@pytest.fixture(scope="function", name="analysis_references_model_data")
def fixture_analysis_references_model_data(
    cache_config_model: CacheConfigModel,
    refgene_bed_file: Path,
    refgene_flat_file: Path,
) -> Dict[str, Path]:
    """Analysis references model data."""
    return {
        "genome_chrom_size": Path(
            cache_config_model.references.genome_chrom_size.file_path
        ),
        "reference_genome": Path(
            cache_config_model.references.reference_genome.file_path
        ),
        "refgene_bed": refgene_bed_file,
        "refgene_flat": refgene_flat_file,
        "refgene_txt": Path(cache_config_model.references.refgene_txt.file_path),
    }


@pytest.fixture(scope="function", name="delly_exclusion_converted_file")
def fixture_delly_exclusion_converted_file(tmp_path: Path) -> Path:
    """Dummy Delly exclusion converted file."""
    reference_file: Path = Path(tmp_path, "genome", "delly_exclusion_converted.tsv")
    reference_file.touch()
    return reference_file


@pytest.fixture(scope="function", name="hg_analysis_references_model_data")
def fixture_hg_analysis_references_model_data(
    cache_config_model: CacheConfigModel,
    analysis_references_model_data: Dict[str, Path],
    delly_exclusion_converted_file: Path,
) -> Dict[str, Path]:
    """Human genome analysis references model data."""
    hg_analysis_references_model_data: Dict[str, Path] = {
        "access_regions": Path(cache_config_model.references.access_regions.file_path),
        "ascat_chr_y_loci": Path(
            cache_config_model.references.ascat_chr_y_loci.file_path
        ),
        "ascat_gc_correction": Path(
            cache_config_model.references.ascat_gc_correction.file_path
        ),
        "clinvar": Path(cache_config_model.references.clinvar.file_path),
        "cosmic": Path(cache_config_model.references.cosmic.file_path),
        "dbsnp": Path(cache_config_model.references.dbsnp.file_path),
        "delly_exclusion": Path(
            cache_config_model.references.delly_exclusion.file_path
        ),
        "delly_exclusion_converted": delly_exclusion_converted_file,
        "delly_mappability": Path(
            cache_config_model.references.delly_mappability.file_path
        ),
        "gnomad_variant": Path(cache_config_model.references.gnomad_variant.file_path),
        "hc_vcf_1kg": Path(cache_config_model.references.hc_vcf_1kg.file_path),
        "known_indel_1kg": Path(
            cache_config_model.references.known_indel_1kg.file_path
        ),
        "mills_1kg": Path(cache_config_model.references.mills_1kg.file_path),
        "rank_score": Path(cache_config_model.references.rank_score.file_path),
        "somalier_sites": Path(cache_config_model.references.somalier_sites.file_path),
        "vcf_1kg": Path(cache_config_model.references.vcf_1kg.file_path),
        "vep_dir": cache_config_model.references_dir,
        "wgs_calling_regions": Path(
            cache_config_model.references.wgs_calling_regions.file_path
        ),
    }
    hg_analysis_references_model_data.update(analysis_references_model_data)
    return hg_analysis_references_model_data


@pytest.fixture(scope="function", name="reference_url")
def fixture_reference_url() -> str:
    """Dummy reference URL."""
    return "gs://gatk-legacy-bundles/b37/reference.vcf.gz"


@pytest.fixture(scope="function", name="reference_file")
def fixture_reference_file(tmp_path: Path) -> Path:
    """Dummy reference file."""
    reference_file: Path = Path(tmp_path, "reference.vcf")
    reference_file.touch()
    return reference_file


@pytest.fixture(scope="function", name="reference_url_model_data")
def fixture_reference_url_model_data(
    reference_url: str, reference_file: Path, cosmic_key: str
) -> Dict[str, Any]:
    """Reference URL model data."""
    return {
        "url": reference_url,
        "file_type": FileType.VCF,
        "gzip": False,
        "file_name": "reference.vcf",
        "dir_name": "variants",
        "file_path": reference_file.as_posix(),
        "secret": cosmic_key,
    }


@pytest.fixture(scope="function", name="references_model_data")
def fixture_references_model_data(
    cache_config_model: CacheConfigModel,
) -> Dict[str, dict]:
    """References model data."""
    return {
        "genome_chrom_size": cache_config_model.references.genome_chrom_size,
        "reference_genome": cache_config_model.references.reference_genome,
        "refgene_sql": cache_config_model.references.refgene_sql,
        "refgene_txt": cache_config_model.references.refgene_txt,
    }


@pytest.fixture(scope="function", name="references_model")
def fixture_references_model(references_model_data: Dict[str, dict]) -> ReferencesModel:
    """Mocked references model."""
    return ReferencesModel(**references_model_data)


@pytest.fixture(scope="function", name="hg_references_model_data")
def fixture_hg_references_model_data(
    cache_config_model: CacheConfigModel,
) -> Dict[str, dict]:
    """Human genome references model data."""
    return dict(cache_config_model.references)


@pytest.fixture(scope="function", name="hg_references_model")
def fixture_hg_references_model(
    hg_references_model_data: Dict[str, dict]
) -> HgReferencesModel:
    """Mocked human genome references model."""
    return HgReferencesModel(**hg_references_model_data)


# @pytest.fixture(scope="function", name="reference_genome_file")
# def fixture_reference_genome_file(tmp_path: Path) -> Path:
#     """Dummy reference file."""
#     reference_genome_file: Path = Path(tmp_path, "reference.fasta")
#     reference_genome_file.touch()
#     return reference_genome_file
#
#
# @pytest.fixture(scope="function", name="refgene_txt_file")
# def fixture_refgene_txt_file(tmp_path: Path) -> Path:
#     """Dummy RefSeq's gene file."""
#     refgene_txt_file: Path = Path(tmp_path, "refGene.txt")
#     refgene_txt_file.touch()
#     return refgene_txt_file
#
#
# @pytest.fixture(scope="function", name="delly_exclusion_file")
# def fixture_delly_exclusion_file(tmp_path: Path) -> Path:
#     """Dummy Delly exclusion file."""
#     delly_exclusion_file: Path = Path(tmp_path, "human.hg19.excl.tsv")
#     delly_exclusion_file.touch()
#     return delly_exclusion_file
#
#
# @pytest.fixture(scope="function", name="delly_exclusion_converted_file")
# def fixture_delly_exclusion_converted_file(tmp_path: Path) -> Path:
#     """Dummy Delly exclusion converted file."""
#     delly_exclusion_converted: Path = Path(tmp_path, "human.hg19.excl_converted.tsv")
#     delly_exclusion_converted.touch()
#     return delly_exclusion_converted
#
#
# @pytest.fixture(scope="function", name="delly_mappability_file")
# def fixture_delly_mappability_file(tmp_path: Path) -> Path:
#     """Dummy Delly mappability file."""
#     delly_mappability_file: Path = Path(tmp_path, "GRCh37.delly.blacklist.gz")
#     delly_mappability_file.touch()
#     return delly_mappability_file
#
#
# @pytest.fixture(scope="function", name="delly_mappability_findex_file")
# def fixture_delly_mappability_findex_file(tmp_path: Path) -> Path:
#     """Dummy Delly mappability findex file."""
#     delly_mappability_findex_file: Path = Path(
#         tmp_path, "GRCh37.delly.blacklist.gz.fai"
#     )
#     delly_mappability_findex_file.touch()
#     return delly_mappability_findex_file
#
#
# @pytest.fixture(scope="function", name="delly_mappability_gindex_file")
# def fixture_delly_mappability_gindex_file(tmp_path: Path) -> Path:
#     """Dummy Delly exclusion file."""
#     delly_mappability_gindex_file: Path = Path(
#         tmp_path, "GRCh37.delly.blacklist.gz.gzi"
#     )
#     delly_mappability_gindex_file.touch()
#     return delly_mappability_gindex_file
#
#
# @pytest.fixture(scope="function", name="gnomad_variant_file")
# def fixture_gnomad_variant_file(tmp_path: Path) -> Path:
#     """Dummy gnomad_variant file."""
#     gnomad_variant_file: Path = Path(tmp_path, "gnomad.genomes.r2.1.1.sites.vcf.bgz")
#     gnomad_variant_file.touch()
#     return gnomad_variant_file
#
#
# @pytest.fixture(scope="function", name="gnomad_variant_index_file")
# def fixture_gnomad_variant_index_file(tmp_path: Path) -> Path:
#     """Dummy gnomad_variant index file."""
#     gnomad_variant_index_file: Path = Path(
#         tmp_path, "gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi"
#     )
#     gnomad_variant_index_file.touch()
#     return gnomad_variant_index_file
#
#
# @pytest.fixture(scope="function", name="known_indel_1kg_file")
# def fixture_known_indel_1kg_file(tmp_path: Path) -> Path:
#     """Dummy 1000 Genome known InDels VCF file."""
#     known_indel_1kg_file: Path = Path(tmp_path, "1kg_known_indels_b37.vcf")
#     known_indel_1kg_file.touch()
#     return known_indel_1kg_file
#
#
# @pytest.fixture(scope="function", name="mills_1kg_file")
# def fixture_mills_1kg_file(tmp_path: Path) -> Path:
#     """Dummy Mills' high confidence InDels VCF file."""
#     mills_1kg_file: Path = Path(tmp_path, "mills_1kg_index.vcf")
#     mills_1kg_file.touch()
#     return mills_1kg_file
#
#
# @pytest.fixture(scope="function", name="hc_vcf_1kg_file")
# def fixture_hc_vcf_1kg_file(tmp_path: Path) -> Path:
#     """Dummy high confidence 1000 Genome VCF file."""
#     hc_vcf_1kg_file: Path = Path(tmp_path, "1kg_phase1_snps_high_confidence_b37.vcf")
#     hc_vcf_1kg_file.touch()
#     return hc_vcf_1kg_file
#
#
# @pytest.fixture(scope="function", name="vcf_1kg_file")
# def fixture_vcf_1kg_file(tmp_path: Path) -> Path:
#     """Dummy 1000 Genome all SNPs file."""
#     vcf_1kg_file: Path = Path(tmp_path, "1k_genome_wgs_p1_v3_all_sites.vcf")
#     vcf_1kg_file.touch()
#     return vcf_1kg_file
