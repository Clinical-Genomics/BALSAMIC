import pytest
import json
import os

from unittest import mock
from distutils.dir_util import copy_tree
from pathlib import Path
from functools import partial

from _pytest.tmpdir import TempPathFactory

from BALSAMIC.constants.workflow_params import VCF_DICT
from click.testing import CliRunner

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
    Path(fastq_dir, f"ACC1_XXXXX_R_1.fastq.gz").touch()
    Path(fastq_dir, f"ACC1_XXXXX_R_2.fastq.gz").touch()
    Path(fastq_dir, f"ACC2_XXXXX_R_1.fastq.gz").touch()
    Path(fastq_dir, f"ACC2_XXXXX_R_2.fastq.gz").touch()
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
def config_files():
    """dict: path of the config files"""
    return {
        "sample": "BALSAMIC/config/sample.json",
        "analysis_paired": "BALSAMIC/config/analysis_paired.json",
        "cluster_json": "BALSAMIC/config/cluster.json",
        "analysis_paired_umi": "BALSAMIC/config/analysis_paired_umi.json",
        "analysis_single": "BALSAMIC/config/analysis_single.json",
        "analysis_single_umi": "BALSAMIC/config/analysis_single_umi.json",
        "panel_bed_file": "tests/test_data/references/panel/panel.bed",
        "background_variant_file": "tests/test_data/references/panel/background_variants.txt",
        "pon_cnn": "tests/test_data/references/panel/test_panel_ponn.cnn",
        "pon_fastq_path": "tests/test_data/fastq/",
    }


@pytest.fixture(scope="session")
def reference():
    """reference json model"""
    return {
        "reference": {
            "reference_genome": "tests/test_data/references/genome/human_g1k_v37_decoy.fasta",
            "dbsnp": "tests/test_data/references/variants/dbsnp_grch37_b138.vcf.gz",
            "1kg_snps_all": "tests/test_data/references/variants/1k_genome_wgs_p1_v3_all_sites.vcf.gz",
            "1kg_snps_high": "tests/test_data/references/variants/1kg_phase1_snps_high_confidence_b37.vcf.gz",
            "1kg_known_indel": "tests/test_data/references/variants/1kg_known_indels_b37.vcf.gz",
            "mills_1kg": "tests/test_data/references/variants/mills_1kg_index.vcf.gz",
            "gnomad_variant": "tests/test_data/reference/variants/gnomad.genomes.r2.1.1.sites.vcf.bgz",
            "cosmic": "tests/test_data/references/variants/cosmic_coding_muts_v89.vcf.gz",
            "vep": "tests/test_data/references/vep/",
            "refflat": "tests/test_data/references/genome/refseq.flat",
            "refGene": "tests/test_data/references/genome/refGene.txt",
            "wgs_calling_interval": "tests/test_data/references/genome/wgs_calling_regions.v1",
            "genome_chrom_size": "tests/test_data/references/genome/hg19.chrom.sizes",
            "exon_bed": "tests/test_data/references/genome/refseq.flat.bed",
            "rankscore": "tests/test_data/references/genome/cancer_rank_model_-v0.1-.ini",
            "access_regions": "tests/test_data/references/genome/access-5k-mappable.hg19.bed",
            "delly_exclusion": "tests/test_data/references/genome/delly_exclusion.tsv",
            "delly_exclusion_converted": "tests/test_data/references/genome/delly_exclusion_converted.tsv",
            "delly_mappability": "tests/test_data/references/genome/delly_mappability.gz",
            "delly_mappability_gindex": "tests/test_data/references/genome/delly_mappability.gz.gzi",
            "delly_mappability_findex": "tests/test_data/references/genome/delly_mappability.fai",
            "ascat_gccorrection": "tests/test_data/references/genome/GRCh37_SnpGcCorrections.tsv",
            "ascat_chryloci": "tests/test_data/references/genome/GRCh37_Y.loci",
            "clinvar": "tests/test_data/references/genome/clinvar.vcf.gz",
            "clinical_snv_observations": "tests/test_data/references/variants/clinical_snv_variants.vcf.gz",
            "clinical_sv_observations": "tests/test_data/references/variants/clinical_sv_variants.vcf.gz",
            "swegen_snv_frequency": "tests/test_data/references/variants/swegen_snv.vcf.gz",
            "swegen_sv_frequency": "tests/test_data/references/variants/swegen_sv.vcf.gz",
            "somalier_sites": "tests/test_data/references/variants/GRCh37.somalier.sites.vcf.gz",
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

    cache_dir = tmp_path_factory.mktemp("balsmic_coche")

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

    # Fill the analysis fastq path folder with the concatenated fastq files
    analysis_fastq_dir = Path(analysis_dir, case_id_tumor_only, "analysis", "fastq")
    analysis_fastq_dir.mkdir(parents=True, exist_ok=True)
    for read in [1, 2]:
        Path(analysis_fastq_dir, f"concatenated_ACC1_R_{read}.fastq.gz").touch()

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
def sample_config():
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
            "S1_R": {
                "file_prefix": "S1_R",
                "type": "tumor",
                "readpair_suffix": ["1", "2"],
            },
            "S2_R": {
                "file_prefix": "S2_R",
                "type": "normal",
                "readpair_suffix": ["1", "2"],
            },
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
        "fastq",
        "concatenated_ACC1_R_{read}.fastq.gz",
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
