import pytest

from functools import partial
from click.testing import CliRunner

from BALSAMIC.commands import cli


@pytest.fixture
def cli_runner():
    """ click - cli testing """
    runner = CliRunner()
    return runner


@pytest.fixture
def invoke_cli(cli_runner):
    """ invoking cli commands with options"""
    return partial(cli_runner.invoke, cli)


@pytest.fixture(scope='session')
def config_files():
    """ dict: path of the config files """
    return {
        "install": "BALSAMIC/config/install.json",
        "reference": "BALSAMIC/config/reference.json",
        "sample": "BALSAMIC/config/sample.json",
        "analysis_paired": "BALSAMIC/config/analysis_paired.json",
        "analysis_paired_umi": "BALSAMIC/config/analysis_paired_umi.json",
        "analysis_single": "BALSAMIC/config/analysis_single.json",
        "analysis_single_umi": "BALSAMIC/config/analysis_single_umi.json",
        "panel_bed_file": "tests/test_data/references/GRCh37/panel/panel.bed",
        "test_reference": "tests/test_data/references/reference.json"
    }


@pytest.fixture(scope='session')
def conda_yaml():
    """
    conda env config file paths
    """
    return {
        "balsamic-base": "BALSAMIC/conda_yaml/BALSAMIC-base.yaml",
        "balsamic-p27": "BALSAMIC/conda_yaml/BALSAMIC-py27.yaml",
        "balsamic-py36": "BALSAMIC/conda_yaml/BALSAMIC-py36.yaml",
    }


@pytest.fixture(scope='session')
def sample_config():
    """
    sample config dict to test workflow utils
    """
    sample_config = {
        "QC": {
            "picard_rmdup": "TRUE",
            "adapter":
            "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
            "min_seq_length": "25"
        },
        "vcf": {
            "manta": {
                "default": [
                    "diploidSV.vcf.gz", "somaticSV.vcf.gz",
                    "candidateSV.vcf.gz", "candidateSmallIndels.vcf.gz"
                ],
                "merged":
                "manta.vcf.gz",
                "mutation":
                "somatic",
                "type":
                "SV"
            },
            "manta_germline": {
                "default": [
                    "diploidSV.vcf.gz", "candidateSV.vcf.gz",
                    "candidateSmallIndels.vcf.gz"
                ],
                "mutation":
                "germline",
                "merged":
                "manta_germline.vcf.gz",
                "type":
                "SV"
            },
            "strelka_germline": {
                "default": ["variants.vcf.gz", "germline.S1.vcf.gz"],
                "mutation": "germline",
                "merged": "strelka_germline.vcf.gz",
                "type": "SNV"
            },
            "strelka": {
                "default": ["somatic.snvs.vcf.gz", "somatic.indels.vcf.gz"],
                "mutation": "somatic",
                "merged": "strelka.vcf.gz",
                "type": "SNV"
            },
            "mutect": {
                "default": "mutect.vcf.gz",
                "mutation": "somatic",
                "merged": "mutect.vcf.gz",
                "type": "SNV"
            },
            "freebayes": {
                "default": "freebayes.vcf.gz",
                "mutation": "germline",
                "merged": "freebayes.vcf.gz",
                "type": "SNV"
            },
            "haplotypecaller": {
                "default": "haplotypecaller.vcf.gz",
                "mutation": "germline",
                "merged": "haplotypecaller.vcf.gz",
                "type": "SNV"
            },
            "vardict": {
                "default": "vardict.vcf.gz",
                "mutation": "somatic",
                "merged": "vardict.vcf.gz",
                "type": "SNV"
            }
        },
        "analysis": {
            "sample_id":
            "id1",
            "analysis_type":
            "paired",
            "analysis_dir":
            "/BALSAMIC/tests/test_data/",
            "fastq_path":
            "/BALSAMIC/tests/test_data/id1/fastq/",
            "script":
            "/BALSAMIC/tests/test_data/id1/scripts/",
            "log":
            "/BALSAMIC/tests/test_data/id1/logs/",
            "result":
            "/BALSAMIC/tests/test_data/id1/analysis/",
            "config_creation_date":
            "yyyy-mm-dd xx",
            "BALSAMIC_version":
            "2.9.8",
            "dag":
            "/BALSAMIC/tests/test_data/id1/id1_analysis.json_BALSAMIC_2.9.8_graph.pdf"
        },
        "samples": {
            "S1_R": {
                "file_prefix": "S1_R",
                "type": "tumor",
                "readpair_suffix": ["1", "2"]
            },
            "S2_R": {
                "file_prefix": "S2_R",
                "type": "normal",
                "readpair_suffix": ["1", "2"]
            }
        },
    }

    return sample_config
