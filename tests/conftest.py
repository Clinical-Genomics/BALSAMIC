#! /usr/bin/python

import pytest
import yaml
import json

from pathlib import Path
from functools import partial
from click.testing import CliRunner

from BALSAMIC.commands.base import cli


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
        "sample": "BALSAMIC/config/sample.json",
        "reference": "tests/test_data/references/reference.json",
        "analysis_paired": "BALSAMIC/config/analysis_paired.json",
        "analysis_paired_umi": "BALSAMIC/config/analysis_paired_umi.json",
        "analysis_single": "BALSAMIC/config/analysis_single.json",
        "analysis_single_umi": "BALSAMIC/config/analysis_single_umi.json",
        "panel_bed_file": "tests/test_data/references/panel/panel.bed",
        "test_reference": "tests/test_data/references/reference.json"
    }


@pytest.fixture(scope='session')
def conda():
    """
    conda env config file paths
    """
    return {
        "balsamic-base": "BALSAMIC/conda/BALSAMIC-base.yaml",
        "balsamic-p27": "BALSAMIC/conda/BALSAMIC-py27.yaml",
        "balsamic-py36": "BALSAMIC/conda/BALSAMIC-py36.yaml",
    }


@pytest.fixture(scope='session')
def BALSAMIC_env(tmp_path_factory):
    """
    Writes BALSAMIC_env.yaml file.
    """
    # create a conda_env directory
    conda_env_path = tmp_path_factory.mktemp("conda_env")

    # create a yaml file inside conda_env_path
    conda_packages = {
        str("env_1"): ["package_1", "package_2"],
        str("env_2"): ["package_4"]
    }

    conda_env_file = conda_env_path / "test_BALSAMIC_env.yaml"

    yaml.dump(conda_packages, open(conda_env_file, 'w'))

    return str(conda_env_file)


@pytest.fixture(scope='session')
def no_write_perm_path(tmp_path_factory):
    """
    Writes BALSAMIC_env.yaml file.
    """
    # create a conda_env directory
    bad_perm_path = tmp_path_factory.mktemp("bad_perm_path")
    
    Path(bad_perm_path).chmod(0o444)

    return str(bad_perm_path)


@pytest.fixture(scope='session')
def sample_fastq(tmp_path_factory):
    """
    create sample fastq files
    """
    tmp_dir = tmp_path_factory.mktemp("fastq")
    fastq_valid = tmp_dir / "S1_R_1.fastq.gz"
    fastq_invalid = tmp_dir / "sample.fastq.gz"

    for file in (fastq_valid, fastq_invalid):
        with open(file, 'w') as fh:
            fh.write("ATCG")
            fh.close()

    return {"fastq_valid": fastq_valid, "fastq_invalid": fastq_invalid}


@pytest.fixture(scope='session')
def install_config(BALSAMIC_env):
    """
    Fixture for install.json
    """
    install_json = {
      "conda_env_yaml": BALSAMIC_env,
      "rule_directory": str(Path('BALSAMIC/snakemake_rules').absolute())
    } 
    
    with open('BALSAMIC/config/install.json', 'w') as install_json_file:
        json.dump(install_json, install_json_file)


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
                "merged": "manta.vcf.gz",
                "mutation": "somatic",
                "type": "SV"
            },
            "manta_germline": {
                "default": [
                    "diploidSV.vcf.gz", "candidateSV.vcf.gz",
                    "candidateSmallIndels.vcf.gz"
                ],
                "mutation": "germline",
                "merged": "manta_germline.vcf.gz",
                "type": "SV"
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
            "sample_id": "id1",
            "analysis_type": "paired",
            "analysis_dir": "tests/test_data/",
            "fastq_path": "tests/test_data/id1/fastq/",
            "script": "tests/test_data/id1/scripts/",
            "log": "tests/test_data/id1/logs/",
            "result": "tests/test_data/id1/analysis/",
            "config_creation_date": "yyyy-mm-dd xx",
            "BALSAMIC_version": "2.9.8",
            "dag": "tests/test_data/id1/id1_analysis.json_BALSAMIC_2.9.8_graph.pdf"
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
        }
    }

    return sample_config
