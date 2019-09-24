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
        "env_py27": ["python", "strelka", "manta", "tabix"],
        "env_py36": [
            "python", "pip", "bcftools", "bwa", "fastqc", "sambamba",
            "samtools", "tabix", "gatk", "picard", "fgbio", "freebayes",
            "vardict", "vardict-java", "ensembl-vep", "cnvkit", "pindel",
            "multiqc", "bedtools", "fastp"
        ]
    }

    conda_env_file = conda_env_path / "test_BALSAMIC_env.yaml"

    yaml.dump(conda_packages, open(conda_env_file, 'w'))

    return str(conda_env_file)


@pytest.fixture(scope='session')
def no_write_perm_path(tmp_path_factory):
    """
    A path with no write permission
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
    fastq_dir = tmp_path_factory.mktemp("fastq")
    fastq_valid = fastq_dir / "S1_R_1.fastq.gz"
    fastq_invalid = fastq_dir / "sample.fastq.gz"

    # dummy tumor fastq file
    tumor_fastq_R_1 = fastq_dir / "tumor_R_1.fastq.gz"
    tumor_fastq_R_2 = fastq_dir / "tumor_R_2.fastq.gz"

    # dummy normal fastq file
    normal_fastq_R_1 = fastq_dir / "normal_R_1.fastq.gz"
    normal_fastq_R_2 = fastq_dir / "normal_R_2.fastq.gz"

    for fastq_file in (fastq_valid, fastq_invalid, tumor_fastq_R_1,
                       tumor_fastq_R_2, normal_fastq_R_1, normal_fastq_R_2):
        fastq_file.touch()

    return {
        "fastq_valid": str(fastq_valid.absolute()),
        "fastq_invalid": str(fastq_invalid.absolute()),
        "tumor": str(tumor_fastq_R_1.absolute()),
        "normal": str(normal_fastq_R_2.absolute())
    }


@pytest.fixture(scope='session')
def analysis_dir(tmp_path_factory):
    """
    Creates and returns analysis directory
    """
    analysis_dir = tmp_path_factory.mktemp('analysis')
    return analysis_dir


@pytest.fixture(scope='session')
def tumor_normal_config(tmp_path_factory, sample_fastq, analysis_dir):
    """
    invokes balsamic config sample -t xxx -n xxx to create sample config
    for tumor-normal
    """
    case_name = 'sample_tumor_normal'
    tumor = sample_fastq['tumor']
    normal = sample_fastq['normal']
    panel_bed_file = 'tests/test_data/references/panel/panel.bed'
    reference_json = 'tests/test_data/references/reference.json'
    sample_config_file_name = 'sample.json'

    runner = CliRunner()
    result = runner.invoke(cli, [
        'config', 'case', '-p', panel_bed_file, '-t',
        str(tumor), '-n',
        str(normal), '--case-id', case_name, '--analysis-dir',
        str(analysis_dir), '--output-config', sample_config_file_name,
        '--reference-config', reference_json
    ])

    return str(analysis_dir / case_name / sample_config_file_name)


@pytest.fixture(scope='session')
def tumor_normal_wgs_config(tmp_path_factory, sample_fastq, analysis_dir):
    """
    invokes balsamic config sample -t xxx -n xxx to create sample config
    for tumor-normal
    """
    case_name = 'sample_tumor_normal'
    tumor = sample_fastq['tumor']
    normal = sample_fastq['normal']
    reference_json = 'tests/test_data/references/reference.json'
    sample_config_file_name = 'sample.json'

    runner = CliRunner()
    result = runner.invoke(cli, [
        'config', 'case', '-t',
        str(tumor), '-n', str(normal), '--case-id', case_name, '--analysis-dir',
        str(analysis_dir), '--output-config', sample_config_file_name,
        '--reference-config', reference_json
    ])

    return str(analysis_dir / case_name / sample_config_file_name)


@pytest.fixture(scope='session')
def tumor_only_config(tmp_path_factory, sample_fastq, analysis_dir):
    """
    invokes balsamic config sample -t xxx to create sample config
    for tumor only
    """
    case_name = 'sample_tumor_only'
    tumor = sample_fastq['tumor']
    panel_bed_file = 'tests/test_data/references/panel/panel.bed'
    reference_json = 'tests/test_data/references/reference.json'
    sample_config_file_name = 'sample.json'

    runner = CliRunner()
    result = runner.invoke(cli, [
        'config', 'case', '-p', panel_bed_file, '-t',
        str(tumor), '--case-id', case_name, '--analysis-dir',
        str(analysis_dir), '--output-config', sample_config_file_name,
        '--reference-config', reference_json
    ])

    return str(analysis_dir / case_name / sample_config_file_name)


@pytest.fixture(scope='session')
def tumor_only_wgs_config(tmp_path_factory, sample_fastq, analysis_dir):
    """
    invokes balsamic config sample -t xxx to create sample config
    for tumor only
    """
    case_name = 'sample_tumor_only'
    tumor = sample_fastq['tumor']
    reference_json = 'tests/test_data/references/reference.json'
    sample_config_file_name = 'sample.json'

    runner = CliRunner()
    result = runner.invoke(cli, ['config', 'case', '-t', str(tumor),
                                 '--case-id', case_name, '--analysis-dir', str(analysis_dir),
                                 '--output-config', sample_config_file_name, '--reference-config',
                                 reference_json])

    return str(analysis_dir / case_name / sample_config_file_name)


@pytest.fixture(scope='session')
def sample_config():
    """
    sample config dict to test workflow utils
    """
    sample_config = {
        "QC": {
            "picard_rmdup": "False",
            "adapter":
            "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
            "min_seq_length": "25",
            "quality_trim": "True",
            "adapter_trim": "False",
            "umi_trim": "True",
            "umi_trim_length": "5"
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
            "dag":
            "tests/test_data/id1/id1_analysis.json_BALSAMIC_2.9.8_graph.pdf"
        },
        "vcf": {
            "manta": {
                "mutation": "somatic",
                "type": "SV"
            },
            "vardict": {
                "mutation": "somatic",
                "type": "SNV"
            },
            "pindel": {
                "mutation": "somatic",
                "type": "SV"
            },
            "strelka": {
                "mutation": "somatic",
                "type": "SNV"
            },
            "mutect": {
                "mutation": "somatic",
                "type": "SNV"
            },
            "tnscope": {
                "mutation": "somatic",
                "type": "SNV"
            },
            "tnsnv": {
                "mutation": "somatic",
                "type": "SNV"
            },
            "tnhaplotyper": {
                "mutation": "somatic",
                "type": "SNV"
            },
            "dnascope": {
                "mutation": "germline",
                "type": "SNV"
            },
            "manta_germline": {
                "mutation": "germline",
                "type": "SV"
            },
            "haplotypecaller": {
                "mutation": "germline",
                "type": "SNV"
            },
            "strelka_germline": {
                "mutation": "germline",
                "type": "SNV"
            }
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


#
#{
#    "analysis": {
#        "case_id":
#        "test_sample",
#        "analysis_type":
#        "paired",
#        "analysis_dir":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/",
#        "fastq_path":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_sample/analysis/fastq/",
#        "script":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_sample/scripts/",
#        "log":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_sample/logs/",
#        "result":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_sample/analysis",
#        "config_creation_date":
#        "2019-09-13 11:15",
#        "BALSAMIC_version":
#        "3.0.0",
#        "dag":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_sample/test_sample.json_BALSAMIC_3.0.0_graph.pdf"
#    },
#    "reference": {
#        "reference_genome":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/genome/human_g1k_v37_decoy.fasta",
#        "dbsnp":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/variants/dbsnp_grch37_b138.vcf.gz",
#        "1kg_snps_all":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/variants/1k_genome_wgs_p1_v3_all_sites.vcf.gz",
#        "1kg_snps_high":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/variants/1kg_phase1_snps_high_confidence_b37.vcf.gz",
#        "mills_1kg":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/variants/mills_1kg_index.vcf.gz",
#        "cosmic":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/variants/cosmic_coding_muts_v89.vcf.gz",
#        "vep":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/vep",
#        "refflat":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/genome/refseq.flat",
#        "exon_bed":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/genome/refseq.flat.bed"
#    },
#    "conda_env_yaml":
#    "/home/proj/long-term-stage/cancer/BALSAMIC/BALSAMIC_env.yaml",
#    "rule_directory": "/home/proj/long-term-stage/cancer/BALSAMIC/BALSAMIC/",
#    "bioinfo_tools": {
#        "bcftools": "1.9.0",
#        "fastqc": "0.11.5",
#        "gatk": "3.8",
#        "sambamba": "0.6.6",
#        "strelka": "2.8.4",
#        "bwa": "0.7.15",
#        "cutadapt": "1.15",
#        "samtools": "1.6",
#        "picard": "2.17.0",
#        "tabix": "0.2.5",
#        "manta": "1.3.0"
#    },
#    "panel": {
#        "capture_kit":
#        "/home/proj/long-term-stage/cancer/BALSAMIC/tests/test_data/references/panel/panel.bed",
#        "chrom": [
#            "11", "5", "14", "12", "17", "9", "19", "2", "7", "20", "13", "1",
#            "22", "6", "10", "4", "15", "18", "8", "16", "3", "21"
#        ]
#    }
#}
