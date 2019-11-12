#! /usr/bin/python

import pytest
import yaml
import json
import os

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
        "sample": "BALSAMIC/config/sample.json",
        "reference": "tests/test_data/references/reference.json",
        "analysis_paired": "BALSAMIC/config/analysis_paired.json",
        "cluster_json": "BALSAMIC/config/cluster.json",
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
def singularity_container(tmp_path_factory):
    """
    Create singularity container
    """
    
    container_dir = tmp_path_factory.mktemp("test_container")
    container_file = container_dir / "singularity_container.simg"

    return str(container_file)


@pytest.fixture(scope='session')
def analysis_dir(tmp_path_factory):
    """
    Creates and returns analysis directory
    """
    analysis_dir = tmp_path_factory.mktemp('analysis')
    return analysis_dir

@pytest.fixture(scope='session')
def snakemake_job_script(tmp_path_factory, tumor_normal_config):
    """
    Creates a dummy snakemake jobscript
    """
    case_name = 'job_submit_test_case'
    with open(tumor_normal_config, 'r') as input_config:
        sample_config = json.load(input_config)
    
    script_dir = tmp_path_factory.mktemp('snakemake_script')
    snakemake_script_file = script_dir / 'example_script.sh' 
    snakemake_script = '''#!/bin/sh
# properties = {"type": "single", "rule": "all", "local": false, "input": ["dummy_path"], "output": ["dummy_path"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 0, "cluster": {"name": "BALSAMIC.all.", "time": "00:15:00", "n": 1, "mail_type": "END", "partition": "core"}}
ls -l # dummy command
'''
    snakemake_script_file.touch()
    with open(snakemake_script_file, 'w') as fn:
        fn.write(snakemake_script)

    return {
      "snakescript": str(snakemake_script_file)
      }

@pytest.fixture(scope='session')
def tumor_normal_config(tmp_path_factory, sample_fastq, analysis_dir, singularity_container):
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
        str(normal), '--case-id', case_name,
        '--singularity', singularity_container,  
        '--analysis-dir',
        str(analysis_dir), '--output-config', sample_config_file_name,
        '--reference-config', reference_json
    ])

    return str(analysis_dir / case_name / sample_config_file_name)


@pytest.fixture(scope='session')
def tumor_normal_wgs_config(tmp_path_factory, sample_fastq, analysis_dir, singularity_container):
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
        str(tumor), '-n',
        str(normal), '--case-id', case_name,
        '--singularity', singularity_container,  
        '--analysis-dir',
        str(analysis_dir), '--output-config', sample_config_file_name,
        '--reference-config', reference_json
    ])

    return str(analysis_dir / case_name / sample_config_file_name)


@pytest.fixture(scope='session')
def tumor_only_config(tmpdir_factory, sample_fastq, singularity_container):
    """
    invokes balsamic config sample -t xxx to create sample config
    for tumor only
    """
    case_name = 'sample_tumor_only'
    analysis_dir = fn = tmpdir_factory.mktemp('analysis', numbered=False)
    tumor = sample_fastq['tumor']
    panel_bed_file = 'tests/test_data/references/panel/panel.bed'
    reference_json = 'tests/test_data/references/reference.json'
    sample_config_file_name = 'sample.json'

    runner = CliRunner()
    result = runner.invoke(cli, [
        'config', 'case', '-p', panel_bed_file, '-t',
        str(tumor), '--case-id', case_name, '--analysis-dir',
        str(analysis_dir), '--output-config', sample_config_file_name,
        '--singularity', singularity_container,  
        '--reference-config', reference_json
    ])

    dummy_log = tmpdir_factory.mktemp('analysis/' + case_name + '/logs', numbered=False).join('example_file').write('not_empty')

    return os.path.join(analysis_dir, case_name, sample_config_file_name)


@pytest.fixture(scope='session')
def tumor_only_wgs_config(tmp_path_factory, sample_fastq, analysis_dir, singularity_container):
    """
    invokes balsamic config sample -t xxx to create sample config
    for tumor only
    """
    case_name = 'sample_tumor_only'
    tumor = sample_fastq['tumor']
    reference_json = 'tests/test_data/references/reference.json'
    sample_config_file_name = 'sample.json'

    runner = CliRunner()
    result = runner.invoke(cli, [
        'config', 'case', '-t',
        str(tumor), '--case-id', case_name, '--analysis-dir',
        str(analysis_dir), '--output-config', sample_config_file_name,
        '--singularity', singularity_container,  
        '--reference-config', reference_json
    ])

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
