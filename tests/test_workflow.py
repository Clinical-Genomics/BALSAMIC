#! /usr/bin/env python3
# syntax=python tabstop=4 expandtab

import snakemake

from BALSAMIC.utils.cli import get_snakefile, get_abs_path


def test_workflow(sample_config):
    # GIVEN a sample config dict and snakefile
    workflows = ['paired', 'qc', 'single']
    config_json = 'tests/test_data/config.json'
    rules = {'rule_directory': get_abs_path(config_json)}

    # WHEN invoking snakemake module with dryrun option
    # THEN it should return true
    for workflow in workflows:
        snakefile = get_snakefile(workflow)
        assert snakemake.snakemake(snakefile,
                                   config=rules,
                                   configfile=config_json,
                                   dryrun=True)

