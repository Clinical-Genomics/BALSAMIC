#! /usr/bin/env python3
# syntax=python tabstop=4 expandtab

import snakemake

from BALSAMIC.utils.cli import get_snakefile, get_abs_path


def test_variant_calling_paired(sample_config):
    # GIVEN a sample config dict and snakefile
    snakefile = get_snakefile('single')
    config_json = 'tests/test_data/config.json'
    rules = {'rule_directory': get_abs_path(config_json)}

    # WHEN invoking snakemake module
    assert snakemake.snakemake(snakefile,
                               config=rules,
                               configfile=config_json,
                               dryrun=True)

