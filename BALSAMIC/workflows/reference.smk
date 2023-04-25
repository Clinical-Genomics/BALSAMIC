# syntax=python tabstop=4 expandtab
# coding: utf-8

import logging
from pathlib import Path

from BALSAMIC.constants.paths import BALSAMIC_DIR
from BALSAMIC.constants.workflow_rules import SNAKEMAKE_RULES
from BALSAMIC.models.cache_models import CacheConfigModel
from BALSAMIC.utils.io import write_finish_file
from BALSAMIC.utils.rule import get_threads

LOG = logging.getLogger(__name__)

# Shell options
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# Balsamic cache configuration model
cache_config: CacheConfigModel = CacheConfigModel.parse_obj(config)

# Rules to include
for rule in SNAKEMAKE_RULES["cache"][cache_config.genome_version]:
    include: Path(BALSAMIC_DIR, rule).as_posix()

LOG.info(f"The rules {SNAKEMAKE_RULES['cache'][cache_config.genome_version]} will be included in the reference workflow")


rule all:
    """Target rule for Balsamic cache generation."""
    input:
        expand(
            f"{cache_config.containers_dir.as_posix()}/{{singularity_image}}.sif",
            singularity_image=cache_config.containers.keys(),
        ),
        expand(
            f"{cache_config.references_dir.as_posix()}/{{reference_file}}",
            reference_file=cache_config.get_reference_paths(),
        ),
    output:
        finish_file=f"{cache_config.references_dir.as_posix()}/reference.finish",
    threads:
        get_threads(cluster_config, "all")
    run:
        write_finish_file(output.finish_file)
