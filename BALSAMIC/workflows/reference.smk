# syntax=python tabstop=4 expandtab
# coding: utf-8
import json
import logging
import os
from pathlib import Path

from BALSAMIC.constants.cache import FileType, Species, VEP_PLUGINS

from BALSAMIC.constants.paths import BALSAMIC_DIR, REFSEQ_SCRIPT_PATH
from BALSAMIC.constants.workflow_rules import SNAKEMAKE_RULES
from BALSAMIC.models.cache_models import CacheConfigModel
from BALSAMIC.utils.io import write_finish_file, write_json
from BALSAMIC.utils.rule import get_threads

LOG = logging.getLogger(__name__)

# Balsamic cache configuration model
cache_config: CacheConfigModel = CacheConfigModel.parse_obj(config)
analysis_references_json: dict = json.loads(
    cache_config.get_analysis_references().json()
)
reference_finish_path: str = Path(
    cache_config.references_dir, "reference.json"
).as_posix()

# Temporary directory and shell options
os.environ["TMPDIR"] = cache_config.references_dir.as_posix()
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

# Rules to include
for rule in SNAKEMAKE_RULES["cache"][cache_config.genome_version]:

    include: Path(BALSAMIC_DIR, rule).as_posix()


LOG.info(
    f"The rules {SNAKEMAKE_RULES['cache'][cache_config.genome_version]} will be included in the reference workflow"
)


rule all:
    """Target rule for Balsamic cache generation."""
    input:
        cache_config.get_container_output_paths(),
        cache_config.get_reference_output_paths(),
    output:
        finish_file=f"{cache_config.references_dir.as_posix()}/reference.finish",
    threads: get_threads(cluster_config, "all")
    run:
        write_json(json_obj=analysis_references_json, path=reference_finish_path)
        write_finish_file(file_path=output.finish_file)
