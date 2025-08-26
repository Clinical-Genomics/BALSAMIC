"""Report deliver CLI command."""
import logging
from pathlib import Path
from typing import Any, Dict, List

import click
import subprocess
import snakemake
from BALSAMIC.constants.constants import FileType

from BALSAMIC.commands.options import (
    OPTION_RULES_TO_DELIVER,
    OPTION_SAMPLE_CONFIG,
)
from BALSAMIC.models.config import ConfigModel
from BALSAMIC.utils.cli import (
    get_snakefile,
    convert_deliverables_tags,
    get_file_extension,
)
from BALSAMIC.utils.delivery import get_multiqc_deliverables
from BALSAMIC.utils.io import read_json, write_json, write_yaml

LOG = logging.getLogger(__name__)


@click.command(
    "deliver", short_help="Create a <case_id>.hk file with output analysis files"
)
@OPTION_RULES_TO_DELIVER
@OPTION_SAMPLE_CONFIG
@click.pass_context
def deliver(
    context: click.Context,
    rules_to_deliver: List[str],
    sample_config: str,
):
    """Report deliver command to generate output analysis files."""
    LOG.info(f"BALSAMIC started with log level {context.obj['log_level']}.")
    LOG.info("Creating <case_id>.hk deliverables file")
    config: Dict[str, Any] = read_json(sample_config)
    config_model: ConfigModel = ConfigModel(**config)
    output_dir: Path = Path(config_model.analysis.result, "delivery_report")
    output_dir.mkdir(exist_ok=True)

    snakefile: Path = get_snakefile(
        analysis_type=config_model.analysis.analysis_type,
        analysis_workflow=config_model.analysis.analysis_workflow,
    )

    LOG.info(f"Delivering analysis workflow: {config_model.analysis.analysis_workflow}")
    hk_file: Path = Path(output_dir, f"{config_model.analysis.case_id}.hk")
    delivery_ready_file: Path = Path(
        output_dir, f"{config_model.analysis.case_id}_delivery_ready.hk"
    )

    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--config",
        f"delivery=True",
        f"rules_to_deliver={','.join(rules_to_deliver)}",
        "--configfile",
        sample_config,
        "--cores",
        "1",
        "--dryrun",
        "--quiet",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.stdout.strip():
        LOG.info("snakemake stdout:\n%s", result.stdout)
    if result.stderr.strip():
        LOG.info("snakemake stderr:\n%s", result.stderr)

    if result.returncode != 0:
        LOG.error("Snakemake failed with code %s", result.returncode)
        raise SystemExit(result.returncode)

    hk_deliverables: List[Dict[str, Any]] = read_json(delivery_ready_file.as_posix())
    hk_deliverables: List[Dict[str, Any]] = convert_deliverables_tags(
        delivery_json=hk_deliverables, sample_config_dict=config
    )

    # Sample configuration file
    hk_deliverables.append(
        {
            "path": Path(sample_config).resolve().as_posix(),
            "step": "case_config",
            "format": get_file_extension(sample_config),
            "tag": ["balsamic-config"],
            "id": config_model.analysis.case_id,
        }
    )

    # DAG
    hk_deliverables.append(
        {
            "path": config_model.analysis.dag,
            "step": "case_config",
            "format": get_file_extension(config_model.analysis.dag),
            "tag": ["balsamic-dag"],
            "id": config_model.analysis.case_id,
        }
    )

    # MultiQC intermediate files
    multiqc_deliverables: List[Dict[str, Any]] = get_multiqc_deliverables(
        case_id=config_model.analysis.case_id,
        multiqc_dir=Path(config_model.analysis.result, "qc", "multiqc_data"),
    )
    hk_deliverables.extend(multiqc_deliverables)

    hk_deliverables: Dict[str, Any] = {"files": hk_deliverables}
    write_json(json_obj=hk_deliverables, path=hk_file.as_posix())
    write_yaml(data=hk_deliverables, file_path=f"{hk_file}.{FileType.YAML}")
    LOG.info(f"Generated analysis deliverables: {hk_file.as_posix()}")
