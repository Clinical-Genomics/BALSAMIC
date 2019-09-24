#!/usr/bin/env python
import os
import click
import json
import logging
from pathlib import Path

LOG = logging.getLogger(__name__)

@click.command("init", short_help="Run this command after installing BALSAMIC")
@click.option("--singularity",
              type=click.Path(),
              required=True,
              help='Download singularity image for BALSAMIC')
@click.pass_context
def initiate(context, singularity):
    """
    Creates install.json, which is required for workflow to function properly 
    
    """

    config_path = Path(__file__).parents[2] / "config"
    config_path = config_path.absolute()

    install_json_file = config_path / "install.json"
    balsamic_env = config_path / "balsamic_env.yaml"
    rule_directory = Path(__file__).parents[2]

    install_json = dict()
     
    install_json["conda_env_yaml"] = balsamic_env.as_posix()
    install_json["rule_directory"] = rule_directory.as_posix() + "/"

    install_json["singularity"] = dict()
    install_json["singularity"]["image"] = Path(singularity).absolute().as_posix()

    LOG.info("Writing install_json file to: %s", install_json_file.as_posix())
    with open(install_json_file, 'w') as json_out:
        json.dump(install_json, json_out, indent=4)
