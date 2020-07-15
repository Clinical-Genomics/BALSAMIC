import BALSAMIC
import sys

from pathlib import Path
"""This file contains consatant variables used during BALSAMIC config generation"""
"""Path to conda folder containing YAML files with verions of software usen un BALSAMIC workflow"""
CONDA_ENV_PATH = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() /
    "conda").as_posix()
"""Path to config YAML file to be accessed by Snakemake"""

CONDA_ENV_YAML = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" /
    "balsamic_env.yaml").as_posix()
""""Path to rule files to be accessed by Snakemake"""
RULE_DIRECTORY = Path(
    sys.modules["BALSAMIC"].__file__).parent.resolve().as_posix() + "/"
