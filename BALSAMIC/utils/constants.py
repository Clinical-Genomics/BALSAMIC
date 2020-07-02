import BALSAMIC
import sys

from pathlib import Path

CONDA_ENV_PATH = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() /
    "conda").as_posix()

CONDA_ENV_YAML = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" /
    "balsamic_env.yaml").as_posix()

RULE_DIRECTORY = Path(
    sys.modules["BALSAMIC"].__file__).parent.resolve().as_posix() + "/"
