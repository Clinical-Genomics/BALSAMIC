import BALSAMIC
import sys

from pathlib import Path

CONDA_ENV_PATH = Path(
    sys.modules["BALSAMIC"].__file__).parent.resolve() / "conda"

CONDA_ENV_YAML = Path(sys.modules["BALSAMIC"].__file__).parent.resolve(
) / "config" / "balsamic_env.yaml"

RULE_DIRECTORY = Path(
    sys.modules["BALSAMIC"].__file__).parent.resolve().as_posix() + "/"


