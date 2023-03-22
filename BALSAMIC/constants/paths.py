"""Common and general balsamic constants."""
import sys
from pathlib import Path

BALSAMIC_DIR: str = Path(sys.modules["BALSAMIC"].__file__).parent.resolve().as_posix()

SCRIPT_DIR: str = Path(BALSAMIC_DIR, "assets", "scripts").as_posix()

CONTAINERS_DIR: str = Path(BALSAMIC_DIR, "containers").as_posix()

SENTIEON_MODELS_DIR: str = Path(BALSAMIC_DIR, "assets", "sentieon_models").as_posix()

SENTIEON_DNASCOPE_DIR: str = Path(
    SENTIEON_MODELS_DIR, "SentieonDNAscopeModelBeta0.4a-201808.05.model"
).as_posix()

SENTIEON_TNSCOPE_DIR: str = Path(
    SENTIEON_MODELS_DIR, "SentieonTNscopeModel_GiAB_HighAF_LowFP-201711.05.model"
).as_posix()
