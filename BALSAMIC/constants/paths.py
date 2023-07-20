"""Balsamic path constants."""
import sys
from pathlib import Path

# Balsamic working directory constants
BALSAMIC_DIR: Path = Path(sys.modules["BALSAMIC"].__file__).parent.resolve()
CONSTANTS_DIR: Path = Path(BALSAMIC_DIR, "constants")
CONTAINERS_DIR: Path = Path(BALSAMIC_DIR, "containers")
SCRIPT_DIR: Path = Path(BALSAMIC_DIR, "assets", "scripts")
REFSEQ_SCRIPT_PATH: Path = Path(SCRIPT_DIR, "refseq_sql.awk")
SCHEDULER_PATH: Path = Path(BALSAMIC_DIR, "utils", "scheduler.py")

# Sentieon specific constants
SENTIEON_MODELS_DIR: Path = Path(BALSAMIC_DIR, "assets", "sentieon_models")
SENTIEON_DNASCOPE_DIR: Path = Path(
    SENTIEON_MODELS_DIR, "SentieonDNAscopeModelBeta0.4a-201808.05.model"
)
SENTIEON_TNSCOPE_DIR: Path = Path(
    SENTIEON_MODELS_DIR, "SentieonTNscopeModel_GiAB_HighAF_LowFP-201711.05.model"
)
