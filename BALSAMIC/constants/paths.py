"""Balsamic path constants."""
import sys
from pathlib import Path

# Balsamic working directory constants
BALSAMIC_DIR: Path = Path(sys.modules["BALSAMIC"].__file__).parent.resolve()
CONSTANTS_DIR: Path = Path(BALSAMIC_DIR, "constants")
CONTAINERS_DIR: Path = Path(BALSAMIC_DIR, "containers")
ASSETS_DIR: Path = Path(BALSAMIC_DIR, "assets")
SCRIPT_DIR: Path = Path(ASSETS_DIR, "scripts")
REFSEQ_SCRIPT_PATH: Path = Path(SCRIPT_DIR, "refseq_sql.awk")
IMMEDIATE_SUBMIT_PATH: Path = Path(SCRIPT_DIR, "immediate_submit.py")

# Sentieon specific constants
SENTIEON_MODELS_DIR: Path = Path(BALSAMIC_DIR, "assets", "sentieon_models")
SENTIEON_DNASCOPE_DIR: Path = Path(
    SENTIEON_MODELS_DIR, "SentieonDNAscopeModelBeta0.4a-201808.05.model"
)
SENTIEON_TNSCOPE_DIR: Path = Path(
    SENTIEON_MODELS_DIR, "SentieonTNscopeModel_GiAB_HighAF_LowFP-201711.05.model"
)

# Singularity container hardcoded paths
CADD_ANNOTATIONS_CONTAINER_DIR = Path("/opt/conda/share/CADD-scripts/data/annotations")

# Pytest paths
TEST_DATA_DIR: Path = Path("tests/test_data")
FASTQ_TEST_INFO: Path = Path(TEST_DATA_DIR, "fastq_test_info.json")
