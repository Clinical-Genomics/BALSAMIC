"""Utility functions for Balsamic analysis."""
import os
from pathlib import Path
from typing import List, Any, Dict

from BALSAMIC.constants.paths import ASSETS_DIR, BALSAMIC_DIR
from BALSAMIC.models.cache import CacheConfig
from BALSAMIC.models.snakemake import SingularityBindPath
from BALSAMIC.utils.cli import get_resolved_fastq_files_directory


def get_singularity_bind_paths(
    sample_config: Dict[str, Any]
) -> List[SingularityBindPath]:
    """Return a list of singularity binding paths for Balsamic analysis."""
    analysis_dir: Path = Path(sample_config["analysis"]["analysis_dir"])
    fastq_dir: Path = Path(
        get_resolved_fastq_files_directory(sample_config["analysis"]["fastq_path"])
    )
    cache_dir: Path = Path(os.path.commonpath(sample_config["reference"].values()))
    singularity_bind_paths: List[SingularityBindPath] = [
        SingularityBindPath(source=fastq_dir, destination=fastq_dir),
        SingularityBindPath(source=ASSETS_DIR, destination=ASSETS_DIR),
        SingularityBindPath(source=analysis_dir, destination=analysis_dir),
        SingularityBindPath(source=cache_dir, destination=cache_dir),
    ]
    if sample_config.get("panel"):
        capture_kit_path: Path = Path(sample_config.get("panel").get("capture_kit"))
        singularity_bind_paths.append(
            SingularityBindPath(source=capture_kit_path, destination=capture_kit_path)
        )
        if sample_config.get("panel").get("pon_cnn"):
            pon_cnn_path: Path = Path(sample_config.get("panel").get("pon_cnn"))
            singularity_bind_paths.append(
                SingularityBindPath(source=pon_cnn_path, destination=pon_cnn_path)
            )
    if sample_config.get("background_variants"):
        background_variants_path: Path = Path(sample_config.get("background_variants"))
        singularity_bind_paths.append(
            SingularityBindPath(
                source=background_variants_path, destination=background_variants_path
            )
        )
    if sample_config.get("reference").get("cadd_annotations"):
        cadd_annotations_path: Path = Path(sample_config.get("reference").get("cadd_annotations"))
        cadd_annotations_dest_path: Path = (
            "/opt/conda/share/CADD-scripts/data/annotations"
        )
        singularity_bind_paths.append(
            SingularityBindPath(
                source=cadd_annotations_path, destination=cadd_annotations_dest_path
            )
        )
    return singularity_bind_paths


def get_cache_singularity_bind_paths(
    cache_config: CacheConfig,
) -> List[SingularityBindPath]:
    """Return a list of singularity binding paths for Balsamic init."""
    return [
        SingularityBindPath(source=BALSAMIC_DIR, destination=BALSAMIC_DIR),
        SingularityBindPath(
            source=cache_config.references_dir, destination=cache_config.references_dir
        ),
    ]
