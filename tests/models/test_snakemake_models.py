"""Tests for the Snakemake model related methods."""
from pathlib import Path
from typing import Dict, Any

import pytest
from pydantic import ValidationError

from BALSAMIC.models.snakemake import SingularityBindPath, Snakemake


def test_singularity_bind_path_model(singularity_bind_path_data: Dict[str, Path]):
    """Test singularity bind path model initialisation."""

    # GIVEN singularity bind path model data

    # WHEN initialising the model
    singularity_bind_path_model: SingularityBindPath = SingularityBindPath(
        **singularity_bind_path_data
    )

    # THEN the model should have been correctly built
    assert singularity_bind_path_model.dict() == singularity_bind_path_data


def test_snakemake_model(
    snakemake_data: Dict[str, Any], snakemake_validated_data: Dict[str, Any]
):
    """Test snakemake model initialisation."""

    # GIVEN a cluster ready snakemake model

    # WHEN initialising the model
    snakemake_model: Snakemake = Snakemake(**snakemake_data)

    # THEN the model should have been correctly built
    assert snakemake_model.dict() == snakemake_validated_data


def test_snakemake_model_empty(
    snakemake_data: Dict[str, Any], snakemake_validated_data: Dict[str, Any]
):
    """Test Snakemake empty model initialisation."""

    # GIVEN no input for the snakemake model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        Snakemake()
