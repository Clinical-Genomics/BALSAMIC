"""Balsamic cache models."""
from typing import Dict, Optional

from pydantic import BaseModel, validator

from BALSAMIC.constants.analysis import GenomeVersion
from BALSAMIC.models.validators import validate_genome_version, validate_dir_path


class SingularityContainersModel(BaseModel):
    """
    Singularity containers model.

    Attributes:
        image_dir: path to the singularity images
        containers: dictionary linking the container names and their dockerhub image paths
    """

    image_dir: str
    containers: Dict[str, str]

    _image_dir = validator("image_dir", allow_reuse=True)(validate_dir_path)


class CacheConfigModel(BaseModel):
    """
    Reference build configuration model.

    Attributes:
        output_dir: output directory for the generated cache
        balsamic_dir: project directory of balsamic
        genome_version: genome version associated with the balsamic cache
        cosmic_key: COSMIC database key
        bioinfo_tools: dictionary of bioinformatics software and their associated containers
        singularity: singularity cache model

    """

    output_dir: str
    balsamic_dir: str
    genome_version: GenomeVersion
    cosmic_key: Optional[str]
    bioinfo_tools: dict
    singularity: SingularityContainersModel

    _dir = validator("output_dir", "balsamic_dir", allow_reuse=True)(validate_dir_path)
    _genome_version = validator("genome_version", allow_reuse=True)(
        validate_genome_version
    )
