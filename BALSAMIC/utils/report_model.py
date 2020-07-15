from pathlib import Path
from datetime import datetime

from pydantic import BaseModel, ValidationError, validator, Field
from pydantic.types import DirectoryPath, FilePath
from typing import Optional, List, Dict



class BalsamicOutputFile(BaseModel):
    """
    Model for each entry in deliverables.hk file.
    """
    path: str = Field(
        description="Path to file in string format"
        )
    date: str = Field(
        default=datetime.today().isoformat(), 
        description="Date of report"
        )
    step: str = Field(
        description="Step of the pipleine which generated the file"
        )
    file_format: str = Field(
        alias="format",
        description="Extension of the file"
        )
    tag : List[str] = Field(
        default=[],
        description="Additional tags to be associated with the file"
        )
    case_id: str =Field(
        alias="id",
        description="case_id associated with analysis"
        )
    status: str = Field(
        description="Completion status"
        )
