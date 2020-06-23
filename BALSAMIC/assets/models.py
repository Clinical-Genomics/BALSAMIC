'''
Contains constants and models for analysis or filtering
'''
from pydantic import BaseModel
from typing import Optional

class NgsFilters(BaseModel):
    value: float
    name: str
    tag: str
    
class VarCallerFilter(BaseModel):
    AD: NgsFilters
    AF_min: Optional[NgsFilters]
    AF_max: Optional[NgsFilters]
    MQ: Optional[NgsFilters]
    DP: NgsFilters
    filter_version: str
    varcaller_name: str
    filter_type: str
    analysis_type: str
    description: str

VARDICT = VarCallerFilter(
    AD=NgsFilters(value=5, name="balsamic_low_tumor_ad", tag="INFO"),
    DP=NgsFilters(value=100, name="balsamic_low_tumor_dp", tag="INFO"),
    MQ=NgsFilters(value=50, name="balsamic_low_mq", tag="INFO"),
    AF_max=NgsFilters(value=1, name="balsamic_af_one", tag="INFO"),
    AF_min=NgsFilters(value=0.02, name="balsamic_low_af", tag="INFO"),
    filter_version = "0.0.1",
    varcaller_name = "VarDict",
    filter_type="general",
    analysis_type="tumor_only",
    description="General purpose filters used for filtering VarDict"
)
