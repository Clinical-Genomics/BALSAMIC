'''
Contains constants and models for analysis or filtering
'''
from typing import Optional
from pydantic import BaseModel


class VCFAttributes(BaseModel):
    """General purpose filter to manage various VCF attributes

    This class handles three parameters for the purpose filtering variants
    based on a tag_values, filter_name, and which field in VCF.

    E.g. AD=VCFAttributes(tag_value=5, filter_name="balsamic_low_tumor_ad", field="INFO")
    A value of 5 from INFO field and filter_name will be balsamic_low_tumor_ad

    Attributes:
      tag_value: float
      filter_name: str
      field: str
    """

    tag_value: float
    filter_name: str
    field: str


class VarCallerFilter(BaseModel):
    """General purpose for variant caller filters

    This class handles attributes and filter for variant callers

    Attributes:
      AD: VCFAttributes (required); minimum allelic depth
      AF_min: VCFAttributes (optional); minimum allelic fraction
      AF_max: VCFAttributes (optional); maximum allelic fraction
      MQ: VCFAttributes (optional); minimum mapping quality
      DP: VCFAttributes (optional); minimum read depth
      varcaller_name: str (required); variant caller name
      filter_type: str (required); filter name for variant caller
      analysis_type: str (required); analysis type e.g. tumor_normal or tumor_only
      description: str (required); comment section for description
    """

    AD: VCFAttributes
    AF_min: Optional[VCFAttributes]
    AF_max: Optional[VCFAttributes]
    MQ: Optional[VCFAttributes]
    DP: VCFAttributes
    varcaller_name: str
    filter_type: str
    analysis_type: str
    description: str


VARDICT = VarCallerFilter(
    AD=VCFAttributes(tag_value=5,
                     filter_name="balsamic_low_tumor_ad",
                     field="INFO"),
    DP=VCFAttributes(tag_value=100,
                     filter_name="balsamic_low_tumor_dp",
                     field="INFO"),
    MQ=VCFAttributes(tag_value=50, filter_name="balsamic_low_mq",
                     field="INFO"),
    AF_max=VCFAttributes(tag_value=1,
                         filter_name="balsamic_af_one",
                         field="INFO"),
    AF_min=VCFAttributes(tag_value=0.02,
                         filter_name="balsamic_low_af",
                         field="INFO"),
    varcaller_name="VarDict",
    filter_type="general",
    analysis_type="tumor_only",
    description="General purpose filters used for filtering VarDict")
