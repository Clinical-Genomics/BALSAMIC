import pytest

from BALSAMIC.utils.models import (
    VCFAttributes, 
    VarCallerFilter, 
    QCModel, 
    VarcallerAttribute,
    VCFModel,
    AnalysisModel,
    SampleInstanceModel,
    BioinfoToolsModel
    )


def test_vcfattributes():
    """test VCFAttributes model for correct validation"""

    # GIVEN a VCF attribute
    dummy_attribute = {
        "tag_value": 5.0,
        "filter_name": "dummy_filter_name",
        "field": "INFO"
    }

    # WHEN building the model
    dummy_attribute_built = VCFAttributes(**dummy_attribute)

    # THEN assert values can be reterived currently
    assert dummy_attribute_built.tag_value == 5.0
    assert dummy_attribute_built.field == "INFO"
    assert dummy_attribute_built.filter_name == "dummy_filter_name"


def test_varcallerfilter():
    """test required VarCallerFilters for being set correctly"""

    # GIVEN a VarCallerFilter
    dummy_varcaller = {
        "AD": {
            "tag_value": 5.0,
            "filter_name": "dummy_alt_depth",
            "field": "INFO"
        },
        "DP": {
            "tag_value": 100.0,
            "filter_name": "dummy_depth",
            "field": "INFO"
        },
        "varcaller_name": "dummy_varcaller",
        "filter_type": "dummy_ffpe_filter",
        "analysis_type": "dummy_tumor_only",
        "description": "dummy description of this filter"
    }

    # WHEN building the model
    dummy_varcaller_filter = VarCallerFilter(**dummy_varcaller)

    # THEN assert required values are set
    assert dummy_varcaller_filter.AD.tag_value == 5.0
    assert dummy_varcaller_filter.DP.tag_value == 100.0
    assert dummy_varcaller_filter.analysis_type == "dummy_tumor_only"



def test_qc_model():
    #GIVEN valid input arguments
    #THEN we can successully create a config dict
    valid_args = {
        "umi_trim": True,
        "min_seq_length": 25,
        "umi_trim_length": 5}
    assert QCModel.parse_obj(valid_args)
    #GIVEN invalid input arguments
    invalid_args = {
        "umi_trim": True,
        "min_seq_length": "not number",
        "umi_trim_length": "but string"}
    #THEN should trigger ValueError
    with pytest.raises(ValueError) as excinfo:
        QCModel.parse_obj(invalid_args)
        assert "Invalid" in excinfo.value

def test_varcaller_attribute():
    #GIVEN valid input arguments
    valid_args={
            "mutation": "somatic",
            "type": "SNV"}
    #THEN we can successully create a config dict
    assert VarcallerAttribute.parse_obj(valid_args)
    #GIVEN invalid input arguments
    invalid_args={
            "mutation": "strange",
            "type": "unacceptable"}
    #THEN should trigger ValueError
    with pytest.raises(ValueError) as excinfo:
        VarcallerAttribute.parse_obj(invalid_args)
        assert "not a valid argument" in excinfo.value

def test_analysis_model():
    #GIVEN valid input arguments
    valid_args = {
        "case_id": "case_id",
        "analysis_type": "paired",
        "sequencing_type": "targeted",
        "analysis_dir": "tests/test_data"
        }
    #THEN we can successully create a config dict
    assert AnalysisModel.parse_obj(valid_args)

    #GIVEN invalid input arguments
    invalid_args = {
        "case_id": "case_id",
        "analysis_type": "odd",
        "sequencing_type": "wrong",
        "analysis_dir": "tests/test_data"
        }
    #THEN should trigger ValueError
    with pytest.raises(ValueError) as excinfo:
        AnalysisModel.parse_obj(invalid_args)
        assert "not supported" in excinfo.value

def test_sample_instance_model():
    #GIVEN valid input arguments
    valid_args = {
        "file_prefix": "S2_R",
        "type": "normal",
        }
    #THEN we can successully create a config dict
    assert SampleInstanceModel.parse_obj(valid_args)

    #GIVEN invalid input arguments
    invalid_args = {
        "file_prefix": "S2_R",
        "type": "fungal",
        }
    #THEN should trigger ValueError
    with pytest.raises(ValueError) as excinfo:
        SampleInstanceModel.parse_obj(invalid_args)
        assert "not supported" in excinfo.value
