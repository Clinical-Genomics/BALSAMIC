"""Tests for Balsamic analysis params models."""
import pytest
from pydantic import ValidationError

from BALSAMIC.models.config import VarcallerAttribute
from BALSAMIC.models.params import (
    ParamsVardict,
    ParamsVEP,
    QCModel,
    UMIParamsCommon,
    UMIParamsConsensuscall,
    UMIParamsTNscope,
    UMIParamsUMIextract,
    VarCallerFilter,
    VCFAttributes,
)


def test_params_vardict():
    """test UMIParamsVardict model for correct validation"""

    # GIVEN vardict params
    test_vardict_params = {
        "allelic_frequency": 0.01,
        "max_pval": 0.5,
        "max_mm": 2,
        "column_info": "-a 1 -b 2 -c 3",
    }

    # WHEN building the model
    test_vardict_built = ParamsVardict(**test_vardict_params)

    # THEN assert values
    assert test_vardict_built.allelic_frequency == 0.01
    assert test_vardict_built.max_pval == 0.5
    assert test_vardict_built.max_mm == 2
    assert test_vardict_built.column_info == "-a 1 -b 2 -c 3"


def test_params_vep():
    """test UMIParamsVEP model for correct validation"""

    # GIVEN vardict params
    test_vep = {"vep_filters": "all defaults params"}

    # WHEN building the model
    test_vep_built = ParamsVEP(**test_vep)

    # THEN assert values
    assert test_vep_built.vep_filters == "all defaults params"


def test_qc_model():
    # GIVEN valid input arguments
    # THEN we can successully create a config dict
    valid_args = {
        "umi_trim": True,
        "min_seq_length": 25,
        "umi_trim_length": 5,
        "n_base_limit": 50,
    }
    assert QCModel.model_validate(valid_args)


def test_vcfattributes():
    """test VCFAttributes model for correct validation"""

    # GIVEN a VCF attribute
    dummy_attribute = {
        "tag_value": 5.0,
        "filter_name": "dummy_filter_name",
        "field": "INFO",
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
        "AD": {"tag_value": 5.0, "filter_name": "dummy_alt_depth", "field": "INFO"},
        "DP": {"tag_value": 100.0, "filter_name": "dummy_depth", "field": "INFO"},
        "pop_freq": {
            "tag_value": 0.005,
            "filter_name": "dummy_pop_freq",
            "field": "INFO",
        },
        "varcaller_name": "dummy_varcaller",
        "filter_type": "dummy_ffpe_filter",
        "analysis_type": "dummy_tumor_only",
        "description": "dummy description of this filter",
    }

    # WHEN building the model
    dummy_varcaller_filter = VarCallerFilter(**dummy_varcaller)

    # THEN assert required values are set
    assert dummy_varcaller_filter.AD.tag_value == 5.0
    assert dummy_varcaller_filter.DP.tag_value == 100.0
    assert dummy_varcaller_filter.analysis_type == "dummy_tumor_only"


def test_varcaller_attribute():
    # GIVEN valid input arguments
    valid_args = {"mutation": "somatic", "mutation_type": "SNV"}
    # THEN we can successully create a config dict
    assert VarcallerAttribute.model_validate(valid_args)
    # GIVEN invalid input arguments
    invalid_args = {"mutation": "strange", "mutation_type": "unacceptable"}
    # THEN should trigger ValueError
    with pytest.raises(ValidationError) as excinfo:
        VarcallerAttribute.model_validate(invalid_args)
    assert "2 validation errors" in str(excinfo.value)


def test_umiparams_common():
    """test UMIParamsCommon model for correct validation"""

    # GIVEN a UMI workflow common params
    test_commonparams = {
        "align_header": "test_header_name",
        "align_intbases": 100,
        "filter_tumor_af": 0.01,
    }
    # WHEN building the model
    test_commonparams_built = UMIParamsCommon(**test_commonparams)
    # THEN assert values
    assert test_commonparams_built.align_header == "test_header_name"
    assert test_commonparams_built.filter_tumor_af == 0.01
    assert test_commonparams_built.align_intbases == 100


def test_umiparams_umiextract():
    """test UMIParamsUMIextract model for correct validation"""
    # GIVEN umiextract params
    test_umiextractparams = {"read_structure": "['mode', 'r1,r2']"}

    # WHEN building the model
    test_umiextractparams_built = UMIParamsUMIextract(**test_umiextractparams)

    # THEN assert values
    assert test_umiextractparams_built.read_structure == "['mode', 'r1,r2']"


def test_umiparams_consensuscall():
    """test UMIParamsConsensuscall model for correct validation"""

    # GIVEN consensuscall params
    test_consensuscall = {
        "align_format": "BAM",
        "filter_minreads": "6,3,3",
        "tag": "XZ",
    }

    # WHEN building the model
    test_consensuscall_built = UMIParamsConsensuscall(**test_consensuscall)

    # THEN assert values
    assert test_consensuscall_built.align_format == "BAM"
    assert test_consensuscall_built.filter_minreads == "6,3,3"
    assert test_consensuscall_built.tag == "XZ"


def test_umiparams_tnscope():
    """test UMIParamsTNscope model for correct validation"""

    # GIVEN tnscope params
    test_tnscope_params = {
        "algo": "algoname",
        "init_tumorLOD": 0.5,
        "min_tumorLOD": 6,
        "error_rate": 5,
        "prunefactor": 3,
        "padding": 30,
        "disable_detect": "abc",
    }

    # WHEN building the model
    test_tnscope_params_built = UMIParamsTNscope(**test_tnscope_params)

    # THEN assert values
    assert test_tnscope_params_built.algo == "algoname"
    assert test_tnscope_params_built.init_tumorLOD == 0.5
    assert test_tnscope_params_built.min_tumorLOD == 6
    assert test_tnscope_params_built.error_rate == 5
    assert test_tnscope_params_built.prunefactor == 3
    assert test_tnscope_params_built.disable_detect == "abc"
    assert test_tnscope_params_built.padding == 30
