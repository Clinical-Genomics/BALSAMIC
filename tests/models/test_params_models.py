"""Tests for Balsamic analysis params models."""
from math import isclose

import pytest
from pydantic import ValidationError

from BALSAMIC.constants.analysis import SequencingType
from BALSAMIC.constants.workflow_params import WORKFLOW_PARAMS
from BALSAMIC.models.params import BalsamicWorkflowConfig
from BALSAMIC.models.config import VarcallerAttribute
from BALSAMIC.models.params import (
    ParamsManta,
    ParamsSentieonWGSMetrics,
    ParamsVEP,
    QCModel,
    UMIParamsCommon,
    UMIParamsConsensuscall,
    ParamsTNscope,
    UMIParamsUMIextract,
    VCFFilter,
    StructuralVariantFilters,
)


def test_params_sentieon_wgs_metrics():
    """Test sentieon wgs metrics settings model for correct validation."""

    # GIVEN Manta params
    test_sentieon_wgs_metrics_params = {
        "min_base_qual": 10,
        "cov_threshold": [50, 100, 150, 200, 250],
    }

    # WHEN building the model
    test_sentieon_wgs_metrics_built = ParamsSentieonWGSMetrics(
        **test_sentieon_wgs_metrics_params
    )

    # THEN values should be correctly populated and parsed into the model
    assert test_sentieon_wgs_metrics_built.min_base_qual == 10
    assert (
        test_sentieon_wgs_metrics_built.cov_threshold
        == "--cov_thresh 50 --cov_thresh 100 --cov_thresh 150 --cov_thresh 200 --cov_thresh 250"
    )


def test_params_manta():
    """Test Manta settings model for correct validation."""

    # GIVEN Manta params
    test_manta_params = {"wgs_settings": "", "tga_settings": "--exome"}

    # WHEN building the model
    test_manta_built = ParamsManta(**test_manta_params)

    # THEN string values should be correctly populated into the model
    assert test_manta_built.tga_settings == "--exome"
    assert test_manta_built.wgs_settings == ""


def test_get_manta_settings_tga():
    """Test get Manta settings based on sequencing type TGA."""

    # GIVEN workflow params
    params = BalsamicWorkflowConfig.model_validate(WORKFLOW_PARAMS)

    # WHEN getting manta settings for TGA
    manta_settings = params.get_manta_settings(SequencingType.TARGETED)

    # THEN manta setting should specify exome argument
    assert manta_settings == "--exome"


def test_get_manta_settings_wgs():
    """Test get Manta settings based on sequencing type WGS."""
    # GIVEN workflow params
    params = BalsamicWorkflowConfig.model_validate(WORKFLOW_PARAMS)

    # WHEN getting manta settings for WGS
    manta_settings = params.get_manta_settings(SequencingType.WGS)

    # THEN manta setting should be empty
    assert manta_settings == ""


def test_params_vep():
    """test UMIParamsVEP model for correct validation"""

    # GIVEN params
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


def test_varcallerfilter():
    """test VCFFilter model for correct validation"""

    # GIVEN a VCF attribute
    dummy_attribute = {
        "tag_value": 5.0,
        "filter_name": "dummy_filter_name",
        "field": "INFO",
    }

    # WHEN building the model
    dummy_attribute_built = VCFFilter(**dummy_attribute)

    # THEN assert required values are set
    assert isclose(dummy_attribute_built.tag_value, 5.0)
    assert dummy_attribute_built.field == "INFO"
    assert dummy_attribute_built.filter_name == "dummy_filter_name"


def test_structuralvariantfilters():
    """test StructuralVariantFilters model for correct validation"""

    # GIVEN a SV VarCallerFilter
    dummy_varcaller = {
        "low_pr_sr_count": {
            "tag_value": 4,
            "filter_name": "low_pr_sr_count",
            "field": "INFO",
        },
        "loqusdb_clinical_sv_freq": {
            "tag_value": 0.02,
            "filter_name": "dummy_pop_freq",
            "field": "INFO",
        },
        "varcaller_name": "dummy_varcaller",
        "filter_type": "dummy_ffpe_filter",
        "analysis_type": "dummy_tumor_only",
        "description": "dummy description of this filter",
    }

    # WHEN building the model
    dummy_varcaller_filter = StructuralVariantFilters(**dummy_varcaller)

    # THEN assert required values are set
    assert isclose(dummy_varcaller_filter.low_pr_sr_count.tag_value, 4)
    assert isclose(dummy_varcaller_filter.loqusdb_clinical_sv_freq.tag_value, 0.02)
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
        "align_intbases": 100,
    }
    # WHEN building the model
    test_commonparams_built = UMIParamsCommon(**test_commonparams)
    # THEN assert values
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


def test_params_tnscope():
    """test ParamsTNscope model for correct validation"""

    # GIVEN tnscope params
    test_tnscope_params = {
        "algo": "algoname",
        "init_tumorLOD": 0.5,
        "min_tumorLOD": 6,
        "filter_tumor_af": 0.01,
        "error_rate": 5,
        "prunefactor": 3,
        "padding": 30,
        "disable_detect": "abc",
        "pcr_model": "NONE",
    }

    # WHEN building the model
    test_tnscope_params_built = ParamsTNscope(**test_tnscope_params)

    # THEN assert values
    assert test_tnscope_params_built.algo == "algoname"
    assert isclose(test_tnscope_params_built.init_tumorLOD, 0.5)
    assert test_tnscope_params_built.min_tumorLOD == 6
    assert isclose(test_tnscope_params_built.filter_tumor_af, 0.01, rel_tol=1e-9)
    assert test_tnscope_params_built.error_rate == 5
    assert test_tnscope_params_built.prunefactor == 3
    assert test_tnscope_params_built.disable_detect == "abc"
    assert test_tnscope_params_built.pcr_model == "NONE"
    assert test_tnscope_params_built.padding == 30
