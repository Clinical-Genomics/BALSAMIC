"""Test module for Balsamic common models."""
import pytest

from BALSAMIC.models.analysis import (
    VCFAttributes,
    VarCallerFilter,
    QCModel,
    VarcallerAttribute,
    AnalysisModel,
    SampleInstanceModel,
    UMIParamsCommon,
    UMIParamsUMIextract,
    UMIParamsConsensuscall,
    UMIParamsTNscope,
    ParamsVardict,
    ParamsVEP,
    AnalysisPonModel,
    FastqInfoModel,
    BalsamicConfigModel,
)
from pydantic import ValidationError
from BALSAMIC.utils.exc import BalsamicError
from typing import List, Dict
from BALSAMIC.constants.analysis import FastqName, SampleType
import os
import copy
from pathlib import Path


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


def test_qc_model():
    # GIVEN valid input arguments
    # THEN we can successully create a config dict
    valid_args = {
        "umi_trim": True,
        "min_seq_length": 25,
        "umi_trim_length": 5,
        "n_base_limit": 50,
    }
    assert QCModel.parse_obj(valid_args)


def test_varcaller_attribute():
    # GIVEN valid input arguments
    valid_args = {"mutation": "somatic", "type": "SNV"}
    # THEN we can successully create a config dict
    assert VarcallerAttribute.parse_obj(valid_args)
    # GIVEN invalid input arguments
    invalid_args = {"mutation": "strange", "type": "unacceptable"}
    # THEN should trigger ValueError
    with pytest.raises(ValidationError) as excinfo:
        VarcallerAttribute.parse_obj(invalid_args)
    assert ("strange is not a valid mutation type" in str(excinfo.value))
    assert ("unacceptable is not not a valid mutation class" in str(excinfo.value))


def test_analysis_model(test_data_dir: str):
    """Test analysis model instantiation."""

    # GIVEN valid input arguments
    valid_args = {
        "case_id": "case_id",
        "gender": "female",
        "analysis_type": "paired",
        "sequencing_type": "targeted",
        "analysis_dir": test_data_dir,
        "fastq_path": test_data_dir,
        "analysis_workflow": "balsamic-umi",
    }

    # THEN we can successfully create a config dict
    assert AnalysisModel.parse_obj(valid_args)

    # GIVEN invalid input arguments
    invalid_args = {
        "case_id": "case_id",
        "gender": "unknown",
        "analysis_type": "odd",
        "sequencing_type": "wrong",
        "analysis_dir": "tests/test_data",
        "analysis_workflow": "umi",
    }

    # THEN should trigger ValueError
    with pytest.raises(ValueError) as excinfo:
        AnalysisModel.parse_obj(invalid_args)
    assert "not supported" in str(excinfo.value)


def test_sample_instance_model(config_dict_w_fastqs):
    """Test sample instance model initialisation."""

    # GIVEN a sample list
    sample_list = config_dict_w_fastqs["samples"]

    # WHEN parsing the sample dictionary
    for idx, sample in enumerate(sample_list):
        sample: SampleInstanceModel = SampleInstanceModel.parse_obj(sample)

        sample_dict_copy = sample_list[idx].copy()
        for fastq_pattern, values in sample_dict_copy["fastq_info"].items():
            values["fwd"] = Path(values["fwd"]).resolve().as_posix()
            values["rev"] = Path(values["rev"]).resolve().as_posix()

        # THEN the sample model should be correctly initialised
        assert sample.dict() == sample_dict_copy


def test_sample_instance_model_sample_type_error(tumor_normal_fastq_info_correct):
    """Test sample instance model error raise."""

    # GIVEN a sample dictionary with an invalid sample type
    samples: List[Dict] = copy.deepcopy(tumor_normal_fastq_info_correct)
    illegal_sample_type: str = "affected"
    tumor_dict = samples[0]
    tumor_dict["type"] = illegal_sample_type

    # WHEN parsing the sample dictionary
    # THEN a ValueError should be triggered
    with pytest.raises(ValueError) as exc:
        SampleInstanceModel.parse_obj(tumor_dict)
        assert (
            f"The provided sample type ({illegal_sample_type}) is not supported in BALSAMIC"
            in exc.value
        )


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


def test_analysis_pon_model(test_data_dir: str):
    """Tests PON model parsing."""

    # GIVEN valid input arguments
    valid_args = {
        "case_id": "case_id",
        "analysis_type": "pon",
        "sequencing_type": "targeted",
        "analysis_dir": test_data_dir,
        "fastq_path": test_data_dir,
        "analysis_workflow": "balsamic",
        "pon_version": "v1",
    }

    # THEN we can successfully create a config dict
    assert AnalysisPonModel.parse_obj(valid_args)

    # GIVEN an invalid version argument
    invalid_args = {
        "case_id": "case_id",
        "analysis_type": "pon",
        "sequencing_type": "targeted",
        "analysis_dir": test_data_dir,
        "fastq_path": test_data_dir,
        "analysis_workflow": "balsamic",
        "pon_version": "v01",
    }

    # THEN should trigger ValueError
    with pytest.raises(ValidationError) as excinfo:
        AnalysisPonModel.parse_obj(invalid_args)

    assert (
        f"The provided version ({invalid_args['pon_version']}) does not follow the defined syntax (v<int>)"
        in str(excinfo.value)
    )


def test_illegal_sample_name(
    tumor_normal_fastq_info_correct, illegal_normal_sample_name
):
    """Test sample instance model detection of illegal sample name containing underscore."""

    # GIVEN a sample dictionary with an invalid sample name
    modify_sample_dict = copy.deepcopy(tumor_normal_fastq_info_correct)
    samples: List[Dict] = modify_sample_dict
    normal_dict = samples[1]
    normal_dict["name"] = illegal_normal_sample_name

    # WHEN parsing the sample dictionary
    # THEN a ValueError should be triggered
    with pytest.raises(ValueError) as exc:
        SampleInstanceModel.parse_obj(normal_dict)

    assert (
        f"Sample name '{illegal_normal_sample_name}' contains an underscore (_). Underscores are not allowed."
        in str(exc.value)
    )


def test_detect_duplicate_fastq_pattern(
    config_w_fastq_dir_for_duplicate_fastq_patterns_model: Dict,
):
    """Test balsamic models ability to detect duplicate assigned fastq patterns."""
    config_dict = config_w_fastq_dir_for_duplicate_fastq_patterns_model
    # Initialize balsamic model
    with pytest.raises(ValueError) as exc:
        BalsamicConfigModel.parse_obj(config_dict)

    assert (
        f"Duplicate FastqPattern(s) found: ACC1_S1_L001_R across multiple samples"
        in str(exc.value)
    )


def test_detection_unassigned_fastq_file(config_tumor_normal_extrafile: Dict):
    """Test instantiating balsamic model with fastq dir containing unassigned fastq-files."""
    # Initialize balsamic model
    with pytest.raises(ValueError) as exc:
        BalsamicConfigModel.parse_obj(config_tumor_normal_extrafile)

    assert f"Fastqs in fastq-dir not assigned to sample config:" in str(exc.value)


def test_get_all_sample_names(balsamic_model):
    """Validate retrieval of all sample names in analysis from BalsamicConfigModel."""
    sample_names = balsamic_model.get_all_sample_names()
    assert ["ACC1", "ACC2"] == sample_names


def test_get_fastq_patterns_by_sample(
    balsamic_model, tumor_sample_name, normal_sample_name
):
    """Validate retrieval of fastq-pattern by sample from BalsamicConfigModel."""

    def compare_fastq_pattern_lists(expected: List[str], found: List[str]):
        assert all(
            fastq_pattern in found for fastq_pattern in expected
        ), "Not all expected fastq patterns found."
        assert len(expected) == len(found), "Not same number of fastq patterns"

    tumor_fastq_patterns_expected = [
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX",
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX",
    ]
    normal_fastq_patterns_expected = [
        "1_171015_HJ7TLDSX5_ACC2_XXXXXX",
        "2_171015_HJ7TLDSX5_ACC2_XXXXXX",
    ]

    fastq_patterns_all_expected = (
        tumor_fastq_patterns_expected + normal_fastq_patterns_expected
    )

    tumor_fastq_patterns = balsamic_model.get_fastq_patterns_by_sample(
        [tumor_sample_name]
    )
    normal_fastq_patterns = balsamic_model.get_fastq_patterns_by_sample(
        [normal_sample_name]
    )
    fastq_patterns_all = balsamic_model.get_fastq_patterns_by_sample(
        [tumor_sample_name, normal_sample_name]
    )

    compare_fastq_pattern_lists(tumor_fastq_patterns_expected, tumor_fastq_patterns)
    compare_fastq_pattern_lists(normal_fastq_patterns_expected, normal_fastq_patterns)
    compare_fastq_pattern_lists(fastq_patterns_all_expected, fastq_patterns_all)


def test_get_all_fastqs_for_sample(balsamic_model, tumor_sample_name):
    """Validate retrieval of fastq-files by sample and fastq-type from BalsamicConfigModel."""

    def compare_fastq_file_lists(expected: List[str], found: List[str]):
        found_file_names = []
        for found_file in found:
            found_file_names.append(os.path.basename(found_file))
        assert all(
            fastq_file in found_file_names for fastq_file in expected
        ), f"Not all expected fastq files found. {expected}: {found_file_names}"
        assert len(expected) == len(found), "Not same number of fastq files"

    fwd_fastq_files_expected = [
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX_1.fastq.gz",
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX_1.fastq.gz",
    ]
    rev_fastq_files_expected = [
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX_2.fastq.gz",
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX_2.fastq.gz",
    ]
    fastq_files_expected = fwd_fastq_files_expected + rev_fastq_files_expected

    normal_fastq = "1_171015_HJ7TLDSX5_ACC2_XXXXXX_1.fastq.gz"

    fwd_fastq_files = balsamic_model.get_all_fastqs_for_sample(
        tumor_sample_name, [FastqName.FWD]
    )
    rev_fastq_files = balsamic_model.get_all_fastqs_for_sample(
        tumor_sample_name, [FastqName.REV]
    )
    fastq_files = balsamic_model.get_all_fastqs_for_sample(
        tumor_sample_name, [FastqName.FWD, FastqName.REV]
    )

    compare_fastq_file_lists(fwd_fastq_files_expected, fwd_fastq_files)
    compare_fastq_file_lists(rev_fastq_files_expected, rev_fastq_files)
    compare_fastq_file_lists(fastq_files_expected, fastq_files)
    assert normal_fastq not in fastq_files


def test_get_all_fastq_names(balsamic_model):
    """Validate retrieval of all fastq-files from BalsamicConfigModel and optional removal of suffix."""

    all_fastqs_expected = [
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX_1.fastq.gz",
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX_2.fastq.gz",
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX_1.fastq.gz",
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX_2.fastq.gz",
        "1_171015_HJ7TLDSX5_ACC2_XXXXXX_1.fastq.gz",
        "1_171015_HJ7TLDSX5_ACC2_XXXXXX_2.fastq.gz",
        "2_171015_HJ7TLDSX5_ACC2_XXXXXX_1.fastq.gz",
        "2_171015_HJ7TLDSX5_ACC2_XXXXXX_2.fastq.gz",
    ]
    all_fastqs_expected_wo_suffix = [
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX_1",
        "1_171015_HJ7TLDSX5_ACC1_XXXXXX_2",
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX_1",
        "2_171015_HJ7TLDSX5_ACC1_XXXXXX_2",
        "1_171015_HJ7TLDSX5_ACC2_XXXXXX_1",
        "1_171015_HJ7TLDSX5_ACC2_XXXXXX_2",
        "2_171015_HJ7TLDSX5_ACC2_XXXXXX_1",
        "2_171015_HJ7TLDSX5_ACC2_XXXXXX_2",
    ]

    all_fastqs = balsamic_model.get_all_fastq_names()
    all_fastqs_wo_suffix = balsamic_model.get_all_fastq_names(remove_suffix=True)

    assert all_fastqs == all_fastqs_expected
    assert all_fastqs_wo_suffix == all_fastqs_expected_wo_suffix


def test_fastq_by_fastq_pattern(balsamic_model):
    """Validate retrieval of fastq-file by fastq-pattern and fastq-type from BalsamicConfigModel."""

    fastq_pattern = "2_171015_HJ7TLDSX5_ACC2_XXXXXX"
    expected_fwd = "2_171015_HJ7TLDSX5_ACC2_XXXXXX_1.fastq.gz"
    expected_rev = "2_171015_HJ7TLDSX5_ACC2_XXXXXX_2.fastq.gz"

    fwd_fastq = balsamic_model.get_fastq_by_fastq_pattern(fastq_pattern, FastqName.FWD)
    rev_fastq = balsamic_model.get_fastq_by_fastq_pattern(fastq_pattern, FastqName.REV)

    assert os.path.basename(fwd_fastq) == expected_fwd
    assert os.path.basename(rev_fastq) == expected_rev


def test_illegal_fastqtype_get_fastq_by_fastq_pattern(balsamic_model):
    """Validate raise exception for illegal fastqtype BalsamicConfigModel."""

    fastq_pattern = "2_171015_HJ7TLDSX5_ACC2_XXXXXX"

    with pytest.raises(ValueError) as excinfo:
        balsamic_model.get_fastq_by_fastq_pattern(fastq_pattern, "forward")
    assert (
        f"fastq_type must be either {FastqName.FWD} or {FastqName.REV}, not: forward"
        in str(excinfo.value)
    )


def test_sample_name_by_type(balsamic_model):
    """Validate retrieval of sample name by sample type from BalsamicConfigModel."""

    tumor_name_expected = "ACC1"
    normal_name_expected = "ACC2"

    # Given sample type
    tumor_name_retrieved = balsamic_model.get_sample_name_by_type(SampleType.TUMOR)
    normal_name_retrieved = balsamic_model.get_sample_name_by_type(SampleType.NORMAL)

    # Then the retrieved sample name should match the expected
    assert tumor_name_retrieved == tumor_name_expected
    assert normal_name_retrieved == normal_name_expected


def test_sample_type_by_name(balsamic_model):
    """Validate retrieval of sample type by sample name from BalsamicConfigModel."""

    tumor_name = "ACC1"
    normal_name = "ACC2"

    # Given sample name
    tumor_type_retrieved = balsamic_model.get_sample_type_by_name(tumor_name)
    normal_type_retrieved = balsamic_model.get_sample_type_by_name(normal_name)

    # Then the retrieved sample type should match the expected
    assert tumor_type_retrieved == SampleType.TUMOR
    assert normal_type_retrieved == SampleType.NORMAL


def test_get_bam_name_per_lane(balsamic_model):
    """Validate retrieval of per lane bam names by sample name."""

    def compare_bam_file_lists(expected: List[str], found: List[str]):
        assert all(
            bam_file in found for bam_file in expected
        ), "Not all expected bam files found."
        assert len(expected) == len(found), "Not same number of bam files"

    # Fastq patterns for ACC2 in config.json
    normal_lane1_fastq_pattern = "1_171015_HJ7TLDSX5_ACC2_XXXXXX"
    normal_lane2_fastq_pattern = "2_171015_HJ7TLDSX5_ACC2_XXXXXX"

    # Given bam_dir path and sample name
    normal_name = "ACC2"
    result_dir = balsamic_model.analysis.result
    bam_dir = os.path.join(result_dir, "bam", "")

    # When retrieving all per lane bam names for sample
    bam_names = balsamic_model.get_bam_name_per_lane(bam_dir, normal_name)

    # Then the bam names for all fastq patterns should be retrieved and match the expected format
    expected_bam_name_lane1 = (
        f"{bam_dir}{normal_name}_align_sort_{normal_lane1_fastq_pattern}.bam"
    )
    expected_bam_name_lane2 = (
        f"{bam_dir}{normal_name}_align_sort_{normal_lane2_fastq_pattern}.bam"
    )
    expected_bam_names = [expected_bam_name_lane1, expected_bam_name_lane2]
    compare_bam_file_lists(expected_bam_names, bam_names)


def test_get_final_bam_name(balsamic_model):
    """Validate retrieval of final bam name by either sample type or sample name."""

    # Given bam_dir path and sample name or sample type
    sample_name = "ACC1"
    sample_type = SampleType.TUMOR
    result_dir = balsamic_model.analysis.result
    bam_dir = os.path.join(result_dir, "bam", "")

    # When retrieving final bam file name by sample name or sample type
    bam_name_sample_name = balsamic_model.get_final_bam_name(
        bam_dir, sample_name=sample_name
    )
    bam_name_sample_type = balsamic_model.get_final_bam_name(
        bam_dir, sample_type=sample_type
    )

    # Then retrieved final bam names should match the expected format and be identical regardless of request parameter
    expected_final_bam_name = f"{bam_dir}{sample_type}.{sample_name}.dedup.realign.bam"
    assert expected_final_bam_name == bam_name_sample_name
    assert bam_name_sample_name == bam_name_sample_type


def test_no_info_error_get_final_bam_name(balsamic_model):
    """Validate raise ValueError by not sample type or sample name."""

    # Given bam_dir path
    result_dir = balsamic_model.analysis.result
    bam_dir = os.path.join(result_dir, "bam", "")

    # When retrieving final bam file name without supplying sample name or sample type
    # ValueError should be raised
    with pytest.raises(ValueError) as excinfo:
        balsamic_model.get_final_bam_name(bam_dir)
    assert (
        f"Either sample_name or sample_type must be provided to get the final bam name."
        in str(excinfo.value)
    )


def test_get_final_bam_name_pon(balsamic_pon_model):
    """Validate retrieval of final bam name for PON by either sample type or sample name."""

    # Given bam_dir path and sample name or sample type
    sample_name = "ACCN6"
    sample_type = SampleType.NORMAL
    result_dir = balsamic_pon_model.analysis.result
    bam_dir = os.path.join(result_dir, "bam", "")

    # When retrieving final bam file name by sample name or sample type
    bam_name_sample_name = balsamic_pon_model.get_final_bam_name(
        bam_dir, sample_name=sample_name
    )

    # Then retrieved final bam names should match the expected format and be identical regardless of request parameter
    expected_final_bam_name = f"{bam_dir}{sample_type}.{sample_name}.dedup.bam"
    assert expected_final_bam_name == bam_name_sample_name
