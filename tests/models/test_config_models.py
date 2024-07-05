import copy
from datetime import datetime
from pathlib import Path
from typing import Dict, List

import pytest
from _pytest.logging import LogCaptureFixture
from pydantic import ValidationError

from BALSAMIC.constants.analysis import FastqName, SampleType, SequencingType
from BALSAMIC.models.config import AnalysisModel, ConfigModel, SampleInstanceModel


def test_analysis_model(
    test_data_dir: Path, timestamp_now: datetime, caplog: LogCaptureFixture
):
    """Test analysis model instantiation."""

    # GIVEN valid input arguments
    valid_args = {
        "case_id": "case_id",
        "gender": "female",
        "analysis_type": "paired",
        "sequencing_type": "targeted",
        "analysis_dir": test_data_dir.as_posix(),
        "fastq_path": test_data_dir.as_posix(),
        "log": test_data_dir.as_posix(),
        "result": test_data_dir.as_posix(),
        "script": test_data_dir.as_posix(),
        "benchmark": test_data_dir.as_posix(),
        "dag": test_data_dir.as_posix(),
        "config_creation_date": str(timestamp_now),
        "analysis_workflow": "balsamic-umi",
    }

    # THEN we can successfully create a config dict
    analysis_model: AnalysisModel = AnalysisModel.model_validate(valid_args)
    assert analysis_model
    assert analysis_model.script == test_data_dir.as_posix()

    # GIVEN invalid input arguments for the log directory
    valid_args["log"] = "/not/a/directory"

    # THEN should trigger ValueError
    with pytest.raises(ValueError):
        AnalysisModel.model_validate(valid_args)

    assert "The supplied directory path /not/a/directory does not exist" in caplog.text


def test_sample_instance_model(config_dict_w_fastqs: Dict):
    """Test sample instance model initialisation."""

    # GIVEN a sample list
    sample_list = config_dict_w_fastqs["samples"]

    # WHEN parsing the sample dictionary
    for idx, sample in enumerate(sample_list):
        sample: SampleInstanceModel = SampleInstanceModel.model_validate(sample)

        sample_dict_copy = sample_list[idx].copy()
        for fastq_pattern, values in sample_dict_copy["fastq_info"].items():
            values["fwd"] = Path(values["fwd"]).resolve().as_posix()
            values["rev"] = Path(values["rev"]).resolve().as_posix()

        # THEN the sample model should be correctly initialised
        assert sample.model_dump(exclude_none=True) == sample_dict_copy


def test_sample_instance_model_sample_type_error(tumor_normal_fastq_info_correct: Dict):
    """Test sample instance model error raise."""

    # GIVEN a sample dictionary with an invalid sample type
    samples: List[Dict] = copy.deepcopy(tumor_normal_fastq_info_correct)
    illegal_sample_type: str = "affected"
    tumor_dict = samples[0]
    tumor_dict["type"] = illegal_sample_type

    # WHEN parsing the sample dictionary
    with pytest.raises(ValueError) as exc:
        SampleInstanceModel.model_validate(tumor_dict)

    # THEN a ValueError should be triggered
    assert "Input should be 'normal' or 'tumor'" in str(exc.value)


def test_analysis_model_for_pon(test_data_dir: Path, timestamp_now: datetime):
    """Tests PON model parsing."""

    # GIVEN valid input arguments
    valid_args = {
        "case_id": "case_id",
        "analysis_type": "pon",
        "sequencing_type": "targeted",
        "analysis_dir": test_data_dir.as_posix(),
        "fastq_path": test_data_dir.as_posix(),
        "log": test_data_dir.as_posix(),
        "result": test_data_dir.as_posix(),
        "script": test_data_dir.as_posix(),
        "benchmark": test_data_dir.as_posix(),
        "dag": test_data_dir.as_posix(),
        "analysis_workflow": "balsamic",
        "config_creation_date": str(timestamp_now),
        "pon_version": "v1",
    }

    # THEN we can successfully create a config dict
    assert AnalysisModel.model_validate(valid_args)

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
        AnalysisModel.model_validate(invalid_args)

    assert (
        f"The provided PON version ({invalid_args['pon_version']}) does not follow the defined syntax (v<int>)"
        in str(excinfo.value)
    )


def test_detect_duplicate_fastq_pattern(
    config_w_fastq_dir_for_duplicate_fastq_patterns_model: Dict,
):
    """Test balsamic models ability to detect duplicate assigned fastq patterns."""
    config_dict = config_w_fastq_dir_for_duplicate_fastq_patterns_model
    # Initialize balsamic model
    with pytest.raises(ValueError) as exc:
        ConfigModel.model_validate(config_dict)

    assert (
        "Duplicate FastqPattern(s) found: ACC1_S1_L001_R across multiple samples"
        in str(exc.value)
    )


def test_detection_unassigned_fastq_file(config_tumor_normal_extrafile: Dict):
    """Test instantiating balsamic model with fastq dir containing unassigned fastq-files."""
    # Initialize balsamic model
    with pytest.raises(ValueError) as exc:
        ConfigModel.model_validate(config_tumor_normal_extrafile)

    assert "Fastqs in fastq-dir not assigned to sample config:" in str(exc.value)


def test_get_all_sample_names(balsamic_model: ConfigModel):
    """Validate retrieval of all sample names in analysis from ConfigModel."""
    sample_names = balsamic_model.get_all_sample_names()
    assert ["ACC1", "ACC2"] == sample_names


def test_get_fastq_patterns_by_sample(
    balsamic_model: ConfigModel, tumor_sample_name: str, normal_sample_name: str
):
    """Validate retrieval of fastq-pattern by sample from ConfigModel."""

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


def test_get_all_fastqs_for_sample(balsamic_model: ConfigModel, tumor_sample_name: str):
    """Validate retrieval of fastq-files by sample and fastq-type from ConfigModel."""

    def compare_fastq_file_lists(expected: List[str], found: List[str]):
        found_file_names = []
        for found_file in found:
            found_file_names.append(Path(found_file).name)
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
        sample_name=tumor_sample_name, fastq_types=[FastqName.FWD]
    )
    rev_fastq_files = balsamic_model.get_all_fastqs_for_sample(
        sample_name=tumor_sample_name, fastq_types=[FastqName.REV]
    )
    fastq_files = balsamic_model.get_all_fastqs_for_sample(
        sample_name=tumor_sample_name
    )

    compare_fastq_file_lists(fwd_fastq_files_expected, fwd_fastq_files)
    compare_fastq_file_lists(rev_fastq_files_expected, rev_fastq_files)
    compare_fastq_file_lists(fastq_files_expected, fastq_files)
    assert normal_fastq not in fastq_files


def test_get_all_fastq_names(balsamic_model: ConfigModel):
    """Validate retrieval of all fastq-files from ConfigModel and optional removal of suffix."""

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


def test_fastq_by_fastq_pattern(balsamic_model: ConfigModel):
    """Validate retrieval of fastq-file by fastq-pattern and fastq-type from ConfigModel."""

    fastq_pattern = "2_171015_HJ7TLDSX5_ACC2_XXXXXX"
    expected_fwd = "2_171015_HJ7TLDSX5_ACC2_XXXXXX_1.fastq.gz"
    expected_rev = "2_171015_HJ7TLDSX5_ACC2_XXXXXX_2.fastq.gz"

    fwd_fastq = balsamic_model.get_fastq_by_fastq_pattern(fastq_pattern, FastqName.FWD)
    rev_fastq = balsamic_model.get_fastq_by_fastq_pattern(fastq_pattern, FastqName.REV)

    assert Path(fwd_fastq).name == expected_fwd
    assert Path(rev_fastq).name == expected_rev


def test_sample_name_by_type(balsamic_model: ConfigModel):
    """Validate retrieval of sample name by sample type from ConfigModel."""

    tumor_name_expected = "ACC1"
    normal_name_expected = "ACC2"

    # Given sample type
    tumor_name_retrieved = balsamic_model.get_sample_name_by_type(SampleType.TUMOR)
    normal_name_retrieved = balsamic_model.get_sample_name_by_type(SampleType.NORMAL)

    # Then the retrieved sample name should match the expected
    assert tumor_name_retrieved == tumor_name_expected
    assert normal_name_retrieved == normal_name_expected


def test_sample_type_by_name(balsamic_model: ConfigModel):
    """Validate retrieval of sample type by sample name from ConfigModel."""

    tumor_name = "ACC1"
    normal_name = "ACC2"

    # Given sample name
    tumor_type_retrieved = balsamic_model.get_sample_type_by_name(tumor_name)
    normal_type_retrieved = balsamic_model.get_sample_type_by_name(normal_name)

    # Then the retrieved sample type should match the expected
    assert tumor_type_retrieved == SampleType.TUMOR
    assert normal_type_retrieved == SampleType.NORMAL


def test_get_bam_name_per_lane(balsamic_model: ConfigModel):
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
    bam_dir = Path(result_dir, "bam", "").as_posix()

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


def test_get_final_bam_name(balsamic_model: ConfigModel):
    """Validate retrieval of final bam name by either sample type or sample name."""

    # Given bam_dir path and sample name or sample type
    sample_name = "ACC1"
    sample_type = SampleType.TUMOR
    result_dir = balsamic_model.analysis.result
    bam_dir = Path(result_dir, "bam", "").as_posix()

    # When retrieving final bam file name by sample name or sample type
    bam_name_sample_name = balsamic_model.get_final_bam_name(
        bam_dir, sample_name=sample_name
    )
    bam_name_sample_type = balsamic_model.get_final_bam_name(
        bam_dir, sample_type=sample_type
    )

    # Then retrieved final bam names should match the expected format and be identical regardless of request parameter
    expected_final_bam_name = f"{bam_dir}{sample_type}.{sample_name}.dedup.fixmate.bam"
    assert expected_final_bam_name == bam_name_sample_name
    assert bam_name_sample_name == bam_name_sample_type

    # WHEN changing sequencing_type to WGS
    balsamic_model.analysis.sequencing_type = SequencingType.WGS
    # Then retrieved final bam names should have realignment suffix
    expected_final_bam_name = f"{bam_dir}{sample_type}.{sample_name}.dedup.realign.bam"
    bam_name_sample_type = balsamic_model.get_final_bam_name(
        bam_dir, sample_type=sample_type
    )
    assert expected_final_bam_name == bam_name_sample_type

    # WHEN submitting custom bam suffix
    bam_name_sample_type = balsamic_model.get_final_bam_name(
        bam_dir, sample_type=sample_type, specified_suffix="dedup.fixmate.qualcapped"
    )
    # Then the bam name should end with the specified suffix
    expected_final_bam_name = f"{bam_dir}{sample_type}.{sample_name}.dedup.fixmate.qualcapped.bam"
    assert expected_final_bam_name == bam_name_sample_type


def test_no_info_error_get_final_bam_name(balsamic_model: ConfigModel):
    """Validate raise ValueError by not sample type or sample name."""

    # Given bam_dir path
    result_dir = balsamic_model.analysis.result
    bam_dir = Path(result_dir, "bam").as_posix()

    # When retrieving final bam file name without supplying sample name or sample type
    # ValueError should be raised
    with pytest.raises(ValueError) as excinfo:
        balsamic_model.get_final_bam_name(bam_dir)
    assert (
        "Either sample_name or sample_type must be provided to get the final bam name."
        in str(excinfo.value)
    )


def test_get_final_bam_name_pon(balsamic_pon_model: ConfigModel):
    """Validate retrieval of final bam name for PON by either sample type or sample name."""

    # Given bam_dir path and sample name or sample type
    sample_name = "ACCN6"
    sample_type = SampleType.NORMAL
    result_dir = balsamic_pon_model.analysis.result
    bam_dir = Path(result_dir, "bam").as_posix()

    # When retrieving final bam file name by sample name or sample type
    bam_name_sample_name = balsamic_pon_model.get_final_bam_name(
        bam_dir, sample_name=sample_name
    )

    # Then retrieved final bam names should match the expected format and be identical regardless of request parameter
    expected_final_bam_name = f"{bam_dir}{sample_type}.{sample_name}.dedup.bam"
    assert expected_final_bam_name == bam_name_sample_name
