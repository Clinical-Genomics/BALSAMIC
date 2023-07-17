import copy
import os
import pytest

from pathlib import Path
from pydantic import ValidationError

from BALSAMIC.utils.models import (
    VCFAttributes,
    VarCallerFilter,
    QCModel,
    VarcallerAttribute,
    AnalysisModel,
    SampleInstanceModel,
    ReferenceUrlsModel,
    ReferenceMeta,
    UMIParamsCommon,
    UMIParamsUMIextract,
    UMIParamsConsensuscall,
    UMIParamsTNscope,
    ParamsVardict,
    ParamsVEP,
    MetricModel,
    MetricConditionModel,
    MetricValidationModel,
    AnalysisPonModel,
)
from tests.conftest import analysis_dir


def test_referencemeta():
    """test ReferenceMeta for correctly building model"""
    # GIVEN a reference model
    reference_files = {
        "basedir": "basedir",
        "reference_genome": {
            "url": "gs://some_path/b37/human_g1k_v37.fasta.gz",
            "file_type": "fasta",
            "gzip": True,
            "genome_version": "hg19",
            "file_name": "genome.fa",
            "dir_name": "genome",
        },
        "dbsnp": {
            "url": "gs://some_path/b37/dbsnp_138.b37.vcf.gz",
            "file_type": "fasta",
            "gzip": True,
            "genome_version": "hg19",
            "file_name": "dbsnp.vcf",
        },
    }

    # WHEN build the model
    build_model = ReferenceMeta.parse_obj(reference_files)

    # THEN model should have correct attributes
    assert build_model.reference_genome.genome_version == "hg19"
    assert build_model.dbsnp.genome_version == "hg19"
    assert build_model.reference_genome.get_output_file == "basedir/genome/genome.fa"


def test_referenceurlsmodel_build_model():
    """test ReferenceUrlsModel for correctly building the model"""
    # GIVEN a reference model
    dummy_output_file = "some_random_file"
    dummy_output_path = "some_path"
    actual_path = Path(dummy_output_path, dummy_output_file).as_posix()

    dummy_reference = {
        "url": "gs://domain/file_name",
        "file_type": "fasta",
        "gzip": True,
        "genome_version": "hg19",
        "file_name": dummy_output_file,
        "dir_name": dummy_output_path,
    }

    # WHEN building the model
    built_model = ReferenceUrlsModel.parse_obj(dummy_reference)

    # THEN model should have correct attributes
    assert built_model.url.scheme == "gs"
    assert built_model.get_output_file == actual_path


def test_referenceurlsmodel_validate_file_type():
    """test ReferenceUrlsModel for validating file type"""
    # GIVEN a reference model
    dummy_output_file = "some_random_file"
    dummy_output_path = "some_path"
    actual_path = Path(dummy_output_path, dummy_output_file).as_posix()

    dummy_reference = {
        "url": "gs://domain/file_name",
        "file_type": "wrong_type",
        "gzip": True,
        "genome_version": "hg19",
        "file_name": dummy_output_file,
        "dir_name": dummy_output_path,
    }

    # WHEN building the model

    # THEN model raise error on validation
    with pytest.raises(ValidationError) as excinfo:
        built_model = ReferenceUrlsModel.parse_obj(dummy_reference)
        assert "not a valid reference file format" in excinfo.value


def test_referenceurlsmodel_write_md5(tmp_path_factory):
    """test ReferenceUrlsModel for writing md5 of the output file"""
    # GIVEN a reference model
    dummy_output_file = "some_random_file"
    dummy_output_path = tmp_path_factory.mktemp("some_path")
    Path(dummy_output_path, dummy_output_file).write_bytes(os.urandom(8196))

    actual_md5_file = Path(dummy_output_path, dummy_output_file + ".md5")

    dummy_reference = {
        "url": "gs://domain/file_name",
        "file_type": "fasta",
        "gzip": True,
        "genome_version": "hg19",
        "file_name": dummy_output_file,
        "dir_name": dummy_output_path.as_posix(),
    }

    # WHEN building the model
    built_model = ReferenceUrlsModel.parse_obj(dummy_reference)

    # THEN when md5 of the file should exist
    built_model.write_md5
    assert actual_md5_file.is_file()


def test_referenceurlsmodel_write_md5_no_output_file(tmp_path_factory):
    """test ReferenceUrlsModel for failing to write md5 if outputfile doesn't exist"""
    # GIVEN a reference model
    dummy_output_file = "some_random_file"
    dummy_output_path = tmp_path_factory.mktemp("some_path")

    actual_md5_file = Path(dummy_output_path, dummy_output_file + ".md5")

    dummy_reference = {
        "url": "gs://domain/file_name",
        "file_type": "fasta",
        "gzip": True,
        "genome_version": "hg19",
        "file_name": dummy_output_file,
        "dir_name": dummy_output_path.as_posix(),
    }

    # WHEN building the model
    built_model = ReferenceUrlsModel.parse_obj(dummy_reference)

    # THEN when md5 of the file should exist
    with pytest.raises(FileNotFoundError) as excinfo:
        built_model.write_md5
        assert "file does not exist" in excinfo.value


def test_referenceurlsmodel_validate_genome_version():
    """test ReferenceUrlsModel for validating genome version"""
    # GIVEN a reference model
    dummy_output_file = "some_random_file"
    dummy_output_path = "some_path"
    actual_path = Path(dummy_output_path, dummy_output_file).as_posix()

    dummy_reference = {
        "url": "gs://domain/file_name",
        "file_type": "fasta",
        "gzip": True,
        "genome_version": "wrong_genome",
        "file_name": dummy_output_file,
        "dir_name": dummy_output_path,
    }

    with pytest.raises(ValidationError) as excinfo:
        # WHEN building the model
        built_model = ReferenceUrlsModel.parse_obj(dummy_reference)

        # THEN model raise error on validation
        assert "not a valid genome version" in excinfo.value


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
    with pytest.raises(ValueError) as excinfo:
        VarcallerAttribute.parse_obj(invalid_args)
        assert "not a valid argument" in excinfo.value


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
        assert "not supported" in excinfo.value


def test_sample_instance_model():
    """Test sample instance model initialisation."""

    # GIVEN a sample dictionary
    tumor_sample: dict = {"ACC1": {"type": "tumor"}}

    # WHEN parsing the sample dictionary
    sample: SampleInstanceModel = SampleInstanceModel.parse_obj(tumor_sample["ACC1"])

    # THEN the sample model should be correctly initialised
    assert sample.dict() == tumor_sample["ACC1"]


def test_sample_instance_model_error():
    """Test sample instance model error raise."""

    # GIVEN a sample dictionary with an invalid sample type
    sample_type: str = "affected"
    samples: dict = {"ACC1": {"type": sample_type}}

    # WHEN parsing the sample dictionary

    # THEN a ValueError should be triggered
    with pytest.raises(ValueError) as exc:
        SampleInstanceModel.parse_obj(samples["ACC1"])
        assert (
            f"The provided sample type ({sample_type}) is not supported in BALSAMIC"
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


def test_metric_condition_model():
    """test MetricConditionModel attributes parsing"""

    # GIVEN input attributes
    metric_condition = {"norm": "gt", "threshold": 1}

    # WHEN building the metric condition model
    metrics_model = MetricConditionModel(**metric_condition)

    # THEN assert retrieved values from the created model
    assert metrics_model.dict().items() == metric_condition.items()


def test_metric_model_pass_validation():
    """test MetricModel attributes parsing"""

    # GIVEN input attributes
    metrics = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1.sorted.mrkdup.hsmetric",
        "name": "MEDIAN_TARGET_COVERAGE",
        "step": "multiqc_picard_HsMetrics",
        "value": 2393.0,
        "condition": {"norm": "gt", "threshold": 1000.0},
    }

    # WHEN building the metric model
    metric_model = MetricModel(**metrics)

    # THEN assert retrieved values from the created model
    assert metric_model.dict().items() == metrics.items()


def test_metric_model_duplication_refactoring():
    """test MetricModel duplications param refactoring"""

    # GIVEN input attributes
    metrics = {
        "header": None,
        "id": "ACC1",
        "input": "ACC1_R_1_fastqc.zip",
        "name": "FastQC_mqc-generalstats-fastqc-percent_duplicates",
        "step": "multiqc_general_stats",
        "value": 21.517800000611373,
        "condition": None,
    }

    # WHEN building the metric model
    metric_model = MetricModel(**metrics)

    # THEN assert retrieved values from the created model
    assert metric_model.name == "PERCENT_DUPLICATION_R1"


def test_metric_model_fail_validation():
    """test MetricModel behaviour for an incorrect input"""

    # GIVEN a non accepted input
    invalid_input = {"header": None, "id": "ACC1"}

    # THEN the model raises an error due to an incomplete input
    with pytest.raises(ValueError) as input_exc:
        MetricModel(**invalid_input)
    assert f"field required" in str(input_exc.value)


def test_metric_validation_model_pass(qc_extracted_metrics):
    """test MetricValidationModel attribute parsing and positive validation"""

    # WHEN building the MetricValidationModel model
    model = MetricValidationModel(metrics=qc_extracted_metrics)

    # THEN assert retrieved values from the created model
    assert model.dict()["metrics"] == qc_extracted_metrics


def test_metric_validation_model_fail(qc_extracted_metrics):
    """test MetricValidationModel for an overly restrictive metric condition"""

    # GIVEN input attributes with a value that does not meet the filtering condition
    metrics = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["value"] = 2.0  # GC_DROPOUT set to 2.0 (failing condition)

    # THEN check that the model filters the metric according to its norm
    with pytest.raises(ValueError) as val_exc:
        MetricValidationModel(metrics=metrics)
    assert (
        f"QC metric {metrics[4]['name']}: {metrics[4]['value']} validation has failed. "
        f"(Condition: {metrics[4]['condition']['norm']} {metrics[4]['condition']['threshold']}, ID: {metrics[4]['id']})"
        in str(val_exc.value)
    )


def test_multiple_metric_validation_model_fail(qc_extracted_metrics):
    """test MetricValidationModel for multiple metrics with failing conditions"""

    # GIVEN input attributes that does not meet the specified conditions
    metrics = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["value"] = 2.0  # GC_DROPOUT set to 2.0 (failing condition)
    metrics[8]["value"] = 0.5  # PCT_TARGET_BASES_500X set to 50% (failing condition)

    # THEN check that the model filters the metrics according to its norm
    with pytest.raises(ValueError) as val_exc:
        MetricValidationModel(metrics=metrics)
    assert "2 validation errors for MetricValidationModel" in str(val_exc.value)
    assert metrics[4]["name"] in str(val_exc.value)
    assert metrics[8]["name"] in str(val_exc.value)


def test_metric_validation_model_norm_fail(qc_extracted_metrics):
    """test MetricValidationModel ValueError raising for an operator that it is not accepted"""

    # GIVEN a metric with an incorrect norm attribute
    metrics = copy.deepcopy(qc_extracted_metrics)
    metrics[4]["condition"]["norm"] = "lower"

    # THEN model raises an error due to a non accepted norm
    try:
        MetricValidationModel(metrics=metrics)
    except KeyError as key_exc:
        assert metrics[4]["condition"]["norm"] in str(key_exc)


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
    with pytest.raises(ValueError) as excinfo:
        AnalysisPonModel.parse_obj(invalid_args)
        assert (
            f"The provided version {invalid_args['pon_version']} does not follow the defined syntax (v<int>)"
            in excinfo.value
        )
