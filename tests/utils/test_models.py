import os
import pytest

from pathlib import Path
from pydantic import ValidationError

from BALSAMIC.utils.models import (
    VCFAttributes, VarCallerFilter, QCModel, VarcallerAttribute, AnalysisModel,
    SampleInstanceModel, ReferenceUrlsModel, ReferenceMeta, UMIworkflowConfig,
    UMIParamsCommon, UMIParamsUMIextract, UMIParamsConsensuscall,
    UMIParamsTNscope, UMIParamsVardict, UMIParamsVEP)


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
            "output_file": "genome.fa",
            "output_path": "genome",
        },
        "dbsnp": {
            "url": "gs://some_path/b37/dbsnp_138.b37.vcf.gz",
            "file_type": "fasta",
            "gzip": True,
            "genome_version": "hg19",
            "output_file": "dbsnp.vcf",
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
        "output_file": dummy_output_file,
        "output_path": dummy_output_path,
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
        "output_file": dummy_output_file,
        "output_path": dummy_output_path,
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
        "output_file": dummy_output_file,
        "output_path": dummy_output_path.as_posix(),
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
        "output_file": dummy_output_file,
        "output_path": dummy_output_path.as_posix(),
    }

    # WHEN building the model
    built_model = ReferenceUrlsModel.parse_obj(dummy_reference)

    # THEN when md5 of the file should exist
    with pytest.raises(FileNotFoundError) as excinfo:
        built_model.write_md5
        assert "file does not exist" in excinfo.value


def test_referenceurlsmodel_validate_genome_version():
    """test ReferenceUrlsModel for validating genome version """
    # GIVEN a reference model
    dummy_output_file = "some_random_file"
    dummy_output_path = "some_path"
    actual_path = Path(dummy_output_path, dummy_output_file).as_posix()

    dummy_reference = {
        "url": "gs://domain/file_name",
        "file_type": "fasta",
        "gzip": True,
        "genome_version": "wrong_genome",
        "output_file": dummy_output_file,
        "output_path": dummy_output_path,
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
        "pop_freq": {
            "tag_value": 0.005,
            "filter_name": "dummy_pop_freq",
            "field": "INFO"
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
    valid_args = {"umi_trim": True, "min_seq_length": 25, "umi_trim_length": 5}
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


def test_analysis_model():
    # GIVEN valid input arguments
    valid_args = {
        "case_id": "case_id",
        "analysis_type": "paired",
        "sequencing_type": "targeted",
        "analysis_dir": "tests/test_data",
    }
    # THEN we can successully create a config dict
    assert AnalysisModel.parse_obj(valid_args)

    # GIVEN invalid input arguments
    invalid_args = {
        "case_id": "case_id",
        "analysis_type": "odd",
        "sequencing_type": "wrong",
        "analysis_dir": "tests/test_data",
    }
    # THEN should trigger ValueError
    with pytest.raises(ValueError) as excinfo:
        AnalysisModel.parse_obj(invalid_args)
        assert "not supported" in excinfo.value


def test_sample_instance_model():
    # GIVEN valid input arguments
    valid_args = {"file_prefix": "S2_R", "type": "normal", "sample_name": "S2"}
    # THEN we can successully create a config dict
    assert SampleInstanceModel.parse_obj(valid_args)

    # GIVEN invalid input arguments
    invalid_args = {
        "file_prefix": "S2_R",
        "type": "fungal",
    }
    # THEN should trigger ValueError
    with pytest.raises(ValueError) as excinfo:
        SampleInstanceModel.parse_obj(invalid_args)
        assert "not supported" in excinfo.value


def test_umiparams_common():
    """ test UMIParamsCommon model for correct validation """

    # GIVEN a UMI workflow common params
    test_commonparams = {
        "align_header": "test_header_name",
        "align_intbases": 100,
        "filter_tumor_af": 0.01
    }
    # WHEN building the model
    test_commonparams_built = UMIParamsCommon(**test_commonparams)
    # THEN assert values
    assert test_commonparams_built.align_header == "test_header_name"
    assert test_commonparams_built.filter_tumor_af == 0.01
    assert test_commonparams_built.align_intbases == 100


def test_umiparams_umiextract():
    """ test UMIParamsUMIextract model for correct validation """
    # GIVEN umiextract params
    test_umiextractparams = {"read_structure": "['mode', 'r1,r2']"}

    # WHEN building the model
    test_umiextractparams_built = UMIParamsUMIextract(**test_umiextractparams)

    # THEN assert values
    assert test_umiextractparams_built.read_structure == "['mode', 'r1,r2']"


def test_umiparams_consensuscall():
    """ test UMIParamsConsensuscall model for correct validation """

    #GIVEN consensuscall params
    test_consensuscall = {
        "align_format": "BAM",
        "filter_minreads": "6,3,3",
        "tag": "XZ"
    }

    #WHEN building the model
    test_consensuscall_built = UMIParamsConsensuscall(**test_consensuscall)

    #THEN assert values
    assert test_consensuscall_built.align_format == "BAM"
    assert test_consensuscall_built.filter_minreads == "6,3,3"
    assert test_consensuscall_built.tag == "XZ"


def test_umiparams_tnscope():
    """ test UMIParamsTNscope model for correct validation """

    #GIVEN tnscope params
    test_tnscope_params = {
        "algo": "algoname",
        "min_tumorLOD": 6,
        "error_rate": 5,
        "prunefactor": 3,
        "disable_detect": "abc"
    }

    #WHEN building the model
    test_tnscope_params_built = UMIParamsTNscope(**test_tnscope_params)

    #THEN assert values
    assert test_tnscope_params_built.algo == "algoname"
    assert test_tnscope_params_built.min_tumorLOD == 6
    assert test_tnscope_params_built.error_rate == 5
    assert test_tnscope_params_built.prunefactor == 3
    assert test_tnscope_params_built.disable_detect == "abc"


def test_umiparams_vardict():
    """ test UMIParamsVardict model for correct validation"""

    #GIVEN vardict params
    test_umivardict = {"vardict_filters": "-a 1 -b 2 -c 5"}

    #WHEN building the model
    test_umivardict_built = UMIParamsVardict(**test_umivardict)

    #THEN assert values
    assert test_umivardict_built.vardict_filters == "-a 1 -b 2 -c 5"


def test_umiparams_vep():
    """ test UMIParamsVEP model for correct validation"""

    #GIVEN vardict params
    test_umivep = {"vep_filters": "all defaults params"}

    #WHEN building the model
    test_umivep_built = UMIParamsVEP(**test_umivep)

    #THEN assert values
    assert test_umivep_built.vep_filters == "all defaults params"
