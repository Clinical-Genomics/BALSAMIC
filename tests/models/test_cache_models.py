"""Test reference cache models."""
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, List

import pytest
from pydantic import ValidationError

from BALSAMIC.commands.init.utils import get_containers
from BALSAMIC.constants.cache import (
    ContainerVersion,
    GenomeVersion,
    REFERENCE_FILES,
    FileType,
    BwaIndexFileType,
)
from BALSAMIC.constants.analysis import BIOINFO_TOOL_ENV
from BALSAMIC.models.cache import (
    CacheConfigModel,
    AnalysisReferencesModel,
    CanFamAnalysisReferencesModel,
    HgAnalysisReferencesModel,
    ReferenceUrlModel,
    ReferencesModel,
    CanFamReferencesModel,
    HgReferencesModel,
)


def test_analysis_references_model(analysis_references_model_data: Dict[str, Path]):
    """Test common analysis references model."""

    # GIVEN an input for the analysis reference model

    # WHEN initialising the model
    model: AnalysisReferencesModel = AnalysisReferencesModel(
        **analysis_references_model_data
    )

    # THEN the model should have been correctly constructed
    assert model.dict() == analysis_references_model_data


def test_analysis_references_model_empty():
    """Test common analysis references model for an empty input."""

    # GIVEN no input for the analysis reference model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        AnalysisReferencesModel()


def test_canfam_analysis_references_model(
    analysis_references_model_data: Dict[str, Path]
):
    """Test canine analysis references model."""

    # GIVEN an input for the canine analysis reference model

    # WHEN initialising the model
    model: CanFamAnalysisReferencesModel = CanFamAnalysisReferencesModel(
        **analysis_references_model_data
    )

    # THEN the model should have been correctly constructed
    assert model.dict() == analysis_references_model_data


def test_canfam_analysis_references_model_empty():
    """Test canine analysis references model for an empty input."""

    # GIVEN no input for the canine analysis reference model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        CanFamAnalysisReferencesModel()


def test_hg_analysis_references_model(
    hg_analysis_references_model_data: Dict[str, Path]
):
    """Test human genome analysis references model."""

    # GIVEN an input for the human genome analysis reference model

    # WHEN initialising the model
    model: HgAnalysisReferencesModel = HgAnalysisReferencesModel(
        **hg_analysis_references_model_data
    )

    # THEN the model should have been correctly constructed
    assert model.dict() == hg_analysis_references_model_data


def test_hg_analysis_references_model_empty():
    """Test human genome analysis references model for an empty input."""

    # GIVEN no input for the human genome analysis reference model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        HgAnalysisReferencesModel()


def test_reference_url_model(reference_url_model_data: Dict[str, Any]):
    """Test references URL model."""

    # GIVEN an input for the reference URL model

    # WHEN initialising the model
    model: ReferenceUrlModel = ReferenceUrlModel(**reference_url_model_data)

    # THEN the model should have been correctly constructed
    assert model.dict() == reference_url_model_data


def test_reference_url_model_empty():
    """Test references URL model for an empty input."""

    # GIVEN no input for the references URL model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        ReferenceUrlModel()


def test_references_model(references_model_data: Dict[str, dict]):
    """Test references model."""

    # GIVEN an input for the reference model

    # WHEN initialising the model
    model: ReferencesModel = ReferencesModel(**references_model_data)

    # THEN the model should have been correctly constructed
    assert model.dict() == references_model_data


def test_references_model_empty():
    """Test references model for an empty input."""

    # GIVEN no input for the references model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        ReferencesModel()


def test_get_reference_genome_files(
    references_model: ReferencesModel, reference_genome_file: Path
):
    """Test reference genome files retrieval."""

    # GIVEN a references model with a specific reference genome file
    references_model.reference_genome.file_path = reference_genome_file.as_posix()

    # GIVEN the expected files to be retrieved
    expected_file_types: set = {FileType.FASTA, FileType.FAI, FileType.DICT}
    expected_file_types.update(BwaIndexFileType)

    # WHEN getting the reference genome files
    reference_genome_files: List[str] = references_model.get_reference_genome_files()

    # THEN the expected reference genome files should be returned
    for file_type in expected_file_types:
        assert file_type in [file.split(".")[-1] for file in reference_genome_files]


def test_get_reference_genome_bwa_index_files(
    references_model: ReferencesModel, reference_genome_file: Path
):
    """Test extraction of reference genome BWA index files."""

    # GIVEN a references model with a specific reference genome file
    references_model.reference_genome.file_path = reference_genome_file.as_posix()

    # GIVEN the expected files to be retrieved
    expected_file_types: set = set(BwaIndexFileType)

    # WHEN getting the reference genome BWA index files
    bwa_index_files: List[str] = references_model.get_reference_genome_bwa_index_files()

    # THEN the expected reference genome BWA index files should be returned
    for file_type in expected_file_types:
        assert file_type in [file.split(".")[-1] for file in bwa_index_files]


def test_get_refgene_files(references_model: ReferencesModel, refgene_txt_file: Path):
    """Test extraction of RefSeq's gene files."""

    # GIVEN a references model with a specific RefSeq's gene TXT file
    references_model.refgene_txt.file_path = refgene_txt_file.as_posix()

    # GIVEN the expected files to be retrieved
    expected_file_types: set = {FileType.TXT, FileType.FLAT, FileType.BED}

    # WHEN getting the RefSeq's gene files
    refegene_files: List[str] = references_model.get_refgene_files()

    # THEN the expected RefSeq's gene files should be returned
    for file_type in expected_file_types:
        assert file_type in [file.split(".")[-1] for file in expected_file_types]


def test_get_refgene_flat_file(
    references_model: ReferencesModel, refgene_txt_file: Path
):
    """Test extraction of RefSeq's gene FLAT file."""

    # GIVEN a references model with a specific RefSeq's gene TXT file
    references_model.refgene_txt.file_path = refgene_txt_file.as_posix()

    # WHEN getting the RefSeq's gene FLAT file
    refegene_flat_file: str = references_model.get_refgene_flat_file()

    # THEN the correctly formatted flat file should be returned
    assert refegene_flat_file.split(".")[-1] == FileType.FLAT


def test_get_refgene_bed_file(
    references_model: ReferencesModel, refgene_txt_file: Path
):
    """Test extraction of RefSeq's gene BED file."""

    # GIVEN a references model with a specific RefSeq's gene TXT file
    references_model.refgene_txt.file_path = refgene_txt_file.as_posix()

    # WHEN getting the RefSeq's gene BED file
    refegene_flat_file: str = references_model.get_refgene_bed_file()

    # THEN the correctly formatted flat file should be returned
    assert refegene_flat_file.split(".")[-1] == FileType.BED


def test_canfam_references_model(references_model_data: Dict[str, dict]):
    """Test canine references model."""

    # GIVEN an input for the canine reference model

    # WHEN initialising the model
    model: CanFamReferencesModel = CanFamReferencesModel(**references_model_data)

    # THEN the model should have been correctly constructed
    assert model.dict() == references_model_data


def test_canfam_references_model_empty():
    """Test canine references model for an empty input."""

    # GIVEN no input for the canine references model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        CanFamReferencesModel()


def test_hg_references_model(hg_references_model_data: Dict[str, dict]):
    """Test human genome references model."""

    # GIVEN an input for the human genome reference model

    # WHEN initialising the model
    model: HgReferencesModel = HgReferencesModel(**hg_references_model_data)

    # THEN the model should have been correctly constructed
    assert model.dict() == hg_references_model_data


def test_hg_references_model_empty():
    """Test human genome references model for an empty input."""

    # GIVEN no input for the human genome references model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        HgReferencesModel()


def test_get_delly_files(
    hg_references_model: HgReferencesModel,
    delly_exclusion_file: Path,
    delly_exclusion_converted_file: Path,
    delly_mappability_file: Path,
    delly_mappability_findex_file: Path,
    delly_mappability_gindex_file: Path,
):
    """Test Delly specific files retrieval."""

    # GIVEN a human genome references model with delly specific files
    hg_references_model.delly_exclusion.file_path = delly_exclusion_file.as_posix()
    hg_references_model.delly_mappability.file_path = delly_mappability_file.as_posix()
    hg_references_model.delly_mappability_findex.file_path = (
        delly_mappability_findex_file.as_posix()
    )
    hg_references_model.delly_mappability_gindex.file_path = (
        delly_mappability_gindex_file.as_posix()
    )

    # WHEN getting the Delly specific reference files
    delly_files: List[str] = hg_references_model.get_delly_files()

    # THEN all the delly reference files should be returned
    assert len(delly_files) == 5
    assert delly_exclusion_file.as_posix() in delly_files
    assert delly_exclusion_converted_file.as_posix() in delly_files
    assert delly_mappability_file.as_posix() in delly_files
    assert delly_mappability_findex_file.as_posix() in delly_files
    assert delly_mappability_gindex_file.as_posix() in delly_files


def test_get_delly_exclusion_converted_file(
    hg_references_model: HgReferencesModel,
    delly_exclusion_file: Path,
    delly_exclusion_converted_file: Path,
):
    """Test get Delly exclusion converted file."""

    # GIVEN a human genome references model with a specific delly exclusion file
    hg_references_model.delly_exclusion.file_path = delly_exclusion_file.as_posix()

    # WHEN getting the Delly exclusion converted file
    converted_file: str = hg_references_model.get_delly_exclusion_converted_file()

    # THEN the returned file should match the expected one
    assert converted_file == delly_exclusion_converted_file.as_posix()


def test_get_gnomad_files(
    hg_references_model: HgReferencesModel,
    gnomad_variant_file: Path,
    gnomad_variant_index_file: Path,
):
    """Test get gnomad reference files."""

    # GIVEN a human genome references model and specific gnomad files
    hg_references_model.gnomad_variant.file_path = gnomad_variant_file.as_posix()
    hg_references_model.gnomad_variant_index.file_path = (
        gnomad_variant_index_file.as_posix()
    )

    # WHEN getting the genomad reference files
    gnomad_files: List[str] = hg_references_model.get_gnomad_files()

    # THEN the genomad files should be returned
    assert len(gnomad_files) == 2
    assert gnomad_variant_file.as_posix() in gnomad_files
    assert gnomad_variant_index_file.as_posix() in gnomad_files


def test_get_1k_genome_files(hg_references_model: HgReferencesModel):
    """Test get 1000 Genome related files."""

    # GIVEN a human genome references model

    # WHEN getting the 1000 Genome files
    genome_files: List[str] = hg_references_model.get_1k_genome_files()

    # THEN


# def test_get_reference_output_files():
#     # GIVEN a reference genome version
#     genome_ver = "hg38"
#     file_type = "fasta"
#
#     # WHEN getting list of valid types
#     fasta_files = get_reference_output_files(REFERENCE_FILES[genome_ver], file_type)
#
#     # THEN it should return list of file
#     assert "Homo_sapiens_assembly38.fasta" in fasta_files
#
#
# def test_referencemeta():
#     """test ReferenceMeta for correctly building model"""
#     # GIVEN a reference model
#     reference_files = {
#         "basedir": "basedir",
#         "reference_genome": {
#             "url": "gs://some_path/b37/human_g1k_v37.fasta.gz",
#             "file_type": "fasta",
#             "gzip": True,
#             "genome_version": "hg19",
#             "output_file": "genome.fa",
#             "output_path": "genome",
#         },
#         "dbsnp": {
#             "url": "gs://some_path/b37/dbsnp_138.b37.vcf.gz",
#             "file_type": "fasta",
#             "gzip": True,
#             "genome_version": "hg19",
#             "output_file": "dbsnp.vcf",
#         },
#     }
#
#     # WHEN build the model
#     build_model = ReferenceMeta.parse_obj(reference_files)
#
#     # THEN model should have correct attributes
#     assert build_model.reference_genome.genome_version == "hg19"
#     assert build_model.dbsnp.genome_version == "hg19"
#     assert build_model.reference_genome.get_output_file == "basedir/genome/genome.fa"
#
#
# def test_referenceurlsmodel_build_model():
#     """test ReferenceUrlsModel for correctly building the model"""
#     # GIVEN a reference model
#     dummy_output_file = "some_random_file"
#     dummy_output_path = "some_path"
#     actual_path = Path(dummy_output_path, dummy_output_file).as_posix()
#
#     dummy_reference = {
#         "url": "gs://domain/file_name",
#         "file_type": "fasta",
#         "gzip": True,
#         "genome_version": "hg19",
#         "output_file": dummy_output_file,
#         "output_path": dummy_output_path,
#     }
#
#     # WHEN building the model
#     built_model = ReferenceUrlsModel.parse_obj(dummy_reference)
#
#     # THEN model should have correct attributes
#     assert built_model.url.scheme == "gs"
#     assert built_model.get_output_file == actual_path
#
#
# def test_referenceurlsmodel_validate_file_type():
#     """test ReferenceUrlsModel for validating file type"""
#     # GIVEN a reference model
#     dummy_output_file = "some_random_file"
#     dummy_output_path = "some_path"
#     actual_path = Path(dummy_output_path, dummy_output_file).as_posix()
#
#     dummy_reference = {
#         "url": "gs://domain/file_name",
#         "file_type": "wrong_type",
#         "gzip": True,
#         "genome_version": "hg19",
#         "output_file": dummy_output_file,
#         "output_path": dummy_output_path,
#     }
#
#     # WHEN building the model
#
#     # THEN model raise error on validation
#     with pytest.raises(ValidationError) as excinfo:
#         built_model = ReferenceUrlsModel.parse_obj(dummy_reference)
#         assert "not a valid reference file format" in excinfo.value
#
#
# def test_referenceurlsmodel_write_md5(tmp_path_factory):
#     """test ReferenceUrlsModel for writing md5 of the output file"""
#     # GIVEN a reference model
#     dummy_output_file = "some_random_file"
#     dummy_output_path = tmp_path_factory.mktemp("some_path")
#     Path(dummy_output_path, dummy_output_file).write_bytes(os.urandom(8196))
#
#     actual_md5_file = Path(dummy_output_path, dummy_output_file + ".md5")
#
#     dummy_reference = {
#         "url": "gs://domain/file_name",
#         "file_type": "fasta",
#         "gzip": True,
#         "genome_version": "hg19",
#         "output_file": dummy_output_file,
#         "output_path": dummy_output_path.as_posix(),
#     }
#
#     # WHEN building the model
#     built_model = ReferenceUrlsModel.parse_obj(dummy_reference)
#
#     # THEN when md5 of the file should exist
#     built_model.write_md5
#     assert actual_md5_file.is_file()
#
#
# def test_referenceurlsmodel_write_md5_no_output_file(tmp_path_factory):
#     """test ReferenceUrlsModel for failing to write md5 if outputfile doesn't exist"""
#     # GIVEN a reference model
#     dummy_output_file = "some_random_file"
#     dummy_output_path = tmp_path_factory.mktemp("some_path")
#
#     actual_md5_file = Path(dummy_output_path, dummy_output_file + ".md5")
#
#     dummy_reference = {
#         "url": "gs://domain/file_name",
#         "file_type": "fasta",
#         "gzip": True,
#         "genome_version": "hg19",
#         "output_file": dummy_output_file,
#         "output_path": dummy_output_path.as_posix(),
#     }
#
#     # WHEN building the model
#     built_model = ReferenceUrlsModel.parse_obj(dummy_reference)
#
#     # THEN when md5 of the file should exist
#     with pytest.raises(FileNotFoundError) as excinfo:
#         built_model.write_md5
#         assert "file does not exist" in excinfo.value
#
#
# def test_referenceurlsmodel_validate_genome_version():
#     """test ReferenceUrlsModel for validating genome version"""
#     # GIVEN a reference model
#     dummy_output_file = "some_random_file"
#     dummy_output_path = "some_path"
#     actual_path = Path(dummy_output_path, dummy_output_file).as_posix()
#
#     dummy_reference = {
#         "url": "gs://domain/file_name",
#         "file_type": "fasta",
#         "gzip": True,
#         "genome_version": "wrong_genome",
#         "output_file": dummy_output_file,
#         "output_path": dummy_output_path,
#     }
#
#     with pytest.raises(ValidationError) as excinfo:
#         # WHEN building the model
#         built_model = ReferenceUrlsModel.parse_obj(dummy_reference)
#
#         # THEN model raise error on validation
#         assert "not a valid genome version" in excinfo.value
