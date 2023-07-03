"""Test reference cache models."""
from pathlib import Path
from typing import Dict, Any, List

import pytest
from _pytest.logging import LogCaptureFixture

from BALSAMIC.constants.cache import (
    GRCHVersion,
    DockerContainers,
    GenomeVersion,
)
from BALSAMIC.constants.constants import FileType, BwaIndexFileType
from pydantic import ValidationError

from BALSAMIC.models.cache import (
    AnalysisReferencesModel,
    CanFamAnalysisReferencesModel,
    HgAnalysisReferencesModel,
    ReferenceUrlModel,
    ReferencesModel,
    CanFamReferencesModel,
    HgReferencesModel,
    CacheAnalysisModel,
    CacheConfigModel,
)
from BALSAMIC.utils.exc import BalsamicError


def test_analysis_references_model(analysis_references_model_data: Dict[str, Path]):
    """Test common analysis references model."""

    # GIVEN an input for the analysis reference model

    # WHEN initialising the model
    model: AnalysisReferencesModel = AnalysisReferencesModel(
        **analysis_references_model_data
    )

    # THEN the model should have been correctly built
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

    # THEN the model should have been correctly built
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

    # THEN the model should have been correctly built
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

    # THEN the model should have been correctly built
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

    # THEN the model should have been correctly built
    assert model.dict() == references_model_data


def test_references_model_empty():
    """Test references model for an empty input."""

    # GIVEN no input for the references model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        ReferencesModel()


def test_get_reference_genome_files(references_model: ReferencesModel):
    """Test reference genome files retrieval."""

    # GIVEN a references model

    # GIVEN the expected files to be retrieved
    expected_file_types: set = {FileType.FASTA, FileType.FAI, FileType.DICT}
    expected_file_types.update(BwaIndexFileType)

    # WHEN getting the reference genome files
    reference_genome_files: List[str] = references_model.get_reference_genome_files()

    # THEN the expected reference genome files should be returned
    assert len(reference_genome_files) == len(expected_file_types)
    for file_type in expected_file_types:
        assert file_type in [file.split(".")[-1] for file in reference_genome_files]


def test_get_reference_genome_bwa_index_files(references_model: ReferencesModel):
    """Test extraction of reference genome BWA index files."""

    # GIVEN a references model

    # GIVEN the expected files to be retrieved
    expected_file_types: set = set(BwaIndexFileType)

    # WHEN getting the reference genome BWA index files
    bwa_index_files: List[str] = references_model.get_reference_genome_bwa_index_files()

    # THEN the expected reference genome BWA index files should be returned
    assert len(bwa_index_files) == len(expected_file_types)
    for file_type in expected_file_types:
        assert file_type in [file.split(".")[-1] for file in bwa_index_files]


def test_get_refgene_files(
    references_model: ReferencesModel, refgene_bed_file: Path, refgene_flat_file: Path
):
    """Test extraction of RefSeq's gene files."""

    # GIVEN a references model and some  mocked RefSeq's gene file

    # WHEN getting the RefSeq's gene files
    refgene_files: List[str] = references_model.get_refgene_files()

    # THEN the expected RefSeq's gene files should be returned
    assert len(refgene_files) == 3
    assert references_model.refgene_txt.file_path in refgene_files
    assert refgene_bed_file.as_posix() in refgene_files
    assert refgene_flat_file.as_posix() in refgene_files


def test_get_refgene_flat_file(
    references_model: ReferencesModel, refgene_flat_file: Path
):
    """Test extraction of RefSeq's gene FLAT file."""

    # GIVEN a references model and a mocked RefSeq's gene FLAT file

    # WHEN getting the RefSeq's gene FLAT file
    refgene_output_file: str = references_model.get_refgene_flat_file()

    # THEN the correctly formatted flat file should be returned
    assert refgene_output_file == refgene_flat_file.as_posix()


def test_get_refgene_bed_file(
    references_model: ReferencesModel, refgene_bed_file: Path
):
    """Test extraction of RefSeq's gene BED file."""

    # GIVEN a references model and a mocked RefSeq's gene BED file

    # WHEN getting the RefSeq's gene BED file
    refgene_output_file: str = references_model.get_refgene_bed_file()

    # THEN the correctly formatted flat file should be returned
    assert refgene_output_file == refgene_bed_file.as_posix()


def test_canfam_references_model(references_model_data: Dict[str, dict]):
    """Test canine references model."""

    # GIVEN an input for the canine reference model

    # WHEN initialising the model
    model: CanFamReferencesModel = CanFamReferencesModel(**references_model_data)

    # THEN the model should have been correctly built
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

    # THEN the model should have been correctly built
    assert model.dict() == hg_references_model_data


def test_hg_references_model_empty():
    """Test human genome references model for an empty input."""

    # GIVEN no input for the human genome references model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        HgReferencesModel()


def test_get_cadd_snv_files(
    hg_references_model: HgReferencesModel, cadd_snv_indexed_file: Path
):
    """Test get CADD SNV reference output files."""

    # GIVEN a human genome references model and a mocked CADD SNV indexed file

    # WHEN getting the CADD specific reference files
    cadd_snv_files: List[str] = hg_references_model.get_cadd_snv_files()

    # THEN all the CADD SNV reference files should be returned
    assert len(cadd_snv_files) == 2
    assert hg_references_model.cadd_snv.file_path in cadd_snv_files
    assert cadd_snv_indexed_file.as_posix() in cadd_snv_files


def test_get_delly_files(
    hg_references_model: HgReferencesModel, delly_exclusion_converted_file: Path
):
    """Test Delly specific files retrieval."""

    # GIVEN a human genome references model and a mocked Delly exclusion converted file

    # WHEN getting the Delly specific reference files
    delly_files: List[str] = hg_references_model.get_delly_files()

    # THEN all the delly reference files should be returned
    assert len(delly_files) == 5
    assert hg_references_model.delly_exclusion.file_path in delly_files
    assert delly_exclusion_converted_file.as_posix() in delly_files
    assert hg_references_model.delly_mappability.file_path in delly_files
    assert hg_references_model.delly_mappability_findex.file_path in delly_files
    assert hg_references_model.delly_mappability_gindex.file_path in delly_files


def test_get_delly_exclusion_converted_file(
    hg_references_model: HgReferencesModel, delly_exclusion_converted_file: Path
):
    """Test get Delly exclusion converted file."""

    # GIVEN a human genome references model and a delly exclusion converted file

    # WHEN getting the Delly exclusion converted file
    converted_file: str = hg_references_model.get_delly_exclusion_converted_file()

    # THEN the returned file should match the expected one
    assert converted_file == delly_exclusion_converted_file.as_posix()


def test_get_gnomad_files(hg_references_model: HgReferencesModel):
    """Test get gnomad reference files."""

    # GIVEN a human genome references model

    # WHEN getting the gnomad reference files
    gnomad_files: List[str] = hg_references_model.get_gnomad_files()

    # THEN the gnomad files should be returned
    assert len(gnomad_files) == 2
    assert hg_references_model.gnomad_variant.file_path in gnomad_files
    assert hg_references_model.gnomad_variant_index.file_path in gnomad_files


def test_get_1k_genome_files(hg_references_model: HgReferencesModel):
    """Test get 1000 Genome related files."""

    # GIVEN a human genome references model

    # WHEN getting the 1k genome files
    genome_files: List[str] = hg_references_model.get_1k_genome_files()

    # THEN the 1k genome files should be returned
    assert len(genome_files) == 4
    assert (
        hg_references_model.known_indel_1kg.file_path + "." + FileType.GZ
        in genome_files
    )
    assert hg_references_model.mills_1kg.file_path + "." + FileType.GZ in genome_files
    assert hg_references_model.hc_vcf_1kg.file_path + "." + FileType.GZ in genome_files
    assert hg_references_model.vcf_1kg.file_path + "." + FileType.GZ in genome_files


def test_cache_analysis_model(cache_analysis_model_data: Dict[str, str]):
    """Test cache analysis model initialisation."""

    # GIVEN an input for the cache analysis model

    # WHEN initialising the model
    model: CacheAnalysisModel = CacheAnalysisModel(**cache_analysis_model_data)

    # THEN the model should have been correctly built
    assert model.dict() == cache_analysis_model_data


def test_cache_analysis_model_empty():
    """Test ache analysis model for an empty input."""

    # GIVEN no input for the cache analysis model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        CacheAnalysisModel()


def test_cache_config_model(
    cache_config_model_data: Dict[str, Any], cache_config_model: CacheConfigModel
):
    """Test cache config model initialisation."""

    # GIVEN an input for the cache config model and a mocked reference model

    # WHEN initialising the model
    model: CacheConfigModel = CacheConfigModel(**cache_config_model_data)

    # THEN the model should have been correctly built
    assert model == cache_config_model


def test_cache_config_model_empty():
    """Test cache config model for an empty input."""

    # GIVEN no input for the cache config model

    # WHEN initialising the model
    with pytest.raises(ValidationError):
        # THEN an empty model should raise a ValidationError
        CacheConfigModel()


def test_cache_config_model_empty_file_path(cache_config_model_data: Dict[str, dict]):
    """Test cache config model reference validation method and file path assignment."""

    # GIVEN a cache config model data with empty file paths

    # WHEN initialising the model
    model: CacheConfigModel = CacheConfigModel(**cache_config_model_data)

    # THEN the file paths should have been assigned
    for reference in model.references:
        assert reference[1].file_path


def test_cache_config_model_empty_cosmic_key(
    cache_config_model_data: Dict[str, dict], cosmic_key: str
):
    """Test cache config model reference validation method and cosmic key assignment."""

    # GIVEN a cache config model data with empty cosmic keys

    # WHEN initialising the model
    model: CacheConfigModel = CacheConfigModel(**cache_config_model_data)

    # THEN a cosmic key should only have been assigned to a cosmic reference file
    for reference in model.references:
        if reference[0] == "cosmic":
            assert reference[1].secret == cosmic_key
            continue
        assert reference[1].secret is None


def test_get_grch_version(cache_config_model: CacheConfigModel):
    """Test extraction of the GRCH format version having a specific genome version."""

    # GIVEN a cache config model

    # WHEN getting the GRCH version
    grch_version: GRCHVersion = cache_config_model.get_grch_version()

    # THEN a correct GRCH format version should be returned
    assert grch_version == GRCHVersion.GRCH37


def test_get_reference_paths(cache_config_model: CacheConfigModel):
    """Test reference path extraction."""

    # GIVEN a cache config model

    # WHEN extracting the list of reference paths
    reference_paths: List[str] = cache_config_model.get_reference_paths()

    # THEN a complete list of reference path should be returned
    assert reference_paths == [
        reference[1].file_path for reference in cache_config_model.references
    ]


def test_get_reference_by_path(cache_config_model: CacheConfigModel):
    """Test reference extraction given its path."""

    # GIVEN a cache config model

    # WHEN getting the reference genome by path
    reference_genome: ReferenceUrlModel = cache_config_model.get_reference_by_path(
        reference_path=cache_config_model.references.reference_genome.file_path
    )

    # THEN the correct reference should be returned
    assert reference_genome == cache_config_model.references.reference_genome


def test_get_reference_by_path_error(
    cache_config_model: CacheConfigModel, invalid_file: Path, caplog: LogCaptureFixture
):
    """Test reference extraction given an invalid path."""

    # GIVEN a cache config model

    # WHEN getting the reference genome by path
    with pytest.raises(BalsamicError):
        cache_config_model.get_reference_by_path(reference_path=invalid_file.as_posix())

    # THEN a Balsamic error should be returned
    assert (
        f"No reference with the provided reference path {invalid_file.as_posix()}"
        in caplog.text
    )


def test_get_reference_paths_by_file_type_and_compression(
    cache_config_model: CacheConfigModel,
):
    """Test reference path extraction by file type and compression."""

    # GIVEN a cache config model

    # WHEN extracting the reference paths by file type and compression status
    reference_paths: List[
        str
    ] = cache_config_model.get_reference_paths_by_file_type_and_compression(
        file_type=FileType.FASTA, compression=True
    )

    # THEN the expected reference path should be returned
    assert reference_paths == [cache_config_model.references.reference_genome.file_path]


def test_get_reference_paths_by_file_type(cache_config_model: CacheConfigModel):
    """Test reference path extraction by file type."""

    # GIVEN a cache config model

    # WHEN extracting the reference paths by file type
    reference_paths: List[str] = cache_config_model.get_reference_paths_by_file_type(
        file_type=FileType.FASTA
    )

    # THEN the TXT file should be returned
    assert reference_paths == [cache_config_model.references.reference_genome.file_path]


def test_get_reference_paths_by_compression(cache_config_model: CacheConfigModel):
    """Test reference path extraction by compression."""

    # GIVEN a cache config model

    # WHEN extracting the reference paths by compression status
    reference_paths: List[str] = cache_config_model.get_reference_paths_by_compression(
        compression=True
    )

    # THEN the expected reference path should be returned
    assert len(reference_paths) == 11
    for reference in [
        cache_config_model.references.ascat_gc_correction.file_path,
        cache_config_model.references.clinvar.file_path,
        cache_config_model.references.cosmic.file_path,
        cache_config_model.references.dbsnp.file_path,
        cache_config_model.references.hc_vcf_1kg.file_path,
        cache_config_model.references.known_indel_1kg.file_path,
        cache_config_model.references.mills_1kg.file_path,
        cache_config_model.references.reference_genome.file_path,
        cache_config_model.references.refgene_txt.file_path,
        cache_config_model.references.somalier_sites.file_path,
        cache_config_model.references.vcf_1kg.file_path,
    ]:
        assert reference in reference_paths


def test_get_compressed_indexed_vcfs(cache_config_model: CacheConfigModel):
    """Test get compressed indexed VCFs."""

    # GIVEN a cache config model

    # WHEN retrieving the compressed and indexed VCFs
    compressed_indexed_vcfs: List[
        str
    ] = cache_config_model.get_compressed_indexed_vcfs()

    # THEN the indexed VCFs should be returned
    assert len(compressed_indexed_vcfs) == 8
    for reference in [
        cache_config_model.references.dbsnp.file_path,
        cache_config_model.references.vcf_1kg.file_path,
        cache_config_model.references.known_indel_1kg.file_path,
        cache_config_model.references.mills_1kg.file_path,
        cache_config_model.references.clinvar.file_path,
        cache_config_model.references.somalier_sites.file_path,
        cache_config_model.references.hc_vcf_1kg.file_path,
        cache_config_model.references.cosmic.file_path,
    ]:
        assert (
            reference + "." + FileType.GZ + "." + FileType.TBI
            in compressed_indexed_vcfs
        )


def test_get_container_output_paths(
    cache_config_model: CacheConfigModel, tmp_path: Path
):
    """Test retrieval of the containers output paths."""

    # GIVEN a cache config model

    # WHEN getting the list of container paths
    container_paths: List[str] = cache_config_model.get_container_output_paths()

    # THEN all the container paths should be returned
    assert len(container_paths) == len(set(DockerContainers))
    for container in set(DockerContainers):
        assert Path(tmp_path, container + "." + FileType.SIF)


def test_get_reference_output_paths(cache_config_model: CacheConfigModel):
    """Test get reference list to be downloaded."""

    # GIVEN a cache config model

    # WHEN retrieving the reference output paths
    reference_output_paths: List[str] = cache_config_model.get_reference_output_paths()

    # THEN all the reference paths should be returned
    assert len(reference_output_paths) == 44


def test_get_canfam_analysis_references(
    cache_config_model: CacheConfigModel,
    analysis_references_model_data: Dict[str, Path],
):
    """Test analysis references retrieval to be used for Balsamic canine analyses."""

    # GIVEN a canine cache config model
    cache_config_model.genome_version = GenomeVersion.CanFam3

    # WHEN getting the analysis references
    analysis_references: CanFamAnalysisReferencesModel = (
        cache_config_model.get_analysis_references()
    )

    # THEN the retrieved analysis references should match the mocked one
    assert type(analysis_references) is CanFamAnalysisReferencesModel
    assert analysis_references.dict() == analysis_references_model_data
