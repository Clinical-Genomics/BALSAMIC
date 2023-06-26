"""Test module for Balsamic cache references."""
import os
from pathlib import Path

import pytest
from pydantic import ValidationError

from BALSAMIC.models.cache_models import ReferenceMeta, ReferenceUrlsModel


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
        ReferenceUrlsModel.parse_obj(dummy_reference)
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
        ReferenceUrlsModel.parse_obj(dummy_reference)

        # THEN model raise error on validation
        assert "not a valid genome version" in excinfo.value
