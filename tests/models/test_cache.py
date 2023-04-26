"""Test reference cache models."""
from datetime import datetime
from pathlib import Path

from BALSAMIC.commands.init.utils import get_containers
from BALSAMIC.constants.cache import ContainerVersion, GenomeVersion
from BALSAMIC.constants.common import BIOINFO_TOOL_ENV
from BALSAMIC.constants.references import REFERENCE_FILES
from BALSAMIC.models.cache_models import CacheConfigModel


def test_cache_config_model(tmp_path):
    """"""

    config: CacheConfigModel = CacheConfigModel(
        analysis={"case_id": "reference.v1"},
        references_dir=tmp_path,
        variants_dir=tmp_path,
        genome_dir=tmp_path,
        vep_dir=tmp_path,
        containers_dir=tmp_path,
        genome_version=GenomeVersion.HG19,
        bioinfo_tools=BIOINFO_TOOL_ENV,
        containers=get_containers(ContainerVersion.RELEASE),
        references=REFERENCE_FILES[GenomeVersion.HG19],
        references_date=str(datetime.now),
    )
    config_path: Path = Path(config.references_dir, "config.json")
    json_obj = config.json(exclude_none=True)

    # with open(config_path, "w") as fn:
    #     json.dump(json.loads(json_obj), fn, indent=4)
    #
    # import os
    #
    # os.system(f"open {config.references_dir}")


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
