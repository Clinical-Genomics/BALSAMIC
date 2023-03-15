from pathlib import Path

import pytest


@pytest.fixture
def input_file():
    return "tests/test_data/vcf_tables/test_input.txt"


@pytest.fixture
def output_file(tmp_path):
    return Path(tmp_path, "test_createVCF_output.vcf.gz").as_posix()


@pytest.fixture
def reference_file():
    return "tests/test_data/vcf_tables/test_reference.vcf.gz"
