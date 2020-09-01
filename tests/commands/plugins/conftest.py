import pytest

@pytest.fixture
def input_file():
    return "tests/test_data/vcf_tables/test_input.txt"

@pytest.fixture
def output_file():
    return "tests/test_data/vcf_tables/test_createVCF_output.vcf.gz"

@pytest.fixture
def reference_file():
    return "tests/test_data/vcf_tables/test_reference.vcf.gz"
