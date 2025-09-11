import vcfpy
import pytest
from click.testing import CliRunner

import BALSAMIC.assets.scripts.rescue_snvs as rv


@pytest.fixture
def vcf_header():
    return vcfpy.Header(
        lines=[
            vcfpy.HeaderLine("fileformat", "VCFv4.2"),
            vcfpy.ContigHeaderLine("1", {}),
        ]
    )


@pytest.fixture
def sample_record(vcf_header):
    return vcfpy.Record(
        CHROM="1",
        POS=100,
        ID=["."],
        REF="A",
        ALT=[vcfpy.Substitution("SNV", "T")],
        QUAL=None,
        FILTER=["q10"],
        INFO={"CLNSIG": "Pathogenic", "ONC": ["oncogenic"]},
        FORMAT=[],
        calls=[],
    )


def test_determine_clinvar_reasons_pathogenic(sample_record):
    reasons = rv.determine_clinvar_reasons(sample_record)
    assert rv.RescueReasons.CLINVAR_PATH in reasons
    assert rv.RescueReasons.CLINVAR_ONC in reasons


def test_record_in_rescue_list_true(sample_record):
    key = ("1", 100, "A", "T")
    assert rv.record_in_rescue_list(sample_record, {key}) is True


def test_record_in_rescue_list_false(sample_record):
    key = ("2", 200, "G", "C")
    assert rv.record_in_rescue_list(sample_record, {key}) is False


def test_process_record_rescued(sample_record):
    key = ("1", 100, "A", "T")
    record = rv.process_record(sample_record, {key})
    assert record.FILTER == ["PASS"]
    assert "RescueFilters" in record.INFO
    assert "RescueStatus" in record.INFO
    assert rv.RescueReasons.RESCUE_LIST in record.INFO["RescueStatus"]


def test_process_record_not_rescued(sample_record):
    sample_record.INFO = {"other": "value1"}
    record = rv.process_record(sample_record, None)
    assert record.FILTER == ["q10"]


def test_update_headers(vcf_header):
    rv.update_headers(vcfpy.Reader.from_path)
    reader = vcfpy.Reader.from_string(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    rv.update_headers(reader)
    assert "RescueStatus" in {h.id for h in reader.header.info}


def test_process_vcf_integration(vcf_header, sample_record, tmp_path):
    # Write input VCF
    in_vcf = tmp_path / "in.vcf"
    out_vcf = tmp_path / "out.vcf"
    with vcfpy.Writer.from_path(str(in_vcf), vcf_header) as writer:
        writer.write_record(sample_record)

    # Run processing
    rv.process_vcf(str(in_vcf), str(out_vcf), {("1", 100, "A", "T")})

    # Read output VCF
    reader = vcfpy.Reader.from_path(str(out_vcf))
    out_record = next(reader)
    assert out_record.FILTER == ["PASS"]
    assert "RescueStatus" in out_record.INFO


def test_cli_runs(tmp_path, vcf_header, sample_record):
    in_vcf = tmp_path / "in.vcf"
    out_vcf = tmp_path / "out.vcf"
    with vcfpy.Writer.from_path(str(in_vcf), vcf_header) as writer:
        writer.write_record(sample_record)

    runner = CliRunner()
    result = runner.invoke(rv.cli, ["--vcf", str(in_vcf), "-o", str(out_vcf)])
    assert result.exit_code == 0
    assert out_vcf.exists()
