import vcfpy
import pytest
from click.testing import CliRunner

import BALSAMIC.assets.scripts.rescue_snvs as rv
from pathlib import Path


@pytest.fixture
def vcf_file(tmp_path: Path) -> Path:
    vcf_file = tmp_path / "a.vcf"
    vcf_file.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinSig">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\trs1\tA\tT\t.\tPASS\t.\n"
        "1\t200\trs2\tG\tC\t.\tLowQ\tCLNSIG=Benign\n"
        "2\t300\trs3\tT\tC,G\t.\tq10\tONC=Oncogenic;CLNSIG=Likely_pathogenic\n"
    )
    return vcf_file


@pytest.fixture
def vcf_header(vcf_file: Path) -> vcfpy.Header:
    reader = vcfpy.Reader.from_path(vcf_file)
    return reader.header


@pytest.fixture
def vcf_records(vcf_file: Path) -> list[vcfpy.Record]:
    reader = vcfpy.Reader.from_path(vcf_file)
    return list(reader)


@pytest.fixture
def vcf_record(vcf_records: list[vcfpy.Record]) -> vcfpy.Record:
    return vcf_records[0]


@pytest.fixture
def vcf_record_rescued(vcf_records: list[vcfpy.Record]) -> vcfpy.Record:
    return vcf_records[2]


@pytest.fixture
def vcf_record_not_rescued(vcf_records: list[vcfpy.Record]) -> vcfpy.Record:
    return vcf_records[1]


def test_determine_clinvar_reasons_pathogenic(vcf_record_rescued):
    reasons = rv.determine_clinvar_reasons(vcf_record_rescued)
    assert rv.RescueReasons.CLINVAR_LIKELY_PATH not in reasons
    assert rv.RescueReasons.CLINVAR_ONC in reasons


def test_record_in_rescue_list_true(vcf_record):
    key = ("1", 100, "A", "T")
    assert rv.record_in_rescue_list(vcf_record, {key}) is True


def test_record_in_rescue_list_false(vcf_record):
    key = ("2", 200, "G", "C")
    assert rv.record_in_rescue_list(vcf_record, {key}) is False


def test_process_record_rescued(vcf_record_rescued):
    key = ("2", 300, "T", "G")
    record = rv.process_record(vcf_record_rescued, {key})
    assert record.FILTER == ["PASS"]
    assert "RescueFilters" in record.INFO
    assert "RescueStatus" in record.INFO
    assert rv.RescueReasons.CLINVAR_ONC in record.INFO["RescueStatus"]
    assert rv.RescueReasons.CLINVAR_LIKELY_PATH not in record.INFO["RescueStatus"]
    assert rv.RescueReasons.RESCUE_LIST in record.INFO["RescueStatus"]


def test_process_record_not_rescued(vcf_record_not_rescued):
    record = rv.process_record(vcf_record_not_rescued, None)
    assert record.FILTER == ["LowQ"]


def test_update_headers(vcf_file):
    reader = vcfpy.Reader.from_path(vcf_file)
    rv.update_headers(reader)
    assert "RescueStatus" in {h for h in reader.header.info_ids()}


def test_process_vcf_integration(vcf_file, vcf_header, vcf_record, tmp_path):
    out_vcf = tmp_path / "out.vcf"
    rv.process_vcf(str(vcf_file), str(out_vcf), {("1", 100, "A", "T")})
    reader = vcfpy.Reader.from_path(str(out_vcf))
    out_records = list(reader)
    assert out_records[0].FILTER == ["PASS"]
    assert "RescueStatus" in out_records[2].INFO


def test_cli_runs(tmp_path, vcf_header, vcf_record):
    in_vcf = tmp_path / "in.vcf"
    out_vcf = tmp_path / "out.vcf"
    with vcfpy.Writer.from_path(str(in_vcf), vcf_header) as writer:
        writer.write_record(vcf_record)

    runner = CliRunner()
    result = runner.invoke(rv.cli, ["--vcf", str(in_vcf), "-o", str(out_vcf)])
    assert result.exit_code == 0
    assert out_vcf.exists()
