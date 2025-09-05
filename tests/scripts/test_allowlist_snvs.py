# tests/test_vcf_allowlist.py
from __future__ import annotations

import gzip
import io
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner


import BALSAMIC.assets.scripts.allowlist_snvs as m


# ---------------------------
# Fixtures / small utilities
# ---------------------------


@pytest.fixture
def tmp_text_vcf(tmp_path: Path) -> Path:
    p = tmp_path / "a.vcf"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinSig">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\trs1\tA\tT\t.\tPASS\t.\n"
        "1\t200\trs2\tG\tC\t.\tLowQ\tCLNSIG=Benign\n"
        "2\t300\trs3\tT\tC,G\t.\tq10\tONC;CLNSIG=Likely_pathogenic\n"
    )
    return p


@pytest.fixture
def tmp_gz_vcf(tmp_path: Path) -> Path:
    p = tmp_path / "b.vcf.gz"
    content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "3\t400\t.\tC\tA\t.\tPASS\t.\n"
        "4\t500\t.\tG\tT\t.\tPASS\t.\n"
    ).encode("utf-8")
    with gzip.open(p, "wb") as fh:
        fh.write(content)
    return p


# ---------------------------
# parse_info / format_info
# ---------------------------


def test_parse_info_empty_and_dot():
    assert m.parse_info("") == {}
    assert m.parse_info(".") == {}


def test_parse_info_mixed_kv_and_flags():
    s = "DP=12;SOMATIC;CLNSIG=Pathogenic;EMPTY="
    d = m.parse_info(s)
    assert d["DP"] == "12"
    assert d["CLNSIG"] == "Pathogenic"
    assert d["SOMATIC"] is None
    # split("=", 1) means value can be empty string
    assert d["EMPTY"] == ""


def test_format_info_sorted_and_roundtrip():
    d = {"ZKEY": "1", "AFLAG": None, "BKEY": "x"}
    out = m.format_info(d)
    # Sorted keys: AFLAG;BKEY=x;ZKEY=1
    assert out.split(";")[0] == "AFLAG"
    assert "BKEY=x" in out
    assert out.endswith("ZKEY=1")
    # roundtrip-ish (flags become None)
    d2 = m.parse_info(out)
    assert d2 == {"AFLAG": None, "BKEY": "x", "ZKEY": "1"}


# ---------------------------
# open_maybe_gzip / open_out_text
# ---------------------------


def test_open_maybe_gzip_plain_and_gz(tmp_text_vcf: Path, tmp_gz_vcf: Path):
    with m.open_maybe_gzip(tmp_text_vcf) as fh:
        first = fh.readline().strip()
        assert first.startswith("##fileformat")
    with m.open_maybe_gzip(tmp_gz_vcf) as fh:
        lines = [ln.strip() for ln in fh]
        assert lines[-1].startswith("4\t500")


def test_open_out_text_stdout_and_file(tmp_path: Path, monkeypatch):
    # Redirect stdout to a buffer so we don't close the real stdout
    buf = io.StringIO()
    monkeypatch.setattr(sys, "stdout", buf)

    # "-" returns sys.stdout
    out = m.open_out_text("-")
    assert out is buf
    out.write("hello\n")
    assert "hello" in buf.getvalue()

    # file path opens a real file
    p = tmp_path / "out.txt"
    with m.open_out_text(p) as fh:
        fh.write("x\n")
    assert p.read_text() == "x\n"


# ---------------------------
# build_allowlist_keyset
# ---------------------------


def test_build_allowlist_keyset_plain(tmp_text_vcf: Path):
    keys = m.build_allowlist_keyset(tmp_text_vcf)
    # expect tuples for each ALT allele
    assert ("1", 100, "A", "T") in keys
    assert ("1", 200, "G", "C") in keys
    assert ("2", 300, "T", "C") in keys
    assert ("2", 300, "T", "G") in keys


def test_build_allowlist_keyset_gz(tmp_gz_vcf: Path):
    keys = m.build_allowlist_keyset(tmp_gz_vcf)
    assert ("3", 400, "C", "A") in keys
    assert ("4", 500, "G", "T") in keys


# ---------------------------
# determine_clinvar_reasons
# ---------------------------


@pytest.mark.parametrize(
    "info,expected",
    [
        (
            {"ONC": None, "CLNSIG": "Pathogenic"},
            [m.CLINVAR_REASON_ONC, m.CLINVAR_REASON_PATH],
        ),
        ({"CLNSIG": "Likely_pathogenic"}, [m.CLINVAR_REASON_LIKELY_PATH]),
        ({"CLNSIG": "Benign"}, []),
        ({}, []),
        ({"ONC": None, "CLNSIG": "Benign"}, [m.CLINVAR_REASON_ONC]),
    ],
)
def test_determine_clinvar_reasons(info, expected):
    assert m.determine_clinvar_reasons(info) == expected


# ---------------------------
# ensure_info_headers
# ---------------------------


def test_ensure_info_headers_inserts_before_chrom():
    headers = [
        "##fileformat=VCFv4.3",
        "##source=test",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    out = m.ensure_info_headers(headers)
    # Two extra lines inserted before #CHROM
    idx = out.index("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    assert out[idx - 2] == m.INFO_ALLOWLISTED_FILTERS_HDR
    assert out[idx - 1] == m.INFO_ALLOWLIST_STATUS_HDR


def test_ensure_info_headers_noop_if_present():
    headers = [
        "##fileformat=VCFv4.3",
        m.INFO_ALLOWLISTED_FILTERS_HDR,
        m.INFO_ALLOWLIST_STATUS_HDR,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    out = m.ensure_info_headers(headers)
    assert out == headers


def test_ensure_info_headers_append_if_no_chrom():
    headers = ["##fileformat=VCFv4.3"]
    out = m.ensure_info_headers(headers)
    # If no #CHROM, insert at end
    assert out[-2] == m.INFO_ALLOWLISTED_FILTERS_HDR
    assert out[-1] == m.INFO_ALLOWLIST_STATUS_HDR


# ---------------------------
# _any_alt_in_allowlist
# ---------------------------


def test__any_alt_in_allowlist_true_false():
    keys = {
        ("1", 10, "A", "T"),
        ("1", 10, "A", "G"),
    }
    assert m._any_alt_in_allowlist(keys, "1", 10, "A", ["C", "G"]) is True
    assert m._any_alt_in_allowlist(keys, "1", 10, "A", ["C"]) is False


# ---------------------------
# _process_record_line
# ---------------------------


def test__process_record_line_no_change():
    line = "1\t100\trs1\tA\tT\t.\tPASS\t."
    out = m._process_record_line(line, allow_keys=None)
    assert out.strip() == line  # unchanged


def test__allowlist_sets_pass_and_moves_filters():
    # FILTER is LowQ -> should be moved to INFO AllowlistedFilters, FILTER set to PASS
    allow = {("1", 200, "G", "C")}
    line = "1\t200\trs2\tG\tC\t.\ttriallelic_site\tCLNSIG=Benign"
    out = m._process_record_line(line, allow_keys=allow).rstrip("\n")

    cols = out.split("\t")
    assert cols[6] == "PASS"
    info = m.parse_info(cols[7])
    assert info["AllowlistedFilters"] == "triallelic_site"
    # AllowlistStatus must include manual reason
    assert m.MANUAL_REASON in info["AllowlistStatus"]


def test__process_record_line_clinvar_reasons_only():
    # ONC + Likely_pathogenic
    line = "2\t300\trs3\tT\tC\t.\tHighOccurrenceFrq\tONC;CLNSIG=Likely_pathogenic"
    out = m._process_record_line(line, allow_keys=None).rstrip("\n")
    cols = out.split("\t")
    # FILTER should be set to PASS and original moved
    assert cols[6] == "PASS"
    info = m.parse_info(cols[7])
    assert info["AllowlistedFilters"] == "HighOccurrenceFrq"
    status = info["AllowlistStatus"].split("|")
    assert m.CLINVAR_REASON_ONC in status
    assert m.CLINVAR_REASON_LIKELY_PATH in status


def test__process_record_line_invalid_pos_or_short_cols_pass_through():
    # len(cols) < 8
    assert (
        m._process_record_line("1\t.\tx\tA\tT\t.\tPASS", None).strip()
        == "1\t.\tx\tA\tT\t.\tPASS"
    )
    # invalid pos
    assert (
        m._process_record_line("1\tNaN\tx\tA\tT\t.\tPASS\t.", None).strip()
        == "1\tNaN\tx\tA\tT\t.\tPASS\t."
    )


# ---------------------------
# process_vcf integration
# ---------------------------


def test_process_vcf_headers_written_once_and_records(tmp_path: Path):
    input_vcf = tmp_path / "in.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\t.\tA\tT\t.\tPASS\t.\n"
        "1\t200\t.\tG\tC\t.\tLowQ\tCLNSIG=Pathogenic\n"
    )
    allow = {("1", 100, "A", "T")}
    out_path = tmp_path / "out.vcf"
    with out_path.open("w", encoding="utf-8", newline="") as out_fh:
        m.process_vcf(allow, input_vcf, out_fh)

    out = out_path.read_text().splitlines()
    # Headers present once + inserted INFO lines
    assert out[0].startswith("##fileformat")
    assert any(l == m.INFO_ALLOWLISTED_FILTERS_HDR for l in out[:5])
    assert any(l == m.INFO_ALLOWLIST_STATUS_HDR for l in out[:5])
    # Two records
    rec1 = out[-2].split("\t")
    rec2 = out[-1].split("\t")

    # rec1 was manually allow-listed, FILTER may remain PASS, but AllowlistStatus must exist
    info1 = m.parse_info(rec1[7])
    assert "AllowlistStatus" in info1
    # rec2 had CLNSIG=Pathogenic and FILTER=LowQ -> FILTER becomes PASS and moved
    assert rec2[6] == "PASS"
    info2 = m.parse_info(rec2[7])
    assert info2["AllowlistedFilters"] == "LowQ"
    assert m.CLINVAR_REASON_PATH in (info2["AllowlistStatus"] or "")


def test_process_vcf_file_with_only_headers(tmp_path: Path):
    input_vcf = tmp_path / "only_headers.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    out_path = tmp_path / "out.vcf"
    with out_path.open("w", encoding="utf-8", newline="") as out_fh:
        m.process_vcf(None, input_vcf, out_fh)
    out = out_path.read_text().splitlines()
    # Headers should be present (plus inserted INFO lines)
    assert any(l == m.INFO_ALLOWLISTED_FILTERS_HDR for l in out)
    assert any(l == m.INFO_ALLOWLIST_STATUS_HDR for l in out)


# ---------------------------
# CLI smoke test (optional)
# ---------------------------


def test_cli_smoke(tmp_text_vcf: Path, tmp_path: Path):
    runner = CliRunner()
    out_path = tmp_path / "cli_out.vcf"
    res = runner.invoke(m.cli, ["--vcf", str(tmp_text_vcf), "-o", str(out_path)])
    assert res.exit_code == 0
    text = out_path.read_text()
    assert "##fileformat" in text
    # INFO headers injected
    assert m.INFO_ALLOWLISTED_FILTERS_HDR in text
    assert m.INFO_ALLOWLIST_STATUS_HDR in text
