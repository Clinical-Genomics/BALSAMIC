# tests/test_allowlist_vcf.py
from __future__ import annotations

import io
from pathlib import Path

import pytest
import vcfpy
from click.testing import CliRunner

from BALSAMIC.assets.scripts.allowlist_snvs import (
    INFO_ALLOWLISTED_FILTERS_ID,
    INFO_ALLOWLIST_STATUS_ID,
    MANUAL_REASON,
    CLINVAR_REASON_ONC,
    CLINVAR_REASON_PATH,
    CLINVAR_REASON_LIKELY_PATH,
    build_allowlist_keyset,
    process_vcf,
    cli,
)

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1,length=1000000>\n"
    '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="ClinVar clinical significance">\n'
    '##INFO=<ID=ONC,Number=0,Type=Flag,Description="Oncogenic flag">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def _write_vcf(path: Path, body: str) -> None:
    path.write_text(VCF_HEADER + body, encoding="utf-8")


def _read_records_from_string(vcf_text: str):
    with io.StringIO(vcf_text) as fh:
        rdr = vcfpy.Reader.from_stream(fh)
        return list(rdr)


def _read_records_from_path(path: Path):
    with vcfpy.Reader.from_path(str(path)) as rdr:
        return list(rdr), rdr.header


@pytest.fixture
def tmp_vcf(tmp_path: Path) -> Path:
    return tmp_path / "in.vcf"


@pytest.fixture
def tmp_allow_vcf(tmp_path: Path) -> Path:
    return tmp_path / "allow.vcf"


@pytest.fixture
def out_vcf(tmp_path: Path) -> Path:
    return tmp_path / "out.vcf"


def test_build_allowlist_keyset_single_and_multiallelic(tmp_allow_vcf: Path):
    body = "1\t100\t.\tA\tT\t.\tPASS\t.\n" "1\t200\t.\tG\tA,C\t.\tPASS\t.\n"
    _write_vcf(tmp_allow_vcf, body)
    keys = build_allowlist_keyset(tmp_allow_vcf)
    assert ("1", 100, "A", "T") in keys
    assert ("1", 200, "G", "A") in keys
    assert ("1", 200, "G", "C") in keys
    assert len(keys) == 3


def test_headers_are_added_once_even_if_no_reasons(tmp_vcf: Path, out_vcf: Path):
    _write_vcf(tmp_vcf, "1\t300\t.\tA\tT\t.\tPASS\t.\n")
    # No allow-list, no ClinVar → still add header lines
    process_vcf(allow_keys=None, in_vcf=tmp_vcf, out_path=out_vcf)
    _, header = _read_records_from_path(out_vcf)

    info_ids = {line.id for line in header.get_lines("INFO")}
    assert INFO_ALLOWLISTED_FILTERS_ID in info_ids
    assert INFO_ALLOWLIST_STATUS_ID in info_ids

    # Ensure not duplicated
    assert (
        sum(
            1
            for line in header.get_lines("INFO")
            if line.id == INFO_ALLOWLISTED_FILTERS_ID
        )
        == 1
    )
    assert (
        sum(
            1
            for line in header.get_lines("INFO")
            if line.id == INFO_ALLOWLIST_STATUS_ID
        )
        == 1
    )


def test_manual_allowlist_moves_named_filters_and_sets_pass(
    tmp_vcf: Path, tmp_allow_vcf: Path, out_vcf: Path
):
    # Input record is filtered with LowQ → should move to INFO when allow-listed
    _write_vcf(tmp_vcf, "1\t100\t.\tA\tT\t.\tLowQ\t.\n")
    _write_vcf(tmp_allow_vcf, "1\t100\t.\tA\tT\t.\tPASS\t.\n")

    keys = build_allowlist_keyset(tmp_allow_vcf)
    process_vcf(keys, tmp_vcf, out_vcf)

    recs, _ = _read_records_from_path(out_vcf)
    rec = recs[0]
    assert rec.FILTER == ["PASS"]
    assert rec.INFO[INFO_ALLOWLISTED_FILTERS_ID] == "LowQ"
    # Manual reason is included
    status = rec.INFO[INFO_ALLOWLIST_STATUS_ID]
    assert MANUAL_REASON in status.split("|")


def test_clinvar_likely_pathogenic_with_pass_filter_does_not_add_allowlistedfilters(
    tmp_vcf: Path, out_vcf: Path
):
    # CLNSIG=Likely_pathogenic, FILTER PASS → no AllowlistedFilters
    body = "1\t250\t.\tT\tC\t.\tPASS\tCLNSIG=Likely_pathogenic\n"
    _write_vcf(tmp_vcf, body)

    process_vcf(None, tmp_vcf, out_vcf)

    recs, _ = _read_records_from_path(out_vcf)
    rec = recs[0]
    assert rec.FILTER == ["PASS"]
    assert INFO_ALLOWLISTED_FILTERS_ID not in rec.INFO
    assert rec.INFO[INFO_ALLOWLIST_STATUS_ID] == CLINVAR_REASON_LIKELY_PATH


def test_no_reasons_leaves_record_unchanged(tmp_vcf: Path, out_vcf: Path):
    # No allow-list and no ClinVar triggers → record unchanged
    body = "1\t300\t.\tA\tG\t.\tq10;LowQ\t.\n"
    _write_vcf(tmp_vcf, body)

    process_vcf(None, tmp_vcf, out_vcf)

    recs, _ = _read_records_from_path(out_vcf)
    rec = recs[0]
    assert rec.FILTER == ["q10", "LowQ"]
    assert INFO_ALLOWLISTED_FILTERS_ID not in rec.INFO
    assert INFO_ALLOWLIST_STATUS_ID not in rec.INFO
