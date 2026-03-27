from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from click.testing import CliRunner

import BALSAMIC.assets.scripts.igv_baf_bedgraph as m


@dataclass
class FakeVariant:
    CHROM: str
    POS: int
    _formats: dict[str, Any]

    def format(self, key: str) -> Any:
        return self._formats[key]


def test_variant_to_record_basic():
    v = FakeVariant(
        CHROM="1",
        POS=100,
        _formats={
            "AD": np.array([[10, 5]], dtype=int),
            "DP": np.array([[20]], dtype=int),
        },
    )
    line = m.variant_to_record(v)
    assert line == "1\t100\t101\t0.25\n"


def test_variant_to_record_returns_none_when_ad_is_none():
    v = FakeVariant(
        CHROM="1",
        POS=100,
        _formats={
            "AD": None,
            "DP": np.array([[20]], dtype=int),
        },
    )
    assert m.variant_to_record(v) is None


def test_convert_vcf_to_bedgraph_writes_expected_rows(tmp_path: Path):
    vcf_file: Path = tmp_path / "in.vcf"
    out_file: Path = tmp_path / "out.bedgraph"

    vcf_content = (
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=1>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
                "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,5:20",  # 0.25
                "1\t101\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:3,3:6",  # 0.5
            ]
        )
        + "\n"
    )
    vcf_file.write_text(vcf_content, encoding="utf-8")

    n_written, n_skipped = m.convert_vcf_to_bedgraph(
        str(vcf_file), str(out_file), track_name=None
    )
    assert n_written == 2
    assert n_skipped == 0

    assert out_file.read_text(encoding="utf-8").splitlines() == [
        "1\t100\t101\t0.25",
        "1\t101\t102\t0.5",
    ]


def test_convert_vcf_to_bedgraph_writes_track_header(tmp_path: Path):
    vcf_file: Path = tmp_path / "in.vcf"
    out_file: Path = tmp_path / "out.bedgraph"

    vcf_content = (
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=1>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
                "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,5:20",
            ]
        )
        + "\n"
    )
    vcf_file.write_text(vcf_content, encoding="utf-8")

    n_written, n_skipped = m.convert_vcf_to_bedgraph(
        str(vcf_file), str(out_file), track_name="MyTrack"
    )
    assert n_written == 1
    assert n_skipped == 0

    lines = out_file.read_text(encoding="utf-8").splitlines()
    assert lines[0] == 'track type=bedGraph name="MyTrack"'
    assert lines[1] == "1\t100\t101\t0.25"


def test_cli_end_to_end(tmp_path: Path):
    runner = CliRunner()
    vcf_file: Path = tmp_path / "in.vcf"
    out_file: Path = tmp_path / "out.bedgraph"

    vcf_content = (
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=1>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
                "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,5:20",
            ]
        )
        + "\n"
    )
    vcf_file.write_text(vcf_content, encoding="utf-8")

    res = runner.invoke(m.cli, [str(vcf_file), str(out_file), "--track-name", "T"])
    assert res.exit_code == 0, res.output
    assert "Done. Wrote 1 rows" in res.output
    assert "(0 variants skipped)" in res.output

    assert (
        out_file.read_text(encoding="utf-8").splitlines()[0]
        == 'track type=bedGraph name="T"'
    )
