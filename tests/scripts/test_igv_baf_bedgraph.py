from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from click.testing import CliRunner

import BALSAMIC.assets.scripts.igv_baf_bedgraph as m


@dataclass
class FakeVariant:
    CHROM: str
    POS: int
    _formats: dict[str, Any]

    def format(self, key: str) -> Any:
        return self._formats[key]


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


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
    vcf = tmp_path / "in.vcf"
    out = tmp_path / "out.bedgraph"

    write_text(
        vcf,
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
        + "\n",
    )

    n = m.convert_vcf_to_bedgraph(str(vcf), str(out), track_name=None)
    assert n == 2

    assert out.read_text(encoding="utf-8").splitlines() == [
        "1\t100\t101\t0.25",
        "1\t101\t102\t0.5",
    ]


def test_convert_vcf_to_bedgraph_writes_track_header(tmp_path: Path):
    vcf = tmp_path / "in.vcf"
    out = tmp_path / "out.bedgraph"

    write_text(
        vcf,
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
        + "\n",
    )

    n = m.convert_vcf_to_bedgraph(str(vcf), str(out), track_name="MyTrack")
    assert n == 1

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == 'track type=bedGraph name="MyTrack"'
    assert lines[1] == "1\t100\t101\t0.25"


def test_cli_end_to_end(tmp_path: Path):
    runner = CliRunner()
    vcf = tmp_path / "in.vcf"
    out = tmp_path / "out.bedgraph"

    write_text(
        vcf,
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
        + "\n",
    )

    res = runner.invoke(m.cli, [str(vcf), str(out), "--track-name", "T"])
    assert res.exit_code == 0, res.output
    assert "Done. Wrote 1 rows" in res.output
    assert (
        out.read_text(encoding="utf-8").splitlines()[0]
        == 'track type=bedGraph name="T"'
    )


def test_bedgraph_dash_output_does_not_close_stdout(tmp_path: Path, capsys):
    vcf = tmp_path / "in.vcf"

    write_text(
        vcf,
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
        + "\n",
    )

    n = m.convert_vcf_to_bedgraph(str(vcf), "-", track_name=None)
    assert n == 1

    print("still-open")
    captured = capsys.readouterr()
    assert "1\t100\t101\t0.25" in captured.out
    assert "still-open" in captured.out
