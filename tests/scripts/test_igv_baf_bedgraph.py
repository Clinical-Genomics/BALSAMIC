from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pytest
from click.testing import CliRunner

import BALSAMIC.assets.scripts.igv_baf_bedgraph as m


@dataclass
class FakeVariant:
    CHROM: str
    POS: int
    REF: str
    ALT: list[str]
    _formats: dict[str, Any]
    gt_depths: Optional[list[int]] = None

    def format(self, key: str) -> Any:
        if key not in self._formats:
            raise KeyError(key)
        return self._formats[key]


@pytest.mark.parametrize(
    "alt,expected",
    [
        (None, True),
        (".", True),
        ("*", True),
        ("<DEL>", True),
        ("<NON_REF>", True),
        ("A", False),
        ("C", False),
        ("C,G", False),  # not symbolic; caller splits ALTs upstream anyway
    ],
)
def test_is_symbolic_or_missing_alt(alt, expected):
    assert m._is_symbolic_or_missing_alt(alt) is expected


@pytest.mark.parametrize(
    "ad_array,expected",
    [
        (np.array([[10, 5]], dtype=int), (10, 5)),
        (np.array([[10, 5, 2]], dtype=int), (10, 7)),  # multi-ALT sum
        (np.array([[0, 1]], dtype=int), (0, 1)),
        (
            np.array([[-2147483648, 5]], dtype=int).reshape(1, 2),
            None,
        ),  # missing sentinel
        (np.array([[10]], dtype=int), None),  # too few alleles
        (np.array([], dtype=int), None),  # wrong shape
    ],
)
def test_sample_ad_dp(ad_array, expected):
    v = FakeVariant(
        CHROM="1",
        POS=100,
        REF="A",
        ALT=["C"],
        _formats={"AD": ad_array},
        gt_depths=[20],
    )
    assert m._sample_ad_dp(v) == expected


@pytest.mark.parametrize(
    "dp_array,gt_depths,expected",
    [
        (np.array([[20]], dtype=int), [99], 20),  # FORMAT/DP preferred
        (
            np.array([[0]], dtype=int),
            [30],
            30,
        ),  # FORMAT/DP not sane -> fallback gt_depths
        (None, [30], 30),  # no FORMAT/DP -> fallback gt_depths
        (None, None, None),  # nothing available
    ],
)
def test_sample_dp(dp_array, gt_depths, expected):
    formats = {}
    if dp_array is not None:
        formats["DP"] = dp_array
    v = FakeVariant(
        CHROM="1",
        POS=100,
        REF="A",
        ALT=["C"],
        _formats=formats,
        gt_depths=gt_depths,
    )
    assert m._sample_dp(v) == expected


def test_variant_to_bedgraph_record_happy_path_with_dp():
    v = FakeVariant(
        CHROM="1",
        POS=100,
        REF="A",
        ALT=["C"],
        _formats={
            "AD": np.array([[10, 5]], dtype=int),
            "DP": np.array([[20]], dtype=int),
        },
        gt_depths=[20],
    )
    rec = m._variant_to_bedgraph_record(v)
    assert rec is not None
    assert rec.chrom == "1"
    assert rec.start == 99
    assert rec.end == 100
    assert rec.value == pytest.approx(5 / 20)


def test_variant_to_bedgraph_record_dp_fallback_to_ad_sum():
    v = FakeVariant(
        CHROM="1",
        POS=100,
        REF="A",
        ALT=["C"],
        _formats={"AD": np.array([[10, 5]], dtype=int)},  # no DP
        gt_depths=None,
    )
    rec = m._variant_to_bedgraph_record(v)
    assert rec is not None
    assert rec.value == pytest.approx(5 / 15)


def test_variant_to_bedgraph_record_span_for_indel_ref_length():
    v = FakeVariant(
        CHROM="1",
        POS=100,
        REF="ATG",
        ALT=["A"],
        _formats={
            "AD": np.array([[10, 5]], dtype=int),
            "DP": np.array([[20]], dtype=int),
        },
        gt_depths=[20],
    )
    rec = m._variant_to_bedgraph_record(v)
    assert rec is not None
    assert rec.start == 99
    assert rec.end == 102  # span = len(REF)=3


@pytest.mark.parametrize(
    "variant",
    [
        # no ALT list
        FakeVariant(
            "1",
            100,
            "A",
            [],
            {"AD": np.array([[10, 5]], dtype=int), "DP": np.array([[20]], dtype=int)},
            [20],
        ),
        # missing/no-call/symbolic ALT
        FakeVariant(
            "1",
            100,
            "A",
            ["."],
            {"AD": np.array([[10, 5]], dtype=int), "DP": np.array([[20]], dtype=int)},
            [20],
        ),
        FakeVariant(
            "1",
            100,
            "A",
            ["*"],
            {"AD": np.array([[10, 5]], dtype=int), "DP": np.array([[20]], dtype=int)},
            [20],
        ),
        FakeVariant(
            "1",
            100,
            "A",
            ["<DEL>"],
            {"AD": np.array([[10, 5]], dtype=int), "DP": np.array([[20]], dtype=int)},
            [20],
        ),
        # missing AD
        FakeVariant("1", 100, "A", ["C"], {"DP": np.array([[20]], dtype=int)}, [20]),
        # dp <= 0 and AD sum <= 0
        FakeVariant(
            "1",
            100,
            "A",
            ["C"],
            {"AD": np.array([[0, 0]], dtype=int), "DP": np.array([[0]], dtype=int)},
            [0],
        ),
        # AF trivial 0 or 1
        FakeVariant(
            "1",
            100,
            "A",
            ["C"],
            {"AD": np.array([[10, 0]], dtype=int), "DP": np.array([[10]], dtype=int)},
            [10],
        ),
        FakeVariant(
            "1",
            100,
            "A",
            ["C"],
            {"AD": np.array([[0, 10]], dtype=int), "DP": np.array([[10]], dtype=int)},
            [10],
        ),
    ],
)
def test_variant_to_bedgraph_record_skips_invalid(variant):
    assert m._variant_to_bedgraph_record(variant) is None


def test_variant_to_bedgraph_record_accepts_multi_alt_by_summing_ad():
    v = FakeVariant(
        CHROM="1",
        POS=100,
        REF="A",
        ALT=["C", "G"],
        _formats={
            "AD": np.array([[10, 5, 2]], dtype=int),
            "DP": np.array([[20]], dtype=int),
        },
        gt_depths=[20],
    )
    rec = m._variant_to_bedgraph_record(v)
    assert rec is not None
    assert rec.value == pytest.approx((5 + 2) / 20)


# -----------------------------------------------------------------------------
# Integration tests (real cyvcf2 reading)
# -----------------------------------------------------------------------------


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def test_convert_vcf_to_bedgraph_writes_expected_rows(tmp_path: Path):
    vcf = tmp_path / "in.vcf"
    out = tmp_path / "out.bedgraph"

    write_text(
        vcf,
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##source=test",
                "##contig=<ID=1>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
                "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,5:20",
                "1\t101\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/0:10,0:10",
                "1\t102\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t1/1:0,10:10",
                "1\t103\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:20",
                "1\t104\t.\tA\tC\t.\t.\t.\tGT:AD\t0/1:10,5",
            ]
        )
        + "\n",
    )

    n = m.convert_vcf_to_bedgraph(
        vcf_path=str(vcf), bedgraph_path=str(out), track_name=None
    )
    assert n == 2

    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert lines == [
        "1\t99\t100\t0.250000",
        "1\t103\t104\t0.333333",
    ]


def test_convert_vcf_to_bedgraph_writes_track_header(tmp_path: Path):
    vcf = tmp_path / "in.vcf"
    out = tmp_path / "out.bedgraph"

    write_text(
        vcf,
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##source=test",
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
    text = out.read_text(encoding="utf-8").splitlines()
    assert text[0] == 'track type=bedGraph name="MyTrack"'
    assert text[1] == "1\t99\t100\t0.250000"


def test_convert_vcf_to_bedgraph_dash_output_does_not_close_stdout(
    tmp_path: Path, capsys
):
    # open_output('-') must NOT close sys.stdout.
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

    n = m.convert_vcf_to_bedgraph(vcf_path=str(vcf), bedgraph_path="-", track_name=None)
    assert n == 1

    # If stdout was closed, this would raise ValueError: I/O operation on closed file.
    print("still-open")

    captured = capsys.readouterr()
    assert "1\t99\t100\t0.250000" in captured.out
    assert "still-open" in captured.out


def test_cli_end_to_end(tmp_path: Path):
    runner = CliRunner()
    vcf = tmp_path / "in.vcf"
    out = tmp_path / "out.bedgraph"

    write_text(
        vcf,
        "\n".join(
            [
                "##fileformat=VCFv4.2",
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
