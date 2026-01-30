import gzip
from pathlib import Path

import pytest
from click.testing import CliRunner

import BALSAMIC.assets.scripts.igv_baf_bedgraph as m

@pytest.mark.parametrize(
    "ad,expected",
    [
        ("10,5", (10, 5)),
        ("10,5,2", (10, 7)),          # multi-ALT sum
        ("0,1", (0, 1)),
        (".", None),
        ("", None),
        ("10", None),                 # not enough fields
        ("10,.", None),               # alt missing -> counts list has len 1
        ("A,2", None),                # non-numeric
        ("10,2,.", (10, 2)),          # ignore '.' after valid counts
    ],
)
def test_parse_ad_counts(ad, expected):
    assert m._parse_ad_counts(ad) == expected


@pytest.mark.parametrize(
    "ref_ad,alt_sum,dp_raw,expected",
    [
        (10, 5, "20", 20),
        (10, 5, ".", 15),
        (10, 5, None, 15),
        (10, 5, "NOPE", 15),
        (0, 0, "0", 0),
    ],
)
def test_compute_dp(ref_ad, alt_sum, dp_raw, expected):
    assert m._compute_dp(ref_ad, alt_sum, dp_raw) == expected


def test_parse_header_line_requires_sample_column():
    # 9 columns (no sample) -> should raise
    with pytest.raises(Exception) as ei:
        m._parse_header_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
    # It's a ClickException; we don't need to import click just for isinstance checks
    assert "Expected exactly one sample column" in str(ei.value)


def test_extract_record_happy_path_with_dp():
    parts = [
        "1", "100", ".", "A", "C", ".", ".", ".", "GT:AD:DP", "0/1:10,5:20"
    ]
    rec = m._extract_bedgraph_record(parts)
    assert rec is not None
    assert rec.chrom == "1"
    assert rec.start == 99
    assert rec.end == 100
    assert rec.value == pytest.approx(5 / 20)


def test_extract_record_dp_fallback_to_ad_sum():
    parts = [
        "1", "100", ".", "A", "C", ".", ".", ".", "GT:AD", "0/1:10,5"
    ]
    rec = m._extract_bedgraph_record(parts)
    assert rec is not None
    assert rec.value == pytest.approx(5 / 15)


def test_extract_record_span_for_indel_ref_length():
    # ref length 3 -> span=3, start=pos-1=99, end=102
    parts = [
        "1", "100", ".", "ATG", "A", ".", ".", ".", "GT:AD:DP", "0/1:10,5:20"
    ]
    rec = m._extract_bedgraph_record(parts)
    assert rec is not None
    assert rec.start == 99
    assert rec.end == 102


@pytest.mark.parametrize(
    "parts",
    [
        # too few columns
        ["1", "100", ".", "A", "C", ".", ".", "."],
        # symbolic / no-call ALT
        ["1", "100", ".", "A", ".", ".", ".", ".", "GT:AD:DP", "0/1:10,5:20"],
        ["1", "100", ".", "A", "*", ".", ".", ".", "GT:AD:DP", "0/1:10,5:20"],
        # missing AD
        ["1", "100", ".", "A", "C", ".", ".", ".", "GT:DP", "0/1:20"],
        ["1", "100", ".", "A", "C", ".", ".", ".", "GT:AD:DP", "0/1:.:20"],
        # dp <= 0
        ["1", "100", ".", "A", "C", ".", ".", ".", "GT:AD:DP", "0/1:0,0:0"],
        # AF trivial 0 or 1
        ["1", "100", ".", "A", "C", ".", ".", ".", "GT:AD:DP", "0/1:10,0:10"],
        ["1", "100", ".", "A", "C", ".", ".", ".", "GT:AD:DP", "0/1:0,10:10"],
        # malformed POS
        ["1", "NOPE", ".", "A", "C", ".", ".", ".", "GT:AD:DP", "0/1:10,5:20"],
    ],
)
def test_extract_record_skips_invalid(parts):
    assert m._extract_bedgraph_record(parts) is None


def test_extract_record_accepts_multi_alt_by_summing_ad():
    parts = [
        "1", "100", ".", "A", "C,G", ".", ".", ".", "GT:AD:DP", "0/1:10,5,2:20"
    ]
    rec = m._extract_bedgraph_record(parts)
    assert rec is not None
    assert rec.value == pytest.approx((5 + 2) / 20)


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
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
                # valid -> AF=5/20=0.25
                "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,5:20",
                # AF=0 -> skipped
                "1\t101\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/0:10,0:10",
                # AF=1 -> skipped
                "1\t102\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t1/1:0,10:10",
                # missing AD -> skipped
                "1\t103\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:20",
                # valid without DP -> DP fallback = 15, AF=5/15=0.333333...
                "1\t104\t.\tA\tC\t.\t.\t.\tGT:AD\t0/1:10,5",
            ]
        )
        + "\n",
    )

    n = m.convert_vcf_to_bedgraph(vcf_path=str(vcf), bedgraph_path=str(out), track_name=None)
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


def test_convert_raises_if_missing_chrom_header(tmp_path: Path):
    vcf = tmp_path / "noheader.vcf"
    out = tmp_path / "out.bedgraph"
    write_text(vcf, "##fileformat=VCFv4.2\n1\t100\t.\tA\tC\t.\t.\t.\tGT:AD\t0/1:10,5\n")

    with pytest.raises(Exception) as ei:
        m.convert_vcf_to_bedgraph(str(vcf), str(out))
    assert "missing #CHROM header" in str(ei.value)


def test_convert_raises_if_header_has_no_sample_column(tmp_path: Path):
    vcf = tmp_path / "badheader.vcf"
    out = tmp_path / "out.bedgraph"
    write_text(
        vcf,
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
                "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD\t0/1:10,5",
            ]
        )
        + "\n",
    )

    with pytest.raises(Exception) as ei:
        m.convert_vcf_to_bedgraph(str(vcf), str(out))
    assert "Expected exactly one sample column" in str(ei.value)


def test_convert_reads_gzipped_input(tmp_path: Path):
    vcfgz = tmp_path / "in.vcf.gz"
    out = tmp_path / "out.bedgraph"

    content = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n"
        "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,5:20\n"
    ).encode("utf-8")

    with gzip.open(vcfgz, "wb") as f:
        f.write(content)

    n = m.convert_vcf_to_bedgraph(str(vcfgz), str(out))
    assert n == 1
    assert out.read_text(encoding="utf-8").strip() == "1\t99\t100\t0.250000"


def test_cli_end_to_end(tmp_path: Path):
    runner = CliRunner()
    vcf = tmp_path / "in.vcf"
    out = tmp_path / "out.bedgraph"

    write_text(
        vcf,
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1",
                "1\t100\t.\tA\tC\t.\t.\t.\tGT:AD:DP\t0/1:10,5:20",
            ]
        )
        + "\n",
    )

    res = runner.invoke(m.cli, [str(vcf), str(out), "--track-name", "T"])
    assert res.exit_code == 0, res.output

    # writes progress to stderr; Click captures it in output by default in testing
    assert "Done. Wrote 1 rows" in res.output

    assert out.read_text(encoding="utf-8").splitlines()[0] == 'track type=bedGraph name="T"'
