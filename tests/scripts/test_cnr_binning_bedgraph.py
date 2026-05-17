import io

import pandas as pd
from click.testing import CliRunner

import BALSAMIC.assets.scripts.igv_cnr_binning_bedgraph as m


def test_parse_skips_headers():
    assert m._parse_denoisedcr_row("@SQ\tSN:1\tLN:248956422") is None
    assert m._parse_denoisedcr_row("CONTIG\tSTART\tEND\tVALUE") is None


def test_parse_requires_4_fields():
    assert m._parse_denoisedcr_row("1\t10\t20") is None
    assert m._parse_denoisedcr_row("") is None
    assert m._parse_denoisedcr_row("   ") is None


def test_parse_invalid_numbers():
    assert m._parse_denoisedcr_row("1\tA\t20\t0.1") is None
    assert m._parse_denoisedcr_row("1\t10\tB\t0.1") is None
    assert m._parse_denoisedcr_row("1\t10\t20\tNOPE") is None


def test_parse_valid_row_whitespace():
    # also checks it uses first 4 fields even if more are present
    assert m._parse_denoisedcr_row("1  10  20  0.5  EXTRA") == ("1", 10, 20, 0.5)


def test_flush_buffer_noop_when_empty_or_no_chrom():
    out = []
    m._flush_buffer(out, "1", [], bin_size=2)
    assert out == []

    out = []
    m._flush_buffer(out, None, [(1, 2, 0.1)], bin_size=2)
    assert out == []


def test_flush_buffer_chunks_and_means():
    out = []
    buf = [
        (0, 10, 1.0),
        (10, 20, 3.0),
        (20, 30, 5.0),
        (30, 40, 7.0),
        (40, 50, 9.0),
    ]
    m._flush_buffer(out, "1", buf, bin_size=2)

    # windows: [0-20 mean 2], [20-40 mean 6], [40-50 mean 9]
    assert out == [
        ("1", 0, 20, 2.0),
        ("1", 20, 40, 6.0),
        ("1", 40, 50, 9.0),
    ]


def test_bin_denoised_segments_multichrom_and_skips_invalid():
    text = "\n".join(
        [
            "@metadata line",
            "CONTIG START END VALUE",
            "1\t0\t10\t1.0",
            "1\t10\t20\t3.0",
            "1\t20\t30\t5.0",
            "BADLINE",
            "1\t30\t40\t7.0",
            "2\t0\t10\t2.0",
            "2\t10\t20\t4.0",
            "2\t20\t30\t6.0",
            "2\t30\t40\t8.0",
            "2\t40\t50\t10.0",
            "2\t50\t60\tNOPE",  # skipped
            "",
        ]
    )
    infile = io.StringIO(text)
    outfile = io.StringIO()

    m.bin_denoised_segments(infile=infile, outfile=outfile, bin_size=2)

    outfile.seek(0)
    df = pd.read_csv(
        outfile,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "value"],
        dtype={"chrom": "string"},
    )

    expected = pd.DataFrame(
        [
            ("1", 0, 20, 2.0),
            ("1", 20, 40, 6.0),
            ("2", 0, 20, 3.0),
            ("2", 20, 40, 7.0),
            ("2", 40, 50, 10.0),
        ],
        columns=["chrom", "start", "end", "value"],
    )

    pd.testing.assert_frame_equal(df, expected, check_dtype=False)


def test_bin_denoised_segments_writes_track_header_when_requested():
    text = (
        "\n".join(
            [
                "1\t0\t10\t1.0",
                "1\t10\t20\t3.0",
            ]
        )
        + "\n"
    )
    infile = io.StringIO(text)
    outfile = io.StringIO()

    m.bin_denoised_segments(
        infile=infile, outfile=outfile, bin_size=2, track_name="MyTrack"
    )

    outfile.seek(0)
    lines = outfile.read().splitlines()
    assert lines[0] == 'track type=bedGraph name="MyTrack"'
    assert lines[1:] == ["1\t0\t20\t2.0"]


def test_bin_denoised_segments_empty_input_produces_empty_output():
    infile = io.StringIO("@hdr\nCONTIG\tSTART\tEND\tVALUE\n")
    outfile = io.StringIO()
    m.bin_denoised_segments(infile=infile, outfile=outfile, bin_size=3)

    outfile.seek(0)
    assert outfile.read().strip() == ""


def test_cli_happy_path_writes_output(tmp_path):
    runner = CliRunner()

    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    in_path.write_text(
        "\n".join(
            [
                "1\t0\t10\t1.0",
                "1\t10\t20\t3.0",
                "1\t20\t30\t5.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    res = runner.invoke(m.cli, [str(in_path), str(out_path), "-b", "2"])
    assert res.exit_code == 0, res.output

    df = pd.read_csv(
        out_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "value"],
        dtype={"chrom": "string"},
    )
    expected = pd.DataFrame(
        [
            ("1", 0, 20, 2.0),
            ("1", 20, 30, 5.0),
        ],
        columns=["chrom", "start", "end", "value"],
    )
    pd.testing.assert_frame_equal(df, expected, check_dtype=False)


def test_cli_happy_path_with_track_name_writes_header_and_output(tmp_path):
    runner = CliRunner()

    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    in_path.write_text(
        "\n".join(
            [
                "1\t0\t10\t1.0",
                "1\t10\t20\t3.0",
                "1\t20\t30\t5.0",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    res = runner.invoke(
        m.cli, [str(in_path), str(out_path), "-b", "2", "--track-name", "T"]
    )
    assert res.exit_code == 0, res.output

    lines = out_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == 'track type=bedGraph name="T"'

    df = pd.read_csv(
        io.StringIO("\n".join(lines[1:]) + "\n"),
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "value"],
        dtype={"chrom": "string"},
    )
    expected = pd.DataFrame(
        [
            ("1", 0, 20, 2.0),
            ("1", 20, 30, 5.0),
        ],
        columns=["chrom", "start", "end", "value"],
    )
    pd.testing.assert_frame_equal(df, expected, check_dtype=False)


def test_cli_rejects_nonpositive_bin_size(tmp_path):
    runner = CliRunner()
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    in_path.write_text("1\t0\t10\t1.0\n", encoding="utf-8")

    res = runner.invoke(m.cli, [str(in_path), str(out_path), "-b", "0"])
    assert res.exit_code != 0
    assert "positive integer" in res.output
