

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from click.testing import CliRunner

CNV_REPORT_DIR = (
    Path(__file__).resolve().parents[3]
    / "BALSAMIC"
    / "assets"
    / "scripts"
    / "cnv_report"
)

sys.path.insert(0, str(CNV_REPORT_DIR))

import cnv_report_plotting
from cnv_qc_report import main

def test_generate_report_smoke(tmp_path, monkeypatch):
    """
    End-to-end smoke test:
    - writes tiny input files
    - runs CLI
    - verifies report and chromosome plot outputs are created

    This is intentionally a sanity check, not a detailed content test.
    """

    # ------------------------------------------------------------------
    # Monkeypatch VCF loading so the test does not depend on cyvcf2 parsing
    # details or a realistic VCF payload. We still provide a real file path
    # because the CLI requires it to exist.
    # ------------------------------------------------------------------
    monkeypatch.setattr(
        cnv_report_plotting,
        "load_vcf_with_vaf",
        lambda vcf_path, chr_order: pd.DataFrame(columns=["CHROM", "POS", "VAF"]),
    )

    case_id = "TEST_CASE"

    # ------------------------------------------------------------------
    # Minimal CNR input
    # Multi-gene bins are not needed here. Keep it simple.
    # ------------------------------------------------------------------
    cnr_path = tmp_path / "sample.cnr"
    cnr_path.write_text(
        "chromosome\tstart\tend\tgene\tlog2\tdepth\n"
        "1\t100\t150\tTP53\t0.20\t100\n"
        "1\t200\t250\tTP53\t0.25\t110\n"
        "1\t300\t350\tEGFR\t-0.30\t120\n"
        "1\t400\t450\tEGFR\t-0.20\t125\n"
        "2\t100\t150\tMYC\t0.40\t130\n"
        "2\t200\t250\tMYC\t0.35\t135\n"
    )

    # ------------------------------------------------------------------
    # Minimal called CNS
    # ------------------------------------------------------------------
    cns_path = tmp_path / "sample.cns"
    cns_path.write_text(
        "chromosome\tstart\tend\tlog2\tbaf\tcn\tcn1\tcn2\tdepth\n"
        "1\t100\t250\t0.25\t0.45\t3\t1\t2\t115\n"
        "1\t300\t450\t-0.25\t0.35\t1\t0\t1\t123\n"
        "2\t100\t250\t0.40\t0.50\t3\t1\t2\t132\n"
    )

    # ------------------------------------------------------------------
    # Minimal raw/init CNS
    # ------------------------------------------------------------------
    cns_init_path = tmp_path / "sample_init.cns"
    cns_init_path.write_text(
        "chromosome\tstart\tend\tlog2\n"
        "1\t100\t250\t0.20\n"
        "1\t300\t450\t-0.20\n"
        "2\t100\t250\t0.35\n"
    )

    # ------------------------------------------------------------------
    # Cytoband file (UCSC-like, no header)
    # chr, chromStart, chromEnd, name, gieStain
    # ------------------------------------------------------------------
    cytoband_path = tmp_path / "cytoband.tsv"
    cytoband_path.write_text(
        "chr1\t0\t300000000\tp36.33\tgneg\n"
        "chr2\t0\t300000000\tp25.3\tgneg\n"
    )

    # ------------------------------------------------------------------
    # Cancer genes file
    # ------------------------------------------------------------------
    cancer_genes_path = tmp_path / "cancer_genes.tsv"
    cancer_genes_path.write_text(
        "GENE\tOCCURRENCE\tSOURCE\n"
        "TP53\t10\tONCOKB\n"
        "EGFR\t8\tONCOKB\n"
        "MYC\t6\tONCOKB\n"
    )

    # ------------------------------------------------------------------
    # Minimal purity summary CSV
    # Safe to include so report summary code always has valid input.
    # ------------------------------------------------------------------
    purity_csv_path = tmp_path / "purity.csv"
    purity_csv_path.write_text(
        "Sampleid,Purity,Ploidy,Sex,Contamination,Flagged,Failed,Comment\n"
        "TEST_CASE,0.45,2.1,male,0.01,FALSE,FALSE,ok\n"
    )

    # ------------------------------------------------------------------
    # Required by click as existing file
    # Contents do not matter because load_vcf_with_vaf is monkeypatched
    # ------------------------------------------------------------------
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text("##fileformat=VCFv4.2\n")

    output_html = tmp_path / "report.html"

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--cnr",
            str(cnr_path),
            "--cns",
            str(cns_path),
            "--cns-init",
            str(cns_init_path),
            "--vcf",
            str(vcf_path),
            "--cytoband",
            str(cytoband_path),
            "--case-id",
            case_id,
            "--output-file",
            str(output_html),
            "--purity-csv",
            str(purity_csv_path),
            "--cancer-genes",
            str(cancer_genes_path),
            "--sex",
            "male",
            "--analysis-type",
            "single",
        ],
    )

    assert result.exit_code == 0, result.output

    assert output_html.exists()
    html = output_html.read_text(encoding="utf-8")
    assert case_id in html
    assert "CNV Report" in html

    chr_plots_dir = tmp_path / f"{case_id}_chr_plots"
    assert chr_plots_dir.exists()

    # At least the chromosomes present in the data should have plots
    assert (chr_plots_dir / "cnv_chr1_segments.png").exists()
    assert (chr_plots_dir / "cnv_chr2_segments.png").exists()