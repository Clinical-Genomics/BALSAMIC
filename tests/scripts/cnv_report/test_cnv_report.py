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


def write_cnr_input(tmp_path: Path):
    """
    Create a realistic CNR file.
    """
    cnr_path = tmp_path / "sample.cnr"
    cnr_path.write_text(
        "chromosome\tstart\tend\tgene\tlog2\tdepth\nweight"
        "1\t100\t150\tTP53\t0.35\t100\n,0.1\n"
        "1\t200\t250\tTP53\t0.40\t101\n,0.1\n"
        "1\t300\t350\tTP53\t0.38\t102\n,0.1\n"
        "1\t400\t450\tTP53\t0.42\t103\n,0.1\n"
        "1\t500\t550\tTP53\t0.37\t104\n,0.1\n"
        "1\t600\t650\tTP53\t0.39\t105\n,0.1\n"
        "1\t700\t750\tTP53\t0.41\t106\n,0.1\n"
        "1\t800\t850\tTP53\t0.36\t107\n,0.1\n"
        "2\t100\t150\tMYC\t0.40\t130\n,0.1\n"
        "2\t200\t250\tMYC\t0.35\t135\n,0.1\n"
    )
    return cnr_path


def write_matching_pon_input(tmp_path: Path) -> Path:
    """
    Create a realistic PON .cnn file.
    """
    pon_path = tmp_path / "sample.cnn"

    pon_path.write_text(
        "chromosome\tstart\tend\tgene\tlog2\tdepth\tgc\trmask\tspread\n"
        "1\t100\t150\tTP53\t0.00\t100\t0.40\t0\t0.05\n"
        "1\t200\t250\tTP53\t0.00\t101\t0.41\t0\t0.05\n"
        "1\t300\t350\tTP53\t0.00\t102\t0.42\t0\t0.05\n"
        "1\t400\t450\tTP53\t0.00\t103\t0.43\t0\t0.05\n"
        "1\t500\t550\tTP53\t0.00\t104\t0.44\t0\t0.05\n"
        "1\t600\t650\tTP53\t0.00\t105\t0.45\t0\t0.05\n"
        "1\t700\t750\tTP53\t0.00\t106\t0.46\t0\t0.05\n"
        "1\t800\t850\tTP53\t0.00\t107\t0.47\t0\t0.05\n"
        "2\t100\t150\tMYC\t0.00\t130\t0.48\t0\t0.05\n"
        "2\t200\t250\tMYC\t0.00\t135\t0.49\t0\t0.05\n"
    )

    return pon_path


def write_cns_input(tmp_path: Path) -> Path:
    """
    Create a realistic called CNVkit .cns file.
    """
    cns_path = tmp_path / "sample.cns"
    cns_path.write_text(
        "chromosome\tstart\tend\tgene\tlog2\tbaf\tcn\tcn1\tcn2\tdepth\tprobes\tweight\n"
        "1\t100\t850\tTP53\t0.39\t0.45\t3\t1\t2\t104\t8\t10.5\n"
        "2\t100\t250\tMYC\t0.38\t0.50\t3\t1\t2\t132\t2\t3.2\n"
    )
    return cns_path


def write_raw_cns(tmp_path: Path) -> Path:
    """
    Create a realistic raw/init CNVkit .cns file.
    """
    cns_init_path = tmp_path / "sample_init.cns"
    cns_init_path.write_text(
        "chromosome\tstart\tend\tgene\tlog2\tprobes\tweight\n"
        "1\t100\t850\tTP53\t0.34\t8\t10.5\n"
        "2\t100\t250\tMYC\t0.33\t2\t3.2\n"
    )
    return cns_init_path


def write_cytoband(tmp_path: Path):
    # ------------------------------------------------------------------
    # Cytoband file (UCSC-like, no header)
    # chr, chromStart, chromEnd, name, gieStain
    # ------------------------------------------------------------------
    cytoband_path = tmp_path / "cytoband.tsv"
    cytoband_path.write_text(
        "chr1\t0\t300000000\tp36.33\tgneg\n" "chr2\t0\t300000000\tp25.3\tgneg\n"
    )
    return cytoband_path


def write_cancer_genes(tmp_path: Path):
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
    return cancer_genes_path


def write_purity(tmp_path: Path):
    # ------------------------------------------------------------------
    # Minimal purity summary CSV
    # ------------------------------------------------------------------
    purity_csv_path = tmp_path / "purity.csv"
    purity_csv_path.write_text(
        "Sampleid,Purity,Ploidy,Sex,Contamination,Flagged,Failed,Comment\n"
        "TEST_CASE_PON,0.45,2.1,male,0.01,FALSE,FALSE,ok\n"
    )
    return purity_csv_path


def write_purecn_loh_regions(tmp_path: Path) -> Path:
    """
    Create a minimal but realistic PureCN LOH regions CSV.

    Includes one real LOH segment on chr2 and a few neutral segments.
    """
    loh_path = tmp_path / "purecn_loh_regions.csv"

    loh_path.write_text(
        "Sampleid,chr,start,end,arm,C,M,type,seg.mean,num.mark,num.snps,M.flagged,maf.expected,maf.observed\n"
        "TUMOR,1,1479333,121535434,p,2,1,,0.073,2160,107,FALSE,0.5,0.483\n"
        "TUMOR,1,124535434,245029152,q,2,1,,0.073,2160,107,FALSE,0.5,0.4833\n"
        "TUMOR,2,1033890,88875337,p,2,1,,-0.15019,897,14,FALSE,0.5,0.423\n"
        "TUMOR,2,88875338,91862363,p,1.5943,0,LOH,-0.2583,15,1,TRUE,0.4060,0.4870\n"
        "TUMOR,2,95326672,237264282,q,2,1,,0.073445,1320,68,FALSE,0.5,0.484\n"
    )

    return loh_path


def test_generate_report_no_pon(tmp_path, monkeypatch):
    """
    End-to-end smoke test:
    - writes tiny input files
    - runs CLI
    - verifies report and chromosome plot outputs are created

    No PON input, and missing PureCN LOH file

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

    cnr_path = write_cnr_input(tmp_path)
    cns_path = write_cns_input(tmp_path)
    cns_init_path = write_raw_cns(tmp_path)
    cytoband_path = write_cytoband(tmp_path)
    cancer_genes_path = write_cancer_genes(tmp_path)
    purity_csv_path = write_purity(tmp_path)

    # ------------------------------------------------------------------
    # LOH regions from PureCN
    # Can not exist sometimes
    # ------------------------------------------------------------------
    purecn_loh_regions = "/not/a/real/path/CNV.somatic.fakecase.purecn.LOHregions.csv"

    # ------------------------------------------------------------------
    # Required by click as existing file
    # Contents do not matter because load_vcf_with_vaf is monkeypatched
    # ------------------------------------------------------------------
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text("##fileformat=VCFv4.2\n")

    output_html = tmp_path / "interactive_report.html"

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
            "--loh-regions",
            purecn_loh_regions,
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

    assert result.exit_code == 0, f"{result.output}\nException: {result.exception}"

    assert output_html.exists()
    html = output_html.read_text(encoding="utf-8")
    assert case_id in html
    assert "CNV Report" in html

    chr_plots_dir = tmp_path / f"{case_id}_chr_plots"
    assert chr_plots_dir.exists()

    # At least the chromosomes present in the data should have plots
    assert (chr_plots_dir / "cnv_chr1_segments.png").exists()
    assert (chr_plots_dir / "cnv_chr2_segments.png").exists()


def test_generate_report_with_pon(tmp_path, monkeypatch):
    """
    End-to-end smoke test with PON enabled.

    This exercises:
    - PON loading
    - gene-region creation
    - PON-based metrics and plotting
    """

    monkeypatch.setattr(
        cnv_report_plotting,
        "load_vcf_with_vaf",
        lambda vcf_path, chr_order: pd.DataFrame(columns=["CHROM", "POS", "VAF"]),
    )

    case_id = "TEST_CASE_PON"

    cnr_path = write_cnr_input(tmp_path)
    pon_path = write_matching_pon_input(tmp_path)
    cns_path = write_cns_input(tmp_path)
    cns_init_path = write_raw_cns(tmp_path)
    cytoband_path = write_cytoband(tmp_path)
    cancer_genes_path = write_cancer_genes(tmp_path)
    purity_csv_path = write_purity(tmp_path)

    # ------------------------------------------------------------------
    # Optional LOH file path that does not exist
    # ------------------------------------------------------------------
    purecn_loh_regions = "/not/a/real/path/CNV.somatic.fakecase.purecn.LOHregions.csv"

    # ------------------------------------------------------------------
    # VCF path exists, parsing is monkeypatched away
    # ------------------------------------------------------------------
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text("##fileformat=VCFv4.2\n")

    output_html = tmp_path / "report_with_pon.html"

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
            "--pon",
            str(pon_path),
            "--loh-regions",
            purecn_loh_regions,
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

    assert result.exit_code == 0, f"{result.output}\nException: {result.exception}"

    assert output_html.exists()
    html = output_html.read_text(encoding="utf-8")

    assert case_id in html
    assert "CNV Report" in html
    assert "Panel of Normal" in html

    chr_plots_dir = tmp_path / f"{case_id}_chr_plots"
    assert chr_plots_dir.exists()
    assert (chr_plots_dir / "cnv_chr1_segments.png").exists()
    assert (chr_plots_dir / "cnv_chr2_segments.png").exists()


def test_generate_report(tmp_path, monkeypatch):
    """
    End-to-end smoke test with everything enabled

    This exercises:
    - PON loading
    - gene-region creation
    - PON-based metrics and plotting
    """

    monkeypatch.setattr(
        cnv_report_plotting,
        "load_vcf_with_vaf",
        lambda vcf_path, chr_order: pd.DataFrame(columns=["CHROM", "POS", "VAF"]),
    )

    case_id = "TEST_CASE_PON"

    cnr_path = write_cnr_input(tmp_path)
    pon_path = write_matching_pon_input(tmp_path)
    cns_path = write_cns_input(tmp_path)
    cns_init_path = write_raw_cns(tmp_path)
    cytoband_path = write_cytoband(tmp_path)
    cancer_genes_path = write_cancer_genes(tmp_path)
    purity_csv_path = write_purity(tmp_path)

    # ------------------------------------------------------------------
    # Optional LOH file path EXISTS
    # ------------------------------------------------------------------
    purecn_loh_regions = write_purecn_loh_regions(tmp_path)

    # ------------------------------------------------------------------
    # VCF path exists, parsing is monkeypatched away
    # ------------------------------------------------------------------
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text("##fileformat=VCFv4.2\n")

    output_html = tmp_path / "report_with_pon_and_loh.html"

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
            "--pon",
            str(pon_path),
            "--loh-regions",
            str(purecn_loh_regions),
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

    assert result.exit_code == 0, f"{result.output}\nException: {result.exception}"

    assert output_html.exists()
    html = output_html.read_text(encoding="utf-8")

    assert case_id in html
    assert "CNV Report" in html
    assert "Panel of Normal" in html

    chr_plots_dir = tmp_path / f"{case_id}_chr_plots"
    assert chr_plots_dir.exists()
    assert (chr_plots_dir / "cnv_chr1_segments.png").exists()
    assert (chr_plots_dir / "cnv_chr2_segments.png").exists()


def test_generate_report_exome(tmp_path, monkeypatch):
    """
    End-to-end smoke test with everything enabled for exome

    This exercises:
    - PON loading
    - gene-region creation
    - PON-based metrics and plotting
    """

    monkeypatch.setattr(
        cnv_report_plotting,
        "load_vcf_with_vaf",
        lambda vcf_path, chr_order: pd.DataFrame(columns=["CHROM", "POS", "VAF"]),
    )

    case_id = "TEST_CASE_PON"

    cnr_path = write_cnr_input(tmp_path)
    pon_path = write_matching_pon_input(tmp_path)
    cns_path = write_cns_input(tmp_path)
    cns_init_path = write_raw_cns(tmp_path)
    cytoband_path = write_cytoband(tmp_path)
    cancer_genes_path = write_cancer_genes(tmp_path)
    purity_csv_path = write_purity(tmp_path)

    # ------------------------------------------------------------------
    # Optional LOH file path EXISTS
    # ------------------------------------------------------------------
    purecn_loh_regions = write_purecn_loh_regions(tmp_path)

    # ------------------------------------------------------------------
    # VCF path exists, parsing is monkeypatched away
    # ------------------------------------------------------------------
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text("##fileformat=VCFv4.2\n")

    output_html = tmp_path / "report_with_pon_and_loh.html"

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
            "--pon",
            str(pon_path),
            "--loh-regions",
            str(purecn_loh_regions),
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
            "--is-exome",
        ],
    )

    assert result.exit_code == 0, f"{result.output}\nException: {result.exception}"

    assert output_html.exists()
    html = output_html.read_text(encoding="utf-8")

    assert case_id in html
    assert "CNV Report" in html
    assert "Panel of Normal" in html

    chr_plots_dir = tmp_path / f"{case_id}_chr_plots"
    assert chr_plots_dir.exists()
    assert (chr_plots_dir / "cnv_chr1_segments.png").exists()
    assert (chr_plots_dir / "cnv_chr2_segments.png").exists()
