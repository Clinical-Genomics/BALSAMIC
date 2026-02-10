from __future__ import annotations
import sys
import base64
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Iterable, Any
from jinja2 import Environment, FileSystemLoader, select_autoescape
import click
import numpy as np
import pandas as pd

from cnv_report_utils import (
    build_gene_segment_table,
    build_gene_chunk_table,
    load_cancer_gene_set,
    compute_summary_metrics,
    read_purecn_summary,
    pdf_first_page_to_png,
)
from cnv_report_plotting import plot_chromosomes

from BALSAMIC.constants.analysis import Gender


CURATED_CANCER_GENES: set[str] = {
    "TP53",  # <--- requested by cust087
    "DLEU1",  # <--- requested by cust087
    "DLEU2",
    "RB1",  # <--- requested by cust087
    "KMT2D",
    "KMT2A",
    "ATM",  # <--- requested by cust087
}

# =============================================================================
# Public API
# =============================================================================


def render_cnv_report_html(
    *,
    df_gene: pd.DataFrame,
    out_html: str | Path,
    df_chunk: pd.DataFrame | None = None,
    df_purecn_summary: pd.DataFrame | None = None,
    df_qc_summary: pd.DataFrame | None = None,
    scatter_png: str | Path | None = None,
    diagram_png: str | Path | None = None,
    chr_plots_dir: str | Path | None = None,
    title: str = "CNV Report",
) -> None:
    """
    Render a standalone CNV QC HTML report using a Jinja2 template + embedded CSS/JS assets.

    Includes:
      - Optional genome-wide CNVkit scatter/diagram PNGs (embedded as base64 data URIs)
      - Optional per-chromosome plot cards (cnv_chr*segments.png) (NO gene-level plots)
      - Gene-level and optional chunk-level DataFrame tables (Pandas to_html)
      - Client-side filtering and plot modal viewer

    Notes:
      - Column index maps are injected as JSON for JS filtering logic.
      - Plot cards are structured data; template loops and renders the HTML.
    """
    out_path = Path(out_html)

    # ---- base64 URIs for genome-wide plots
    scatter_data_uri = _png_to_data_uri(scatter_png)
    diagram_data_uri = _png_to_data_uri(diagram_png)

    # ---- CNV/LOH sets used to badge chromosome plots
    cnv_chr_set, _cnv_gene_set_unused, chr_col = _compute_cnv_sets(df_gene)

    # ---- Collect chromosome plot cards only
    plot_groups = _collect_chr_plot_groups(
        chr_plots_dir=chr_plots_dir,
        chr_col=chr_col,
        cnv_chr_set=cnv_chr_set,
    )

    # ---- Summary tables (optional)
    purecn_summary_html = ""
    if df_purecn_summary is not None and not df_purecn_summary.empty:
        purecn_summary_html = df_purecn_summary.to_html(
            index=False, border=0, classes="dataframe", table_id="purecn-summary-table"
        )

    qc_summary_html = ""
    if df_qc_summary is not None and not df_qc_summary.empty:
        qc_display = df_qc_summary.copy()
        for col in qc_display.columns:
            if np.issubdtype(qc_display[col].dtype, np.floating):
                qc_display[col] = qc_display[col].round(3)
        qc_summary_html = qc_display.to_html(
            index=False, border=0, classes="dataframe", table_id="qc-summary-table"
        )

    # ---- Main gene table
    gene_table_html = df_gene.to_html(
        index=False, border=0, classes="dataframe", table_id="report-table"
    )

    # ---- Optional chunk table
    has_chunk_table = df_chunk is not None and not df_chunk.empty
    chunk_table_html = ""
    if has_chunk_table:
        chunk_table_html = df_chunk.to_html(
            index=False, border=0, classes="dataframe", table_id="chunk-table"
        )

    # ---- Column index JSON maps for JS filtering
    col_idx_gene_json = json.dumps(
        {name: idx for idx, name in enumerate(df_gene.columns)}
    )
    col_idx_chunk_json = (
        json.dumps({name: idx for idx, name in enumerate(df_chunk.columns)})
        if has_chunk_table and df_chunk is not None
        else "{}"
    )

    # ---- Load template + assets
    base_dir = Path(__file__).resolve().parent
    template_dir = base_dir / "templates"
    assets_dir = base_dir / "assets"

    css_text = _read_text(assets_dir / "cnv_report.css")
    js_text = _read_text(assets_dir / "cnv_report.js")

    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("cnv_report.html.j2")

    html = template.render(
        title=title,
        css_text=css_text,
        js_text=js_text,
        # Tables (already HTML)
        purecn_summary_html=purecn_summary_html,
        qc_summary_html=qc_summary_html,
        gene_table_html=gene_table_html,
        has_chunk_table=has_chunk_table,
        chunk_table_html=chunk_table_html,
        # Genome-wide plots
        scatter_data_uri=scatter_data_uri,
        diagram_data_uri=diagram_data_uri,
        # Chromosome plot cards only
        plot_groups=plot_groups,
        # JSON for JS
        col_idx_gene_json=col_idx_gene_json,
        col_idx_chunk_json=col_idx_chunk_json,
    )

    out_path.write_text(html, encoding="utf-8")


# =============================================================================
# Helpers
# =============================================================================


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.is_file() else ""


def _png_to_data_uri(png_path: str | Path | None) -> str | None:
    if png_path is None:
        return None
    p = Path(png_path)
    if not p.is_file():
        return None
    b64 = base64.b64encode(p.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def _chr_sort_key_from_stem(stem: str) -> tuple[int, int]:
    label = stem
    if "chr" in stem:
        label = stem.split("chr", 1)[1]
    if "_" in label:
        label = label.split("_", 1)[0]
    if label.isdigit():
        return (0, int(label))
    special = {"X": 23, "Y": 24}
    return (1, special.get(label.upper(), 25))


def _extract_chr_label_from_stem(stem: str) -> str:
    chr_label = stem
    if "chr" in stem:
        chr_label = stem.split("chr", 1)[1]
    if "_" in chr_label:
        chr_label = chr_label.split("_", 1)[0]
    return chr_label


def _as_upper_str(x: Any) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def _is_true_str(x: Any) -> bool:
    return _as_upper_str(x) == "TRUE"


def _is_amp_del(x: Any) -> bool:
    return _as_upper_str(x) in {"DELETION", "AMPLIFICATION"}


def _compute_cnv_sets(df: pd.DataFrame) -> tuple[set[str], set[str], str | None]:
    """
    Return:
      - cnv_chr_set: chromosomes that have LOH/CNV calls
      - cnv_gene_set: genes that have LOH/CNV calls (computed but not used for plots here)
      - chr_col: which column name to use ("chr" or "chromosome") if present
    """
    df_chr = df.copy()

    if "chr" in df_chr.columns:
        chr_col: str | None = "chr"
    elif "chromosome" in df_chr.columns:
        chr_col = "chromosome"
    else:
        chr_col = None

    cnv_chr_set: set[str] = set()
    cnv_gene_set: set[str] = set()

    if chr_col is None:
        return cnv_chr_set, cnv_gene_set, chr_col

    if "gene.symbol" in df_chr.columns:
        df_chr = df_chr[~df_chr["gene.symbol"].isin(["Antitarget", "-"])]

    cnv_mask = pd.Series(False, index=df_chr.index)

    if "loh_flag" in df_chr.columns:
        cnv_mask |= df_chr["loh_flag"].apply(_is_true_str)
    if "cnvkit_cnv_call" in df_chr.columns:
        cnv_mask |= df_chr["cnvkit_cnv_call"].apply(_is_amp_del)
    if "purecn_cnv_call" in df_chr.columns:
        cnv_mask |= df_chr["purecn_cnv_call"].apply(_is_amp_del)

    cnv_chr_set = set(df_chr.loc[cnv_mask, chr_col].astype(str).tolist())

    if "gene.symbol" in df_chr.columns:
        cnv_gene_set = set(df_chr.loc[cnv_mask, "gene.symbol"].astype(str).tolist())

    return cnv_chr_set, cnv_gene_set, chr_col


@dataclass(frozen=True)
class PlotCard:
    title: str
    data_uri: str
    stem: str
    has_cnv: bool


def _collect_chr_plot_groups(
    *,
    chr_plots_dir: str | Path | None,
    chr_col: str | None,
    cnv_chr_set: set[str],
) -> dict[str, list[dict[str, Any]]]:
    """
    Return structured chromosome plot card groups for templating:

      {
        "chr_with_cnv": [ {title, data_uri, stem, has_cnv}, ... ],
        "chr_no_cnv":   [ ... ],
      }

    Only chromosome-level plots are scanned. Gene-zoom plots are ignored.
    """
    groups: dict[str, list[PlotCard]] = {
        "chr_with_cnv": [],
        "chr_no_cnv": [],
    }

    if chr_plots_dir is None:
        return {k: [] for k in groups}

    qc_dir = Path(chr_plots_dir)
    if not qc_dir.is_dir():
        return {k: [] for k in groups}

    png_files = sorted(
        list(qc_dir.glob("cnv_chr*segments.png")),
        key=lambda p: _chr_sort_key_from_stem(p.stem),
    )

    for png_file in png_files:
        stem = png_file.stem

        # Ignore gene-zoom plots (anything with "_genes_" in the stem)
        if "_genes_" in stem:
            continue

        data_uri = _png_to_data_uri(png_file)
        if not data_uri:
            continue

        chr_label = _extract_chr_label_from_stem(stem)
        has_cnv = bool(chr_col is not None and chr_label in cnv_chr_set)
        title = f"Chr {chr_label}"

        card = PlotCard(title=title, data_uri=data_uri, stem=stem, has_cnv=has_cnv)
        (groups["chr_with_cnv"] if has_cnv else groups["chr_no_cnv"]).append(card)

    return {k: [c.__dict__ for c in v] for k, v in groups.items()}


@click.command()
@click.option(
    "--loh-genes",
    type=click.Path(exists=True),
    required=True,  # PureCN genes CSV (can be empty but must exist)
    help="PureCN LOH genes CSV (…_genes.csv).",
)
@click.option(
    "--cnr",
    type=click.Path(exists=True),
    required=True,
    help="CNVkit tumor .cnr file.",
)
@click.option(
    "--cns",
    type=click.Path(exists=True),
    required=True,
    help="CNVkit tumor .cns file.",
)
@click.option(
    "--cns-init",
    type=click.Path(exists=True),
    required=True,
    help="CNVkit tumor .cns initial file. (not adjusted for purity and ploidy)",
)
@click.option(
    "--pon",
    type=click.Path(exists=True),
    required=False,  # <- NOW OPTIONAL
    help="CNVkit PON .cnn file (optional; if absent, plots and table are built without PON spread).",
)
@click.option(
    "--vcf",
    type=click.Path(exists=True),
    required=True,
    help="VCF with germline variants for BAF.",
)
@click.option(
    "--refgene",
    type=click.Path(exists=True),
    required=True,
    help="refgene.flat file.",
)
@click.option(
    "--cytoband",
    type=click.Path(exists=True),
    required=True,
    help="cytoBand file for genome build.",
)
@click.option(
    "--case-id",
    type=str,
    required=True,
    help="Case ID for labels / outputs.",
)
@click.option(
    "--cust-case-id",
    type=str,
    required=False,
    help="Cust Case ID for labels / outputs.",
)
@click.option(
    "--cnvkit-scatter",
    type=click.Path(exists=True),
    required=False,
    help="cnvkit scatter PDF file.",
)
@click.option(
    "--cnvkit-diagram",
    type=click.Path(exists=True),
    required=False,
    help="cnvkit diagram PDF file.",
)
@click.option(
    "--out-prefix",
    type=click.Path(),
    required=True,
    help="Output prefix (e.g. /path/to/sample_cnv_qc.html).",
)
@click.option(
    "--purity-csv",
    type=click.Path(exists=True),
    required=False,
    help="PureCN sample summary CSV (Purity, Ploidy, etc., optional).",
)
@click.option(
    "--cancer-genes",
    type=click.Path(exists=True),
    required=False,
    help="Cancer gene list TSV (OncoKB/CGC aggregate, optional).",
)
@click.option(
    "--is-exome",
    is_flag=True,
    default=False,
    help="Treat assay as exome (use more aggressive non-driver bin compression, etc.).",
)
@click.option(
    "--sex",
    type=click.Choice([Gender.FEMALE, Gender.MALE]),
    default=None,
)
def main(
    loh_genes: str,
    cnr: str,
    cns: str,
    cns_init: str,
    pon: Optional[str],
    vcf: str,
    refgene: str,
    cytoband: str,
    case_id: str,
    cust_case_id: Optional[str],
    cnvkit_scatter: Optional[str],
    cnvkit_diagram: Optional[str],
    out_prefix: str,
    purity_csv: Optional[str],
    cancer_genes: Optional[str],
    is_exome: bool,
    sex: Optional[Gender],
):
    """
    Build CNV QC plots + HTML report from CNVkit outputs, optionally
    annotated with PureCN segments/genes and PON spread.
    """

    out_prefix = Path(out_prefix)
    outdir = out_prefix.parent
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------------
    # Convert CNVkit PDFs to PNG
    # ---------------
    scatter_png_path = outdir / f"cnvkit_scatter_{case_id}.png"
    diagram_png_path = outdir / f"cnvkit_diagram_{case_id}.png"

    if cnvkit_scatter and Path(cnvkit_scatter).is_file():
        pdf_first_page_to_png(cnvkit_scatter, scatter_png_path)
    if cnvkit_diagram and Path(cnvkit_diagram).is_file():
        pdf_first_page_to_png(cnvkit_diagram, diagram_png_path)

    # ----------------------------
    # Load cancer gene list (optional)
    # ----------------------------
    if cancer_genes:
        cancer_gene_set = load_cancer_gene_set(
            cancer_genes,
            min_occurrence=3,
            only_annotated=True,
        )
    else:
        cancer_gene_set = set()
    cancer_gene_set |= CURATED_CANCER_GENES

    # ----------------------------
    # Create per-gene CNV table
    # ----------------------------
    # pon can be None (no PON file); build_gene_segment_table handles that
    genes_df = build_gene_segment_table(
        cnr_path=cnr,
        cns_path=cns,
        cns_init_path=cns_init,
        cancer_genes=cancer_gene_set,
        refgene_path=refgene,
        loh_path=loh_genes,
        cytoband_path=cytoband,
        sex=sex,
        pon_path=pon,
    )

    # ----------------------------
    # Create per chunk table
    # ----------------------------
    chunks_df = build_gene_chunk_table(
        cnr_path=cnr,
        cns_path=cns,
        cns_init_path=cns_init,
        cancer_genes=cancer_gene_set,
        refgene_path=refgene,
        loh_path=loh_genes,
        cytoband_path=cytoband,
        sex=sex,
        pon_path=pon,
    )

    # --- 1) PureCN summary (from purity_csv) ---
    purecn_summary_df = read_purecn_summary(purity_csv)

    # --- 2) Extra QC / CNV metrics (DLR, PON spread, chunk stats) ---
    qc_summary_df = compute_summary_metrics(
        cnr_path=cnr,
        cnn_path=pon,  # CNVkit .cnn PON, or None
    )

    # ----------------------------
    # Generate per-chromosome PNG plots in outdir
    # ----------------------------
    chr_plots_dir = outdir / f"{case_id}_chr_plots"
    chr_plots_dir.mkdir(exist_ok=True, parents=True)
    plot_case_id = cust_case_id if cust_case_id else case_id

    # Different compression for exome vs panel
    neutral_target_factor = 0.1 if is_exome else 0.6
    highlight_only_cancer = True if is_exome else False

    # pon may be None; plot_chromosomes is written to handle that
    plot_chromosomes(
        cnr_path=cnr,
        vcf_path=vcf,
        gdf=genes_df,
        gchunk=chunks_df,
        outdir=chr_plots_dir,
        case_id=plot_case_id,
        pon_path=pon,
        include_y=True,
        neutral_target_factor=neutral_target_factor,
        highlight_only_cancer=highlight_only_cancer,
    )

    if is_exome:
        if "is_cancer_gene" in genes_df.columns:
            genes_df = genes_df[genes_df["is_cancer_gene"].fillna(False).astype(bool)]
        if (
            chunks_df is not None
            and not chunks_df.empty
            and "is_cancer_gene" in chunks_df.columns
        ):
            chunks_df = chunks_df[
                chunks_df["is_cancer_gene"].fillna(False).astype(bool)
            ]

    # ----------------------------
    # HTML report
    # ----------------------------
    out_html_path = str(out_prefix)

    render_cnv_report_html(
        df_gene=genes_df,
        df_chunk=chunks_df,  # optional
        df_purecn_summary=purecn_summary_df,  # optional
        df_qc_summary=qc_summary_df,  # optional
        scatter_png=scatter_png_path,  # optional
        diagram_png=diagram_png_path,  # optional
        chr_plots_dir=chr_plots_dir,
        out_html=out_html_path,
        title=f"CNV Report – {case_id}",
    )

    click.echo(f"[CNV QC] Finished report for {case_id}: {out_html_path}")


if __name__ == "__main__":
    sys.exit(main())
