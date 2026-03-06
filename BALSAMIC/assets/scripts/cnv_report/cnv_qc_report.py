from __future__ import annotations
import sys
import base64
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Any
from jinja2 import Environment, FileSystemLoader, select_autoescape
import click
import numpy as np
import pandas as pd

from cnv_summary_metrics import read_purecn_summary, compute_summary_metrics
from cnv_report_utils import pdf_first_page_to_png
from cnv_io import (
    load_cancer_gene_set,
    load_cnr_bins,
    load_pon_bins,
    load_refgene_exons,
    load_purecn_segments,
    load_cnvkit_segments_with_raw,
    load_cytobands,
)

from cnv_tables import (
    build_gene_chunk_table,
    build_segment_table,
)
from cnv_report_plotting import plot_chromosomes
from cnv_constants import (
    GENE_TABLE_SPEC,
    SEGMENT_TABLE_SPEC,
    TableSpec,
    PURECN_WARNING_TEXT,
)

from BALSAMIC.constants.analysis import Gender, AnalysisType


def render_cnv_report_html(
    *,
    df_segments: pd.DataFrame,
    out_html: str | Path,
    df_chunk: pd.DataFrame | None = None,
    df_purecn_summary: pd.DataFrame | None = None,
    df_qc_summary: pd.DataFrame | None = None,
    scatter_png: str | Path | None = None,
    diagram_png: str | Path | None = None,
    chr_plots_dir: str | Path | None = None,
    title: str = "CNV Report",
    normalisation_method: str,
) -> None:
    """
    Render a standalone CNV QC HTML report using a Jinja2 template + embedded CSS/JS assets.
    """
    out_path = Path(out_html)

    # ---- base64 URIs for genome-wide plots
    scatter_data_uri = _png_to_data_uri(scatter_png)
    diagram_data_uri = _png_to_data_uri(diagram_png)

    # ---- CNV/LOH sets used to badge chromosome plots
    cnv_chr_set = _compute_cnv_sets(df_segments)

    # ---- Collect chromosome plot cards only
    plot_groups = _collect_chr_plot_groups(
        chr_plots_dir=chr_plots_dir,
        cnv_chr_set=cnv_chr_set,
    )

    # ---- PureCN summary table + warning state
    purecn_summary_html = ""
    purecn_failed = False
    purecn_failed_warning_text = ""

    if df_purecn_summary is not None and not df_purecn_summary.empty:
        purecn_display = df_purecn_summary.copy()

        if "Comment" in purecn_display.columns:
            failed_mask = (
                purecn_display["Comment"]
                .astype("string")
                .str.contains("FAILED PURITY ESTIMATION", case=False, na=False)
            )

        purecn_failed = bool(failed_mask.any())

        if purecn_failed:
            if "Purity" in purecn_display.columns:
                purecn_display.loc[failed_mask, "Purity"] = purecn_display.loc[
                    failed_mask, "Purity"
                ].map(lambda x: f"{x} (default fallback 20% purity)")

            if "Ploidy" in purecn_display.columns:
                purecn_display.loc[failed_mask, "Ploidy"] = purecn_display.loc[
                    failed_mask, "Ploidy"
                ].map(lambda x: f"{x} (default fallback 2 ploidy)")

            purecn_failed_warning_text = PURECN_WARNING_TEXT

        purecn_summary_html = purecn_display.to_html(
            index=False,
            border=0,
            classes="dataframe",
            table_id="purecn-summary-table",
        )

    # ---- QC summary table
    qc_summary_html = ""
    if df_qc_summary is not None and not df_qc_summary.empty:
        qc_display = df_qc_summary.copy()
        for col in qc_display.columns:
            if np.issubdtype(qc_display[col].dtype, np.floating):
                qc_display[col] = qc_display[col].round(3)
        qc_summary_html = qc_display.to_html(
            index=False, border=0, classes="dataframe", table_id="qc-summary-table"
        )

    # ---- Add column glossary
    segment_column_glossary_html = build_column_glossary_html(
        table_df=df_segments,
        spec=SEGMENT_TABLE_SPEC,
    )

    chunk_column_glossary_html = build_column_glossary_html(
        table_df=df_chunk,
        spec=GENE_TABLE_SPEC,
    )

    # ---- Main segments table
    segments_table_html = df_for_html(df_segments).to_html(
        index=False, border=0, classes="dataframe", table_id="report-table", na_rep=""
    )

    # ---- Optional chunk table
    has_chunk_table = df_chunk is not None and not df_chunk.empty
    chunk_table_html = ""
    if has_chunk_table:
        chunk_table_html = df_for_html(df_chunk).to_html(
            index=False,
            border=0,
            classes="dataframe",
            table_id="chunk-table",
            na_rep="",
        )

    # ---- Column index JSON maps for JS filtering
    col_idx_segments_json = json.dumps(
        {name: idx for idx, name in enumerate(df_segments.columns)}
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
        normalisation_method=normalisation_method,
        css_text=css_text,
        js_text=js_text,
        segment_column_glossary_html=segment_column_glossary_html,
        chunk_column_glossary_html=chunk_column_glossary_html,
        purecn_summary_html=purecn_summary_html,
        qc_summary_html=qc_summary_html,
        segments_table_html=segments_table_html,
        has_chunk_table=has_chunk_table,
        chunk_table_html=chunk_table_html,
        scatter_data_uri=scatter_data_uri,
        diagram_data_uri=diagram_data_uri,
        plot_groups=plot_groups,
        col_idx_segments_json=col_idx_segments_json,
        col_idx_chunk_json=col_idx_chunk_json,
        purecn_failed=purecn_failed,
        purecn_warning_text=purecn_failed_warning_text,
    )

    out_path.write_text(html, encoding="utf-8")


# =============================================================================
# Helpers
# =============================================================================


def df_for_html(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a display-only copy where all missing values render as empty cells.
    Does not modify the original DataFrame.
    """
    out = df.copy().astype("object")

    # Replace true missing values
    out = out.where(pd.notna(out), "")

    # Remove literal "<NA>" strings if they exist
    out[out == "<NA>"] = ""

    return out


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


def _compute_cnv_sets(df: pd.DataFrame) -> set[str]:
    """
    Return:
      - cnv_chr_set: chromosomes that have LOH/CNV calls
      - cnv_gene_set: genes that have LOH/CNV calls (computed but not used for plots here)
      - chr_col: which column name to use ("chr" or "chromosome") if present
    """
    df_chr = df.copy()

    cnv_mask = pd.Series(False, index=df_chr.index)

    if "purecn_type" in df_chr.columns:
        cnv_mask |= df_chr["purecn_type"].notna()
    if "cnv_call" in df_chr.columns:
        cnv_mask |= df_chr["cnv_call"].apply(_is_amp_del)

    cnv_chr_set = set(df_chr.loc[cnv_mask, "chr"].astype(str).tolist())

    return cnv_chr_set


def _ordered_glossary_columns(
    *,
    df_chunk: pd.DataFrame,
    spec: TableSpec,
) -> list[str]:
    """Return glossary column names in TableSpec order, then extras."""
    # columns we actually need to describe (present in either table)
    present: set[str] = set(df_chunk.columns)

    # 1) spec order (only those that exist)
    ordered = [c for c in spec.column_order if c in present]

    # 2) extras (present but not in spec), append deterministically
    extras = sorted(present - set(spec.column_order))
    return ordered + extras


def build_column_glossary_html(
    *,
    table_df: pd.DataFrame | None,
    spec: TableSpec,
) -> str:
    """
    Build an HTML table describing columns in the same order as TableSpec.column_order,
    then any extra columns found in df_gene/df_chunk appended at the end.
    """
    cols = _ordered_glossary_columns(df_chunk=table_df, spec=spec)

    rows = []
    for c in cols:
        desc = spec.descriptions.get(c, "")
        rows.append({"column": c, "description": desc})

    glossary_df = pd.DataFrame(rows, columns=["column", "description"])

    # Use to_html for consistent styling; pick classes you already use
    return glossary_df.to_html(
        index=False,
        border=0,
        classes="dataframe glossary-table",
        table_id="column-glossary-table",
        escape=False,
    )


@dataclass(frozen=True)
class PlotCard:
    title: str
    data_uri: str
    stem: str
    has_cnv: bool


def _collect_chr_plot_groups(
    *,
    chr_plots_dir: str | Path | None,
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
        has_cnv = bool(chr_label in cnv_chr_set)
        title = f"Chr {chr_label}"

        card = PlotCard(title=title, data_uri=data_uri, stem=stem, has_cnv=has_cnv)
        (groups["chr_with_cnv"] if has_cnv else groups["chr_no_cnv"]).append(card)

    return {k: [c.__dict__ for c in v] for k, v in groups.items()}


@click.command()
@click.option(
    "--loh-regions",
    type=click.Path(exists=False),
    required=False,
    help="PureCN LOH regions CSV",
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
    required=False,
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
    "--output-file",
    type=click.Path(),
    required=True,
    help="Output (e.g. /path/to/sample_cnv_qc.html).",
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
    required=True,
    help="Sample sex",
)
@click.option(
    "--analysis-type",
    type=str,
    required=True,
    help="Paired / Single",
)
def main(
    loh_regions: str,
    cnr: str,
    cns: str,
    cns_init: str,
    pon: str | None,
    vcf: str,
    refgene: str,
    cytoband: str,
    case_id: str,
    cnvkit_scatter: str | None,
    cnvkit_diagram: str | None,
    output_file: str,
    purity_csv: str | None,
    cancer_genes: str | None,
    is_exome: bool,
    sex: Gender,
    analysis_type: str,
):
    """
    Build CNV QC plots + HTML report from CNVkit outputs, optionally
    annotated with PureCN segments/genes and PON spread.
    """

    if pon:
        normalisation_method = f"Panel of Normal: {Path(pon).stem}"
    elif analysis_type == AnalysisType.PAIRED:
        normalisation_method = "Matched Normal"
    else:
        normalisation_method = "Flat reference (tumor only)"

    out_prefix = Path(output_file)
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

    cancer_gene_set = load_cancer_gene_set(
        cancer_genes,
        min_occurrence=1,
    )

    # ----------------------------
    # Read files into dataframes
    # ----------------------------
    cnr_df = load_cnr_bins(cnr)
    cns_df = load_cnvkit_segments_with_raw(cns, cns_init)
    loh_regions_df = None
    if Path(loh_regions).is_file():
        loh_regions_df = load_purecn_segments(loh_regions)
    pon_df = None
    if pon:
        pon_df = load_pon_bins(pon)
    cytoband_df = load_cytobands(cytoband)

    # exon_map: Dict[Tuple[str, str], dict] = load_refgene_exons(refgene)

    # ----------------------------
    # Create segment table
    # ----------------------------
    segments_df = build_segment_table(
        cnr_df=cnr_df,
        cns_df=cns_df,
        cytoband_df=cytoband_df,
        sex=sex,
        cancer_genes=cancer_gene_set,
        is_exome=is_exome,
        loh_regions_df=loh_regions_df,
    )

    # ----------------------------
    # Create per chunk table
    # ----------------------------
    chunks_df = build_gene_chunk_table(
        cnr_df=cnr_df,
        cns_df=cns_df,
        sex=sex,
        pon_df=pon_df,
        cancer_genes=cancer_gene_set,
        loh_regions_df=loh_regions_df,
    )

    # --- 1) PureCN summary (from purity_csv) ---
    purecn_summary_df = read_purecn_summary(purity_csv)

    # --- 2) Extra QC / CNV metrics (DLR, PON spread, chunk stats) ---
    qc_summary_df = compute_summary_metrics(
        cnr_df=cnr_df,
        pon_df=pon_df,
    )

    # ----------------------------
    # Generate per-chromosome PNG plots in outdir
    # ----------------------------
    chr_plots_dir = outdir / f"{case_id}_chr_plots"
    chr_plots_dir.mkdir(exist_ok=True, parents=True)

    # Different compression for exome vs panel
    neutral_target_factor = 0.07 if is_exome else 0.4
    backbone_factor = 0.07 if is_exome else 0.4
    highlight_only_cancer = True if is_exome else False

    # pon may be None; plot_chromosomes is written to handle that
    plot_chromosomes(
        cnr_df=cnr_df,
        vcf_path=vcf,
        segments_df=segments_df,
        gchunk=chunks_df,
        outdir=chr_plots_dir,
        case_id=case_id,
        pon_df=pon_df,
        backbone_factor=backbone_factor,
        neutral_target_factor=neutral_target_factor,
        highlight_only_cancer=highlight_only_cancer,
    )

    if is_exome:
        chunks_df = chunks_df[chunks_df["is_cancer_gene"].fillna(False).astype(bool)]

    # ----------------------------
    # HTML report
    # ----------------------------
    out_html_path = str(output_file)

    render_cnv_report_html(
        df_segments=segments_df,
        df_chunk=chunks_df,  # optional
        df_purecn_summary=purecn_summary_df,  # optional
        df_qc_summary=qc_summary_df,  # optional
        scatter_png=scatter_png_path,  # optional
        diagram_png=diagram_png_path,  # optional
        chr_plots_dir=chr_plots_dir,
        out_html=out_html_path,
        title=f"CNV Report – {case_id}",
        normalisation_method=normalisation_method,
    )

    click.echo(f"[CNV QC] Finished report for {case_id}: {out_html_path}")


if __name__ == "__main__":
    sys.exit(main())
