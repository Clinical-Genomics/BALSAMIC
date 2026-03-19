from __future__ import annotations
import sys
import base64
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from jinja2 import Environment, FileSystemLoader, select_autoescape
import click
import numpy as np
import pandas as pd
import re
from datetime import datetime


from cnv_summary_metrics import read_purecn_summary, compute_summary_metrics
from cnv_report_utils import pdf_first_page_to_png
from cnv_io import (
    load_cancer_gene_set,
    load_cnr_bins,
    load_pon_bins,
    load_purecn_segments,
    load_cnvkit_segments_with_raw,
    load_cytobands,
)

from cnv_tables import (
    build_generegion_table,
    build_segment_table,
)
from cnv_report_plotting import plot_chromosomes
from cnv_constants import (
    GENE_TABLE_SPEC,
    SEGMENT_TABLE_SPEC,
    TableSpec,
    PURECN_WARNING_TEXT,
    PANEL_PLOT_CONFIG,
    EXOME_PLOT_CONFIG,
)

from BALSAMIC.constants.analysis import Gender, AnalysisType


def render_cnv_report_html(
    *,
    df_segments: pd.DataFrame,
    out_html: str | Path,
    df_regions: pd.DataFrame,
    df_purecn_summary: pd.DataFrame,
    df_qc_summary: pd.DataFrame,
    scatter_png: str,
    diagram_png: str,
    chr_plots_dir: str,
    title: str = "CNV Report",
    normalisation_method: str,
    sample_sex: str,
) -> None:
    """
    Render a standalone CNV QC HTML report using a Jinja2 template with embedded assets.
    """
    out_path = Path(out_html)

    scatter_data_uri = _png_to_data_uri(scatter_png)
    diagram_data_uri = _png_to_data_uri(diagram_png)

    cnv_chr_set = _compute_cnv_sets(df_segments)
    plot_groups = _collect_chr_plot_groups(
        chr_plots_dir=chr_plots_dir,
        cnv_chr_set=cnv_chr_set,
    )

    (
        purecn_summary_html,
        purecn_failed,
        purecn_failed_warning_text,
    ) = _build_purecn_summary_html(df_purecn_summary)
    qc_summary_html = _build_qc_summary_html(df_qc_summary)

    (
        df_segments_display,
        segments_table_html,
        segment_column_glossary_html,
    ) = _build_display_table_html(
        df_segments,
        SEGMENT_TABLE_SPEC,
        table_id="report-table",
    )

    (
        df_regions_display,
        region_table_html,
        region_column_glossary_html,
    ) = _build_display_table_html(
        df_regions,
        GENE_TABLE_SPEC,
        table_id="region-table",
    )

    col_idx_segments_json = _column_index_json(df_segments_display)
    col_idx_regions_json = _column_index_json(df_regions_display)

    env, css_text, js_text = _load_report_template_assets()
    template = env.get_template("cnv_report.html.j2")

    html = template.render(
        title=title,
        normalisation_method=normalisation_method,
        sample_sex=sample_sex,
        creation_time=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        css_text=css_text,
        js_text=js_text,
        segment_column_glossary_html=segment_column_glossary_html,
        region_column_glossary_html=region_column_glossary_html,
        purecn_summary_html=purecn_summary_html,
        qc_summary_html=qc_summary_html,
        segments_table_html=segments_table_html,
        region_table_html=region_table_html,
        scatter_data_uri=scatter_data_uri,
        diagram_data_uri=diagram_data_uri,
        plot_groups=plot_groups,
        col_idx_segments_json=col_idx_segments_json,
        col_idx_regions_json=col_idx_regions_json,
        purecn_failed=purecn_failed,
        purecn_warning_text=purecn_failed_warning_text,
    )

    out_path.write_text(html, encoding="utf-8")


def _load_report_template_assets() -> tuple[Environment, str, str]:
    """
    Load Jinja environment and embedded CSS/JS assets for the CNV report.
    """
    base_dir = Path(__file__).resolve().parent
    template_dir = base_dir / "templates"
    assets_dir = base_dir / "assets"

    css_text = _read_text(assets_dir / "cnv_report.css")
    js_text = _read_text(assets_dir / "cnv_report.js")

    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(["html", "xml"]),
    )
    return env, css_text, js_text


def _column_index_json(df: pd.DataFrame) -> str:
    """
    Return JSON mapping of display column name to column index.
    """
    return json.dumps({name: idx for idx, name in enumerate(df.columns)})


def _build_display_table_html(
    df: pd.DataFrame,
    spec,
    *,
    table_id: str,
) -> tuple[pd.DataFrame, str, str]:
    """
    Build display DataFrame, HTML table, and glossary HTML for one report table.

    Returns
    -------
    display_df
        Renamed/formatted DataFrame used for rendering.
    table_html
        HTML representation of the display DataFrame.
    glossary_html
        HTML glossary table for the display columns.
    """
    glossary_html = build_column_glossary_html(
        table_df=df,
        spec=spec,
        table_id=f"{table_id}-column-glossary-table",
    )

    display_df = rename_for_display(df, spec)

    table_html = df_for_html(display_df).to_html(
        index=False,
        border=0,
        classes="dataframe",
        table_id=table_id,
        na_rep="",
    )

    return display_df, table_html, glossary_html


def _build_qc_summary_html(df_qc_summary: pd.DataFrame | None) -> str:
    """
    Build QC summary HTML table.
    """
    if df_qc_summary is None or df_qc_summary.empty:
        return ""

    qc_display = df_qc_summary.copy()
    for col in qc_display.columns:
        if np.issubdtype(qc_display[col].dtype, np.floating):
            qc_display[col] = qc_display[col].round(3)

    return qc_display.to_html(
        index=False,
        border=0,
        classes="dataframe",
        table_id="qc-summary-table",
    )


def _build_purecn_summary_html(
    df_purecn_summary: pd.DataFrame | None,
) -> tuple[str, bool, str]:
    """
    Build PureCN summary HTML and failure warning state.

    Returns
    -------
    purecn_summary_html
        HTML table for the PureCN summary, or empty string if unavailable.
    purecn_failed
        True if a failed purity estimation was detected.
    purecn_warning_text
        Warning text to show in the report when PureCN failed.
    """
    if df_purecn_summary is None or df_purecn_summary.empty:
        return "", False, ""

    purecn_display = df_purecn_summary.copy()

    failed_mask = pd.Series(False, index=purecn_display.index)
    if "Comment" in purecn_display.columns:
        failed_mask = (
            purecn_display["Comment"]
            .astype("string")
            .str.contains("FAILED PURITY ESTIMATION", case=False, na=False)
        )

    purecn_failed = bool(failed_mask.any())
    purecn_warning_text = PURECN_WARNING_TEXT if purecn_failed else ""

    if purecn_failed:
        if "Purity" in purecn_display.columns:
            purecn_display.loc[failed_mask, "Purity"] = purecn_display.loc[
                failed_mask, "Purity"
            ].map(lambda x: f"{x} (default fallback 20% purity)")

        if "Ploidy" in purecn_display.columns:
            purecn_display.loc[failed_mask, "Ploidy"] = purecn_display.loc[
                failed_mask, "Ploidy"
            ].map(lambda x: f"{x} (default fallback 2 ploidy)")

    purecn_summary_html = purecn_display.to_html(
        index=False,
        border=0,
        classes="dataframe",
        table_id="purecn-summary-table",
    )
    return purecn_summary_html, purecn_failed, purecn_warning_text


def rename_for_display(df: pd.DataFrame, spec: TableSpec) -> pd.DataFrame:
    """Return a display copy with spec-based column renames applied."""
    return df.rename(columns=spec.renames or {})


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
    m = re.search(r"chr([A-Za-z0-9]+)", stem)
    return m.group(1) if m else stem


def _as_upper_str(x: Any) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


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


def build_column_glossary_html(
    *,
    table_df: pd.DataFrame | None,
    spec: TableSpec,
    table_id: str,
) -> str:
    """
    Build an HTML glossary table in TableSpec order, then extras.
    Uses raw column names for lookup and display labels from spec.renames.
    """
    if table_df is None or table_df.empty:
        return ""

    present = set(table_df.columns)
    renames = spec.renames or {}
    descriptions = spec.descriptions or {}

    ordered_raw = [c for c in spec.column_order if c in present]
    extras_raw = sorted(present - set(spec.column_order))
    cols_raw = ordered_raw + extras_raw

    rows = []
    for raw_col in cols_raw:
        rows.append(
            {
                "column": renames.get(raw_col, raw_col),
                "description": descriptions.get(raw_col, ""),
            }
        )

    glossary_df = pd.DataFrame(rows, columns=["column", "description"])

    return glossary_df.to_html(
        index=False,
        border=0,
        classes="dataframe glossary-table",
        table_id=table_id,
        escape=False,
    )


@dataclass(frozen=True)
class PlotCard:
    """
    Data container for a single chromosome plot card used in the CNV report template.

    Each instance represents one PNG plot and provides the minimal information
    required by the Jinja2 template to render it in the chromosome plot grid.

    Attributes
    ----------
    title
        Human-readable title displayed above the plot (e.g. "Chr 12").

    data_uri
        Base64-encoded PNG embedded as a data URI so the HTML report is fully
        self-contained.

    stem
        Original filename stem of the plot (used for stable IDs or debugging).

    has_cnv
        Whether the chromosome contains a detected CNV/LOH event. Used to
        group plots into "with CNV" and "without CNV" sections in the report.
    """

    title: str
    data_uri: str
    stem: str
    has_cnv: bool


from pathlib import Path
from typing import Any


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
    """
    empty_groups: dict[str, list[dict[str, Any]]] = {
        "chr_with_cnv": [],
        "chr_no_cnv": [],
    }

    if chr_plots_dir is None:
        return empty_groups

    qc_dir = Path(chr_plots_dir)
    if not qc_dir.is_dir():
        return empty_groups

    groups: dict[str, list[PlotCard]] = {
        "chr_with_cnv": [],
        "chr_no_cnv": [],
    }

    png_files = sorted(
        qc_dir.glob("cnv_chr*segments.png"),
        key=lambda p: _chr_sort_key_from_stem(p.stem),
    )

    for png_file in png_files:
        stem = png_file.stem
        data_uri = _png_to_data_uri(png_file)
        if not data_uri:
            continue

        chr_label = _extract_chr_label_from_stem(stem)
        has_cnv = chr_label in cnv_chr_set
        title = f"Chr {chr_label}"

        card = PlotCard(
            title=title,
            data_uri=data_uri,
            stem=stem,
            has_cnv=has_cnv,
        )
        groups["chr_with_cnv" if has_cnv else "chr_no_cnv"].append(card)

    return {
        group_name: [card.__dict__ for card in cards]
        for group_name, cards in groups.items()
    }


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
    loh_segments_df = None
    if Path(loh_regions).is_file():
        loh_segments_df = load_purecn_segments(loh_regions)
    pon_df = None
    if pon:
        pon_df = load_pon_bins(pon)
    cytoband_df = load_cytobands(cytoband)

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
        loh_segments_df=loh_segments_df,
    )

    # ----------------------------
    # Create gene region table
    # ----------------------------
    generegions_df = build_generegion_table(
        cnr_df=cnr_df,
        cns_df=cns_df,
        sex=sex,
        pon_df=pon_df,
        cancer_genes=cancer_gene_set,
        loh_segments_df=loh_segments_df,
    )

    # --- 1) PureCN summary (from purity_csv) ---
    purecn_summary_df = read_purecn_summary(purity_csv)

    # --- 2) Extra QC / CNV metrics (DLR, PON spread) ---
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
    plot_config = EXOME_PLOT_CONFIG if is_exome else PANEL_PLOT_CONFIG

    # pon may be None; plot_chromosomes is written to handle that
    plot_chromosomes(
        cnr_df=cnr_df,
        vcf_path=vcf,
        segments_df=segments_df,
        generegions_df=generegions_df,
        outdir=chr_plots_dir,
        case_id=case_id,
        plot_config=plot_config,
        pon_df=pon_df,
    )

    if is_exome:
        generegions_df = generegions_df[
            generegions_df["is_cancer_gene"].fillna(False).astype(bool)
        ]

    if pon == None:
        generegions_df = generegions_df.drop(
            columns=[c for c in generegions_df.columns if c.startswith("pon")]
        )

    # ----------------------------
    # HTML report
    # ----------------------------
    out_html_path = str(output_file)

    render_cnv_report_html(
        df_segments=segments_df,
        df_regions=generegions_df,
        df_purecn_summary=purecn_summary_df,
        df_qc_summary=qc_summary_df,
        scatter_png=scatter_png_path,
        diagram_png=diagram_png_path,
        chr_plots_dir=chr_plots_dir,
        out_html=out_html_path,
        title=f"CNV Report – {case_id}",
        normalisation_method=normalisation_method,
        sample_sex=sex,
    )

    click.echo(f"[CNV QC] Finished report for {case_id}: {out_html_path}")


if __name__ == "__main__":
    sys.exit(main())
