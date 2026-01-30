from __future__ import annotations
import sys
import base64
import json
from pathlib import Path
from typing import Optional

import click
import numpy as np
import pandas as pd

from cnv_report_utils import (
    plot_chromosomes,
    build_gene_segment_table,
    build_gene_chunk_table,
    load_cancer_gene_set,
    compute_summary_metrics,
    pdf_first_page_to_png,
)
from BALSAMIC.constants.analysis import Gender

# =============================================================================
# Small helpers to reduce duplication and keep csv_to_html_table readable
# =============================================================================

def _png_to_data_uri(png_path: str | Path | None) -> str | None:
    """Read a PNG file and return a data URI, or None if missing."""
    if png_path is None:
        return None
    path = Path(png_path)
    if not path.is_file():
        return None
    png_bytes = path.read_bytes()
    png_base64 = base64.b64encode(png_bytes).decode("ascii")
    return f"data:image/png;base64,{png_base64}"


def _chr_sort_key_from_stem(stem: str) -> tuple[int, int]:
    """Sort by chr number (1..22) then X/Y, with unknowns last."""
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


def _extract_gene_names_from_stem(stem: str) -> list[str]:
    """Extract gene list from a gene-zoom plot stem, else empty list."""
    if "_genes_" not in stem:
        return []
    after_genes = stem.split("_genes_", 1)[1]
    gene_part = after_genes
    if gene_part.endswith("_segments"):
        gene_part = gene_part[: -len("_segments")]
    return [g for g in gene_part.split("_") if g]


def _as_upper_str(x) -> str:
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def _is_true_str(x) -> bool:
    return _as_upper_str(x) == "TRUE"


def _is_amp_del(x) -> bool:
    return _as_upper_str(x) in {"DELETION", "AMPLIFICATION"}


def _compute_cnv_sets(df: pd.DataFrame) -> tuple[set[str], set[str], str | None]:
    """Return (cnv_chr_set, cnv_gene_set, chr_col)."""
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


def _build_genome_plots_html(scatter_data_uri: str | None, diagram_data_uri: str | None) -> str:
    if not (scatter_data_uri or diagram_data_uri):
        return ""

    parts: list[str] = []
    parts.append("<h2>Genome-wide CNVkit plots</h2>")
    parts.append('<div style="margin-top:10px;">')

    if scatter_data_uri:
        parts.append(
            """
        <h3>CNVkit scatter</h3>
        <div style=\"border:1px solid #ccc; padding:10px; margin-top:10px; width:100%;\">
          <img
            src=\"{scatter_src}\"
            alt=\"CNVkit scatter PNG plot\"
            class=\"clickable-plot\"
            style=\"max-width:100%; height:auto; display:block; cursor:pointer;\"
          />
        </div>
            """.format(
                scatter_src=scatter_data_uri
            )
        )

    if diagram_data_uri:
        parts.append(
            """
        <h3>CNVkit diagram</h3>
        <div style=\"border:1px solid #ccc; padding:10px; margin-top:10px; width:100%;\">
          <img
            src=\"{diagram_src}\"
            alt=\"CNVkit diagram PNG plot\"
            class=\"clickable-plot\"
            style=\"max-width:100%; height:auto; display:block; cursor:pointer;\"
          />
        </div>
            """.format(
                diagram_src=diagram_data_uri
            )
        )

    parts.append("</div>")
    return "\n".join(parts)


def _build_plot_card_html(
    *,
    title: str,
    img_data_uri: str,
    stem: str,
    has_cnv: bool,
) -> str:
    badge_html = '<span class="plot-badge">CNV</span>' if has_cnv else ""
    return f"""
                <div class=\"plot-card {'plot-card-cnv' if has_cnv else 'plot-card-normal'}\">
                  <div class=\"plot-header\">
                    <span class=\"plot-title\">{title}</span>
                    {badge_html}
                  </div>
                  <img
                    src=\"{img_data_uri}\"
                    alt=\"QC plot {stem}\"
                    class=\"clickable-plot\"
                    style=\"max-width:100%; height:auto; display:block; cursor:pointer;\"
                  />
                </div>
    """


def _collect_qc_plot_blocks(
    qc_dir: Path,
    chr_col: str | None,
    cnv_chr_set: set[str],
    cnv_gene_set: set[str],
) -> tuple[list[str], list[str], list[str], list[str]]:
    """Return (chr_with, chr_no, gene_with, gene_no) blocks."""
    chr_plot_blocks_with_cnv: list[str] = []
    chr_plot_blocks_no_cnv: list[str] = []
    gene_plot_blocks_with_cnv: list[str] = []
    gene_plot_blocks_no_cnv: list[str] = []

    png_files = sorted(list(qc_dir.glob("cnv_chr*segments.png")), key=lambda p: _chr_sort_key_from_stem(p.stem))

    for png_file in png_files:
        img_data_uri = _png_to_data_uri(png_file)
        if not img_data_uri:
            continue

        stem = png_file.stem
        chr_label = _extract_chr_label_from_stem(stem)
        gene_names = _extract_gene_names_from_stem(stem)
        is_gene_plot = bool(gene_names)

        if not is_gene_plot:
            has_cnv = bool(chr_col is not None and chr_label in cnv_chr_set)
            title = f"Chr {chr_label}"
        else:
            has_cnv = any(g in cnv_gene_set for g in gene_names)
            title = f"Chr {chr_label} – {', '.join(gene_names)}"

        block = _build_plot_card_html(title=title, img_data_uri=img_data_uri, stem=stem, has_cnv=has_cnv)

        if is_gene_plot:
            (gene_plot_blocks_with_cnv if has_cnv else gene_plot_blocks_no_cnv).append(block)
        else:
            (chr_plot_blocks_with_cnv if has_cnv else chr_plot_blocks_no_cnv).append(block)

    return chr_plot_blocks_with_cnv, chr_plot_blocks_no_cnv, gene_plot_blocks_with_cnv, gene_plot_blocks_no_cnv


def _read_purecn_summary(purity_csv: str | Path | None) -> pd.DataFrame | None:
    if not purity_csv:
        return None
    df = pd.read_csv(purity_csv)
    wanted_cols = [
        "Sampleid",
        "Purity",
        "Ploidy",
        "Sex",
        "Contamination",
        "Flagged",
        "Failed",
        "Curated",
        "Comment",
    ]
    keep = [c for c in wanted_cols if c in df.columns]
    return df[keep] if keep else df

def csv_to_html_table(
    df: pd.DataFrame,                     # gene-level table
    out_html: str | Path,
    scatter_png: str | Path | None = None,
    diagram_png: str | Path | None = None,
    chr_plots_dir: str | Path | None = None,
    purecn_summary_df: pd.DataFrame | None = None,
    qc_summary_df: pd.DataFrame | None = None,
    chunk_df: pd.DataFrame | None = None,
) -> None:
    out_path = Path(out_html)

    scatter_png_data_uri = _png_to_data_uri(scatter_png)
    diagram_png_data_uri = _png_to_data_uri(diagram_png)

    genome_plots_html = _build_genome_plots_html(scatter_png_data_uri, diagram_png_data_uri)

    cnv_chr_set, cnv_gene_set, chr_col = _compute_cnv_sets(df)

    chr_plot_blocks_with_cnv: list[str] = []
    chr_plot_blocks_no_cnv: list[str] = []
    gene_plot_blocks_with_cnv: list[str] = []
    gene_plot_blocks_no_cnv: list[str] = []

    if chr_plots_dir is not None:
        qc_dir = Path(chr_plots_dir)
        if qc_dir.is_dir():
            (
                chr_plot_blocks_with_cnv,
                chr_plot_blocks_no_cnv,
                gene_plot_blocks_with_cnv,
                gene_plot_blocks_no_cnv,
            ) = _collect_qc_plot_blocks(qc_dir, chr_col, cnv_chr_set, cnv_gene_set)

    # Build visible chr & gene plots HTML
    qc_plots_html = ""

    # Chromosome-level plots
    if chr_plot_blocks_with_cnv or chr_plot_blocks_no_cnv:
        section_html = """
        <h2>Per-chromosome CNV plots</h2>
        <p class="muted">
          These plots show log2 coverage, PON spread (if available), BAF and segments/genes
          per chromosome. Plots containing CNV / LOH genes are highlighted.
          Click any plot to view it in a larger overlay.
        </p>
        """
        if chr_plot_blocks_with_cnv:
            section_html += "<h3>Chromosomes with CNV / LOH</h3>\n"
            section_html += (
                '<div class="plot-grid">'
                + "\n".join(chr_plot_blocks_with_cnv)
                + "</div>"
            )
        if chr_plot_blocks_no_cnv:
            section_html += "<h3>Chromosomes without CNV / LOH</h3>\n"
            section_html += (
                '<div class="plot-grid">' + "\n".join(chr_plot_blocks_no_cnv) + "</div>"
            )
        qc_plots_html += section_html

    # Gene-level plots
    if gene_plot_blocks_with_cnv or gene_plot_blocks_no_cnv:
        section_html = """
        <h2>Gene-level CNV plots</h2>
        <p class="muted">
          Zoomed plots centred on specific genes or gene intervals.
          Plots containing CNV / LOH in any of the genes are highlighted.
          Click any plot to view it in a larger overlay.
        </p>
        """
        if gene_plot_blocks_with_cnv:
            section_html += "<h3>Regions with CNV / LOH</h3>\n"
            section_html += (
                '<div class="plot-grid">'
                + "\n".join(gene_plot_blocks_with_cnv)
                + "</div>"
            )
        if gene_plot_blocks_no_cnv:
            section_html += "<h3>Regions without CNV / LOH</h3>\n"
            section_html += (
                '<div class="plot-grid">'
                + "\n".join(gene_plot_blocks_no_cnv) + "</div>"
            )
        qc_plots_html += section_html

    # ---------------- PureCN sample summary table ----------------
    purecn_summary_html = ""
    if purecn_summary_df is not None and not purecn_summary_df.empty:
        purecn_summary_html = purecn_summary_df.to_html(
            index=False,
            border=0,
            classes="dataframe",
            table_id="purecn-summary-table",
        )

    # ---------------- Extra QC / CNV summary table ----------------
    qc_summary_html = ""
    if qc_summary_df is not None and not qc_summary_df.empty:
        qc_display = qc_summary_df.copy()
        # Round floats for readability
        for col in qc_display.columns:
            if np.issubdtype(qc_display[col].dtype, np.floating):
                qc_display[col] = qc_display[col].round(3)

        qc_summary_html = qc_display.to_html(
            index=False,
            border=0,
            classes="dataframe",
            table_id="qc-summary-table",
        )

    # ---------------- Gene-level CNV/LOH table ----------------
    gene_table_html = df.to_html(
        index=False,
        border=0,
        classes="dataframe",
        table_id="report-table",   # gene-level
    )

    # ---------------- Chunk-level table (optional) ----------------
    if chunk_df is not None and not chunk_df.empty:
        chunk_table_html = chunk_df.to_html(
            index=False,
            border=0,
            classes="dataframe",
            table_id="chunk-table",
        )
        # Chunk-level controls (includes PON + split options)
        chunk_section_html = f"""
        <h2>Within-gene PON-driven chunks (optional detail)</h2>
        <p class="muted">
          These rows represent sub-gene chunks defined by PON-based deviation (log2 vs panel-of-normals),
          for genes in the main table. Filters below apply only to this chunk-level table.
        </p>
        <button id="toggle-chunk-table" type="button">Show chunk-level table</button>
        <div id="chunk-table-container" style="display:none; margin-top:10px;">
          <div class="controls">
            <div class="control-group">
              <label>
                <input type="checkbox" id="hide-non-cnv-chunk">
                Show only LOH / CNV chunks (loh_flag TRUE or CNV call ≠ NEUTRAL)
              </label>
              <label style="margin-left:12px;">
                <input type="checkbox" id="show-pon-cnv-chunk">
                Also include chunks with PON-based CNV call (pon_cnv_call ≠ "")
              </label>
              <label style="margin-left:12px;">
                <input type="checkbox" id="only-cancer-genes-chunk">
                Show only cancer genes (is_cancer_gene TRUE)
              </label>
              <label style="margin-left:12px;">
                <input type="checkbox" id="only-split-genes-chunk">
                Show only split genes (is_gene_split TRUE)
              </label>
            </div>
            <div class="control-group">
              <label>
                Min n.targets:
                <input type="number" id="min-targets-chunk" min="0" step="1" value="0">
              </label>
            </div>
            <div class="control-group">
              <label>
                Gene filter (comma-separated gene.symbol list):
                <input type="text" id="gene-filter-chunk" placeholder="TP53,EGFR,MYC">
              </label>
            </div>
          </div>
          {chunk_table_html}
        </div>
        """
    else:
        chunk_section_html = ""  # no chunks

    # Column indices for JS
    col_idx_map = {name: idx for idx, name in enumerate(df.columns)}
    col_idx_json = json.dumps(col_idx_map)

    if chunk_df is not None and not chunk_df.empty:
        chunk_col_idx_map = {name: idx for idx, name in enumerate(chunk_df.columns)}
        chunk_col_idx_json = json.dumps(chunk_col_idx_map)
    else:
        chunk_col_idx_json = "{}"

    html = f"""<!doctype html>
    <html lang="en">
    <head>
      <meta charset="utf-8">
      <title>CNV Report</title>
      <style>
        body {{
          font-family: system-ui, sans-serif;
          margin: 24px;
        }}
        h1, h2, h3 {{
          font-weight: 600;
        }}
        .controls {{
          margin-bottom: 12px;
          padding: 8px 10px;
          border: 1px solid #ddd;
          background: #fafafa;
        }}
        .control-group {{
          margin-bottom: 6px;
        }}
        table.dataframe {{
          border-collapse: collapse;
          width: 100%;
          font-size: 14px;
          margin-bottom: 16px;
        }}
        table.dataframe th, table.dataframe td {{
          border: 1px solid #ddd;
          padding: 6px 8px;
        }}
        table.dataframe th {{
          background: #f4f4f4;
          text-align: left;
          position: sticky;
          top: 0;
          z-index: 1;
        }}
        table.dataframe tr:nth-child(even) {{
          background: #fafafa;
        }}
        .muted {{
          color: #666;
          font-size: 13px;
        }}
        input[type="text"], input[type="number"] {{
          padding: 4px 6px;
          font-size: 13px;
        }}
        button {{
          padding: 4px 10px;
          font-size: 13px;
          cursor: pointer;
        }}
        label {{
          font-size: 13px;
        }}
        .plot-grid {{
          display: flex;
          flex-wrap: wrap;
          gap: 16px;
          margin-top: 8px;
        }}
        .plot-card {{
          border-radius: 4px;
          border: 1px solid #ccc;
          padding: 8px;
          background: #fff;
          flex: 0 0 320px;
          max-width: 380px;
        }}
        .plot-card-cnv {{
          border-color: #c0392b;
          box-shadow: 0 0 4px rgba(192, 57, 43, 0.4);
        }}
        .plot-header {{
          display: flex;
          align-items: center;
          justify-content: space-between;
          margin-bottom: 6px;
        }}
        .plot-title {{
          font-weight: 600;
          font-size: 14px;
        }}
        .plot-badge {{
          display: inline-block;
          padding: 2px 6px;
          border-radius: 12px;
          background: #c0392b;
          color: #fff;
          font-size: 11px;
          font-weight: 600;
        }}
        .modal-overlay {{
          display: none;
          position: fixed;
          z-index: 9999;
          left: 0;
          top: 0;
          width: 100%;
          height: 100%;
          background: rgba(0, 0, 0, 0.85);
        }}
        .modal-content-wrapper {{
          position: absolute;
          top: 50%;
          left: 50%;
          transform: translate(-50%, -50%);
          max-width: 98vw;
          max-height: 98vh;
        }}
        .modal-image {{
          max-width: 98vw;
          max-height: 98vh;
          width: auto;
          height: auto;
          display: block;
          border-radius: 4px;
          box-shadow: 0 0 12px rgba(0,0,0,0.7);
        }}
        .modal-close {{
          position: absolute;
          top: -32px;
          right: 0;
          color: #fff;
          font-size: 24px;
          font-weight: 700;
          cursor: pointer;
        }}
      </style>
    </head>
    <body>
      <h1>CNV Report</h1>

      <h2>PureCN sample summary</h2>
      {purecn_summary_html}

      <h2>CNV / QC summary</h2>
      <p class="muted">
        Additional metrics derived from CNVkit CNR/CNN and gene-level segments,
        including DLR-like noise, standard deviation of adjacent Target Log2-ratios, and PON spread summaries.
      </p>
      {qc_summary_html}

      {genome_plots_html}

      {qc_plots_html}

      <h1>LOH / CNV Genes (gene-level)</h1>
      <button id="toggle-gene-table" type="button">Hide gene-level table</button>
      <div id="gene-table-container" style="margin-top:10px;">
        <div class="controls">
          <div class="control-group">
            <label>
              <input type="checkbox" id="hide-non-cnv">
              Show only LOH / CNV genes (loh_flag TRUE or CNV call ≠ NEUTRAL)
            </label>
            <label style="margin-left:12px;">
              <input type="checkbox" id="only-cancer-genes">
              Show only cancer genes (is_cancer_gene TRUE)
            </label>
          </div>
          <div class="control-group">
            <label>
              Min n.targets:
              <input type="number" id="min-targets" min="0" step="1" value="0">
            </label>
          </div>
          <div class="control-group">
            <label>
              Gene filter (comma-separated gene.symbol list):
              <input type="text" id="gene-filter" placeholder="TP53,EGFR,MYC">
            </label>
          </div>
        </div>

        {gene_table_html}
      </div>

      {chunk_section_html}

      <div id="plot-modal" class="modal-overlay">
        <div class="modal-content-wrapper">
          <span class="modal-close" id="plot-modal-close">&times;</span>
          <img id="plot-modal-image" class="modal-image" src="" alt="Enlarged plot">
        </div>
      </div>

      <script>
        (function() {{
          const geneTable  = document.getElementById("report-table");
          const chunkTable = document.getElementById("chunk-table");

          // Gene-level controls (no PON / split)
          const gHideNonCnv  = document.getElementById("hide-non-cnv");
          const gOnlyCancer  = document.getElementById("only-cancer-genes");
          const gGeneInput   = document.getElementById("gene-filter");
          const gMinTargets  = document.getElementById("min-targets");

          // Chunk-level controls (includes PON + split)
          const cHideNonCnv    = document.getElementById("hide-non-cnv-chunk");
          const cIncludePonCnv = document.getElementById("show-pon-cnv-chunk");
          const cOnlyCancer    = document.getElementById("only-cancer-genes-chunk");
          const cOnlySplit     = document.getElementById("only-split-genes-chunk");
          const cGeneInput     = document.getElementById("gene-filter-chunk");
          const cMinTargets    = document.getElementById("min-targets-chunk");

          const colIndexGene  = {col_idx_json};
          const colIndexChunk = {chunk_col_idx_json};

          function normalizeCellText(td) {{
            return (td && td.textContent ? td.textContent : "")
              .trim()
              .toLowerCase();
          }}

          function isCnvRow(row, colIndex) {{
            const lohColIdx    = (colIndex["loh_flag"] ?? -1);
            const cnvkitColIdx = (colIndex["cnvkit_cnv_call"] ?? -1);
            const purecnColIdx = (colIndex["purecn_cnv_call"] ?? -1);

            let isCnv = false;

            if (lohColIdx !== -1) {{
              const lohCell = row.cells[lohColIdx];
              const lohVal = normalizeCellText(lohCell);
              if (lohVal === "true") {{
                isCnv = true;
              }}
            }}

            if (!isCnv && cnvkitColIdx !== -1) {{
              const cCell = row.cells[cnvkitColIdx];
              const cVal = normalizeCellText(cCell);
              if (cVal === "deletion" || cVal === "amplification") {{
                isCnv = true;
              }}
            }}

            if (!isCnv && purecnColIdx !== -1) {{
              const pCell = row.cells[purecnColIdx];
              const pVal = normalizeCellText(pCell);
              if (pVal === "deletion" || pVal === "amplification") {{
                isCnv = true;
              }}
            }}

            return isCnv;
          }}

          function isCancerGeneRow(row, colIndex) {{
            const cancerColIdx = (colIndex["is_cancer_gene"] ?? -1);
            if (cancerColIdx === -1) return false;
            const cell = row.cells[cancerColIdx];
            const val = normalizeCellText(cell);
            return (val === "true" || val === "1" || val === "yes");
          }}

          function isSplitGeneRow(row, colIndex) {{
            const splitColIdx = (colIndex["is_gene_split"] ?? -1);
            if (splitColIdx === -1) return false;
            const cell = row.cells[splitColIdx];
            const val = normalizeCellText(cell);
            return (val === "true" || val === "1" || val === "yes");
          }}

          function parseGeneList(value) {{
            return (value || "")
              .split(",")
              .map(g => g.trim().toLowerCase())
              .filter(g => g.length > 0);
          }}

          function filterTable(table, colIndex, cfg) {{
            if (!table) return;

            const geneColIdx     = (colIndex["gene.symbol"] ?? -1);
            const nTargetsColIdx = (colIndex["n.targets"] ?? -1);
            const ponCnvColIdx   = (colIndex["pon_cnv_call"] ?? -1);

            const body = table.tBodies[0];
            if (!body) return;
            const rows = Array.from(body.rows);

            for (const row of rows) {{
              let hideRow = false;

              if (!hideRow && cfg.hideNonCnv) {{
                const cnv = isCnvRow(row, colIndex);

                let isPonCnv = false;
                if (ponCnvColIdx !== -1) {{
                  const ponCell = row.cells[ponCnvColIdx];
                  const ponVal = normalizeCellText(ponCell);
                  if (ponVal === "amplification" || ponVal === "deletion") {{
                    isPonCnv = true;
                  }}
                }}

                const treatedAsCnv = cnv || (cfg.includePonCnv && isPonCnv);
                if (!treatedAsCnv) {{
                  hideRow = true;
                }}
              }}

              if (!hideRow && cfg.onlyCancer) {{
                if (!isCancerGeneRow(row, colIndex)) {{
                  hideRow = true;
                }}
              }}

              if (!hideRow && cfg.onlySplit) {{
                if (!isSplitGeneRow(row, colIndex)) {{
                  hideRow = true;
                }}
              }}

              if (!hideRow && nTargetsColIdx !== -1 && cfg.minTargets > 0) {{
                const nCell = row.cells[nTargetsColIdx];
                const nText = normalizeCellText(nCell);
                if (nText.length > 0) {{
                  const nVal = parseFloat(nText);
                  if (!Number.isNaN(nVal) && nVal < cfg.minTargets) {{
                    hideRow = true;
                  }}
                }}
              }}

              if (!hideRow && geneColIdx !== -1 && cfg.geneList.length > 0) {{
                const geneCell = row.cells[geneColIdx];
                const geneText = normalizeCellText(geneCell);
                const matches = cfg.geneList.some(g => geneText === g);
                if (!matches) {{
                  hideRow = true;
                }}
              }}

              row.style.display = hideRow ? "none" : "";
            }}
          }}

          function filterGeneTable() {{
            if (!geneTable) return;

            const cfg = {{
              hideNonCnv:  gHideNonCnv && gHideNonCnv.checked,
              includePonCnv: false,  // gene-level ignores PON-based CNV override
              onlyCancer:  gOnlyCancer && gOnlyCancer.checked,
              onlySplit:   false,    // gene-level ignores split-gene filter
              geneList:    parseGeneList(gGeneInput && gGeneInput.value),
              minTargets:  gMinTargets ? (parseInt(gMinTargets.value, 10) || 0) : 0,
            }};
            filterTable(geneTable, colIndexGene, cfg);
          }}

          function filterChunkTable() {{
            if (!chunkTable) return;

            const cfg = {{
              hideNonCnv:   cHideNonCnv && cHideNonCnv.checked,
              includePonCnv: cIncludePonCnv && cIncludePonCnv.checked,
              onlyCancer:   cOnlyCancer && cOnlyCancer.checked,
              onlySplit:    cOnlySplit && cOnlySplit.checked,
              geneList:     parseGeneList(cGeneInput && cGeneInput.value),
              minTargets:   cMinTargets ? (parseInt(cMinTargets.value, 10) || 0) : 0,
            }};
            filterTable(chunkTable, colIndexChunk, cfg);
          }}

          function applyFilter() {{
            filterGeneTable();
            filterChunkTable();
          }}

          // Hook gene-level controls
          if (gHideNonCnv)  gHideNonCnv.addEventListener("change", applyFilter);
          if (gOnlyCancer)  gOnlyCancer.addEventListener("change", applyFilter);
          if (gGeneInput)   gGeneInput.addEventListener("input", applyFilter);
          if (gMinTargets)  gMinTargets.addEventListener("input", applyFilter);

          // Hook chunk-level controls
          if (cHideNonCnv)    cHideNonCnv.addEventListener("change", applyFilter);
          if (cIncludePonCnv) cIncludePonCnv.addEventListener("change", applyFilter);
          if (cOnlyCancer)    cOnlyCancer.addEventListener("change", applyFilter);
          if (cOnlySplit)     cOnlySplit.addEventListener("change", applyFilter);
          if (cGeneInput)     cGeneInput.addEventListener("input", applyFilter);
          if (cMinTargets)    cMinTargets.addEventListener("input", applyFilter);

          // Initial filter
          applyFilter();

          // --- Gene-level table toggle ---
          const geneToggleBtn   = document.getElementById("toggle-gene-table");
          const geneContainer   = document.getElementById("gene-table-container");
          if (geneToggleBtn && geneContainer) {{
            let geneVisible = true;
            geneToggleBtn.addEventListener("click", () => {{
              geneVisible = !geneVisible;
              geneContainer.style.display = geneVisible ? "block" : "none";
              geneToggleBtn.textContent = geneVisible
                ? "Hide gene-level table"
                : "Show gene-level table";
            }});
          }}

          // --- Chunk-level table toggle ---
          const chunkToggleBtn   = document.getElementById("toggle-chunk-table");
          const chunkContainer   = document.getElementById("chunk-table-container");
          if (chunkToggleBtn && chunkContainer) {{
            let chunkVisible = false;
            chunkToggleBtn.addEventListener("click", () => {{
              chunkVisible = !chunkVisible;
              chunkContainer.style.display = chunkVisible ? "block" : "none";
              chunkToggleBtn.textContent = chunkVisible
                ? "Hide chunk-level table"
                : "Show chunk-level table";
            }});
          }}

          // --- Modal logic for enlarged plots ---
          const modal       = document.getElementById("plot-modal");
          const modalImg    = document.getElementById("plot-modal-image");
          const modalClose  = document.getElementById("plot-modal-close");
          const clickableImgs = document.querySelectorAll(".clickable-plot");

          function openModal(src) {{
            if (!modal || !modalImg) return;
            modalImg.src = src;
            modal.style.display = "block";
          }}

          function closeModal() {{
            if (!modal) return;
            modal.style.display = "none";
            if (modalImg) {{
              modalImg.src = "";
            }}
          }}

          clickableImgs.forEach(img => {{
            img.addEventListener("click", () => openModal(img.src));
          }});

          if (modalClose) {{
            modalClose.addEventListener("click", closeModal);
          }}

          if (modal) {{
            modal.addEventListener("click", (event) => {{
              if (event.target === modal) {{
                closeModal();
              }}
            }});
          }}

          document.addEventListener("keydown", (event) => {{
            if (event.key === "Escape") {{
              closeModal();
            }}
          }});
        }})();
      </script>
    </body>
    </html>
    """
    out_path.write_text(html, encoding="utf-8")


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
    default=False,
    help="Sex for how to handle sex chromosomes.",
)
def main(
    loh_genes,
    cnr,
    cns,
    pon,
    vcf,
    refgene,
    cytoband,
    case_id,
    cust_case_id,
    cnvkit_scatter,
    cnvkit_diagram,
    out_prefix,
    purity_csv,
    cancer_genes,
    is_exome,
    sex,
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
    cancer_gene_set = None
    if cancer_genes:
        # You can tune min_occurrence depending on exome/panel if you like
        min_occ = 3
        cancer_gene_set = load_cancer_gene_set(
            cancer_genes,
            min_occurrence=min_occ,
            only_annotated=True,
        )

    # ----------------------------
    # Create per-gene CNV table
    # ----------------------------
    # pon can be None (no PON file); build_gene_segment_table handles that
    genes_df = build_gene_segment_table(
        cnr_path=cnr,
        cns_path=cns,
        cancer_genes=cancer_gene_set,
        refgene_path=refgene,
        transcript_selection="longest_tx",
        loh_path=loh_genes,
        cytoband_path=cytoband,
        sex=sex,
        pon_path=pon,
    )

    # ----------------------------
    # Create per chunk table
    # ----------------------------
    chunk_df = build_gene_chunk_table(cnr_path=cnr, gene_seg_df=genes_df, pon_path=pon)


    # --- 1) PureCN summary (from purity_csv) ---
    purecn_summary_df = _read_purecn_summary(purity_csv)

    # --- 2) Extra QC / CNV metrics (DLR, PON spread, chunk stats) ---
    qc_summary_df = compute_summary_metrics(
        cnr_path=cnr,
        cnn_path=pon,  # CNVkit .cnn PON, or None
        gene_seg_df=genes_df,  # from build_gene_segment_table
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
        gene_seg_df=genes_df,
        gene_chunk_df=chunk_df,
        outdir=chr_plots_dir,
        case_id=plot_case_id,
        pon_path=pon,
        include_y=True,
        neutral_target_factor=neutral_target_factor,
        highlight_only_cancer=highlight_only_cancer,
    )

    plot_chromosomes(
        cnr_path=cnr,
        vcf_path=vcf,
        gene_seg_df=genes_df,
        gene_chunk_df=chunk_df,
        outdir=chr_plots_dir,
        case_id=plot_case_id,
        pon_path=pon,
        include_y=True,
        neutral_target_factor=neutral_target_factor,
        highlight_only_cancer=highlight_only_cancer,
        focus_genes=["RB1", "DLEU1", "TP53", "KMT2D"],
        focus_padding_bp=100_000,
    )

    if is_exome:
        if "is_cancer_gene" in genes_df.columns:
            genes_df = genes_df[genes_df["is_cancer_gene"].fillna(False).astype(bool)]
        if chunk_df is not None and not chunk_df.empty and "is_cancer_gene" in chunk_df.columns:
            chunk_df = chunk_df[chunk_df["is_cancer_gene"].fillna(False).astype(bool)]

    # ----------------------------
    # HTML report
    # ----------------------------
    out_html = str(out_prefix)
    csv_to_html_table(
        df=genes_df,
        out_html=out_html,
        scatter_png=str(scatter_png_path),
        diagram_png=str(diagram_png_path),
        chr_plots_dir=chr_plots_dir,
        purecn_summary_df=purecn_summary_df,
        qc_summary_df=qc_summary_df,
        chunk_df=chunk_df,
    )

    click.echo(f"[CNV QC] Finished report for {case_id}: {out_html}")


if __name__ == "__main__":
    sys.exit(main())
