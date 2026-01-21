import sys
import click
import base64
import json
import pandas as pd
from pathlib import Path
from cnv_report_utils import plot_chromosomes, build_gene_segment_table, pdf_first_page_to_png, load_cancer_gene_set
from BALSAMIC.constants.analysis import Gender
def csv_to_html_table(
    df: pd.DataFrame,
    out_html: str,
    scatter_png: str,
    diagram_png: str,
    chr_plots_dir: str | Path | None = None,
    summary_df: pd.DataFrame | None = None,
) -> None:
    out_path = Path(out_html)

    # --- PNG → Base64 (CNVkit scatter) ---
    png_bytes = Path(scatter_png).read_bytes()
    png_base64 = base64.b64encode(png_bytes).decode("ascii")
    scatter_png_data_uri = f"data:image/png;base64,{png_base64}"

    # --- PNG → Base64 (CNVkit diagram) ---
    png_bytes = Path(diagram_png).read_bytes()
    png_base64 = base64.b64encode(png_bytes).decode("ascii")
    diagram_png_data_uri = f"data:image/png;base64,{png_base64}"

    # ------------------------------------------------------------------
    # Determine which chromosomes & genes have CNV / LOH
    # ------------------------------------------------------------------
    df_chr = df.copy()
    if "chr" in df_chr.columns:
        chr_col = "chr"
    elif "chromosome" in df_chr.columns:
        chr_col = "chromosome"
    else:
        chr_col = None

    cnv_chr_set: set[str] = set()
    cnv_gene_set: set[str] = set()

    if chr_col is not None:
        if "gene.symbol" in df_chr.columns:
            df_chr = df_chr[~df_chr["gene.symbol"].isin(["Antitarget", "-"])]

        cnv_mask = pd.Series(False, index=df_chr.index)

        # LOH flag
        if "loh_flag" in df_chr.columns:
            loh_str = df_chr["loh_flag"].astype(str).str.upper()
            cnv_mask |= (loh_str == "TRUE")

        # CNVkit call
        if "cnvkit_cnv_call" in df_chr.columns:
            cnvkit_str = df_chr["cnvkit_cnv_call"].astype(str).str.upper()
            cnv_mask |= cnvkit_str.isin(["DELETION", "AMPLIFICATION"])

        # PureCN call
        if "purecn_cnv_call" in df_chr.columns:
            purecn_str = df_chr["purecn_cnv_call"].astype(str).str.upper()
            cnv_mask |= purecn_str.isin(["DELETION", "AMPLIFICATION"])

        cnv_chr_set = set(df_chr.loc[cnv_mask, chr_col].astype(str).tolist())

        if "gene.symbol" in df_chr.columns:
            cnv_gene_set = set(
                df_chr.loc[cnv_mask, "gene.symbol"].astype(str).tolist()
            )

    # -------------------------------------------------------------
    # Per-chromosome & gene-level PNGs (sorted, CNV highlighted)
    # -------------------------------------------------------------
    chr_plot_blocks_with_cnv: list[str] = []
    chr_plot_blocks_no_cnv: list[str] = []

    gene_plot_blocks_with_cnv: list[str] = []
    gene_plot_blocks_no_cnv: list[str] = []

    if chr_plots_dir is not None:
        qc_dir = Path(chr_plots_dir)
        if qc_dir.is_dir():
            png_files = list(qc_dir.glob("cnv_chr*segments.png"))

            def chr_sort_key(path: Path):
                stem = path.stem  # e.g. cnv_chr13_genes_RB1_DLEU1_segments
                label = stem
                if "chr" in stem:
                    label = stem.split("chr", 1)[1]
                if "_" in label:
                    label = label.split("_", 1)[0]
                if label.isdigit():
                    return (0, int(label))
                special = {"X": 23, "Y": 24}
                return (1, special.get(label.upper(), 25))

            png_files = sorted(png_files, key=chr_sort_key)

            for png_file in png_files:
                img_bytes = png_file.read_bytes()
                img_b64 = base64.b64encode(img_bytes).decode("ascii")
                img_data_uri = f"data:image/png;base64,{img_b64}"

                stem = png_file.stem
                # Extract chromosome label
                chr_label = stem
                if "chr" in stem:
                    chr_label = stem.split("chr", 1)[1]
                if "_" in chr_label:
                    chr_label = chr_label.split("_", 1)[0]

                # Detect whether this is a gene-level zoom plot
                # Expected pattern: cnv_chr{chr}_genes_{GENE1[_GENE2…]}_segments
                is_gene_plot = "_genes_" in stem

                gene_names: list[str] = []
                gene_title = ""
                if is_gene_plot:
                    # e.g. "cnv_chr13_genes_RB1_DLEU1_segments"
                    after_genes = stem.split("_genes_", 1)[1]  # "RB1_DLEU1_segments"
                    gene_part = after_genes
                    if gene_part.endswith("_segments"):
                        gene_part = gene_part[: -len("_segments")]
                    # Genes were joined by "_"
                    gene_names = [g for g in gene_part.split("_") if g]
                    gene_title = ", ".join(gene_names)

                # Does this plot contain any CNV / LOH?
                has_cnv = False
                if not is_gene_plot:
                    # chromosome-level plot
                    if chr_col is not None and chr_label in cnv_chr_set:
                        has_cnv = True
                else:
                    # gene-level plot: any gene in cnv_gene_set?
                    if any(g in cnv_gene_set for g in gene_names):
                        has_cnv = True

                # Build plot card with highlighting
                if not is_gene_plot:
                    title = f"Chr {chr_label}"
                else:
                    title = f"Chr {chr_label} – {gene_title}"

                badge_html = '<span class="plot-badge">CNV</span>' if has_cnv else ""

                block = f"""
                <div class="plot-card {'plot-card-cnv' if has_cnv else 'plot-card-normal'}">
                  <div class="plot-header">
                    <span class="plot-title">{title}</span>
                    {badge_html}
                  </div>
                  <img
                    src="{img_data_uri}"
                    alt="QC plot {stem}"
                    style="max-width:100%; height:auto; display:block;"
                  />
                </div>
                """

                if is_gene_plot:
                    if has_cnv:
                        gene_plot_blocks_with_cnv.append(block)
                    else:
                        gene_plot_blocks_no_cnv.append(block)
                else:
                    if has_cnv:
                        chr_plot_blocks_with_cnv.append(block)
                    else:
                        chr_plot_blocks_no_cnv.append(block)

    # Build chromosome & gene plots HTML (all collapsible)
    qc_plots_html = ""

    # Chromosome-level plots
    if chr_plot_blocks_with_cnv or chr_plot_blocks_no_cnv:
        section_html = """
        <details>
          <summary class="collapsible-header">
            Per-chromosome CNV plots
          </summary>
          <p class="muted">
            These plots show log2 coverage, PON spread (if available), BAF and segments/genes
            per chromosome. Plots containing CNV / LOH genes are highlighted.
          </p>
        """
        if chr_plot_blocks_with_cnv:
            section_html += "<h3>Chromosomes with CNV / LOH</h3>\n"
            section_html += '<div class="plot-grid">' + "\n".join(chr_plot_blocks_with_cnv) + "</div>"
        if chr_plot_blocks_no_cnv:
            section_html += """
            <details style="margin-top:10px;">
              <summary class="collapsible-subheader">
                Show chromosomes without CNV / LOH
              </summary>
            """
            section_html += '<div class="plot-grid">' + "\n".join(chr_plot_blocks_no_cnv) + "</div>"
            section_html += "</details>"
        section_html += "</details>"
        qc_plots_html += section_html

    # Gene-level plots
    if gene_plot_blocks_with_cnv or gene_plot_blocks_no_cnv:
        section_html = """
        <details>
          <summary class="collapsible-header">
            Gene-level CNV plots
          </summary>
          <p class="muted">
            Zoomed plots centred on specific genes or gene intervals.
            Plots containing CNV / LOH in any of the genes are highlighted.
          </p>
        """
        if gene_plot_blocks_with_cnv:
            section_html += "<h3>Regions with CNV / LOH</h3>\n"
            section_html += '<div class="plot-grid">' + "\n".join(gene_plot_blocks_with_cnv) + "</div>"
        if gene_plot_blocks_no_cnv:
            section_html += """
            <details style="margin-top:10px;">
              <summary class="collapsible-subheader">
                Show regions without CNV / LOH
              </summary>
            """
            section_html += '<div class="plot-grid">' + "\n".join(gene_plot_blocks_no_cnv) + "</div>"
            section_html += "</details>"
        section_html += "</details>"
        qc_plots_html += section_html

    # ---------------- Summary table ----------------
    summary_html = ""
    if summary_df is not None and not summary_df.empty:
        summary_html = summary_df.to_html(
            index=False,
            border=0,
            classes="dataframe",
            table_id="summary-table",
        )

    # ---------------- CNV Gene table ----------------
    table_html = df.to_html(
        index=False,
        border=0,
        classes="dataframe",
        table_id="report-table",
    )

    # Column indices for JS
    col_idx_map = {name: idx for idx, name in enumerate(df.columns)}
    col_idx_json = json.dumps(col_idx_map)

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
        /* Plot layout & highlighting */
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
          flex: 1 1 320px;
          max-width: 520px;
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
        .collapsible-header {{
          cursor: pointer;
          font-weight: 600;
          font-size: 15px;
        }}
        .collapsible-subheader {{
          cursor: pointer;
          font-weight: 500;
          font-size: 14px;
        }}
        details {{
          margin-top: 16px;
        }}
      </style>
    </head>
    <body>
      <h1>CNV Report</h1>

      <h2>Sample summary</h2>
      {summary_html}

      <details open>
        <summary class="collapsible-header">
          Genome-wide CNVkit plots
        </summary>
        <div style="margin-top:10px;">
          <h3>CNVkit scatter</h3>
          <div style="border:1px solid #ccc; padding:10px; margin-top:10px; width:100%;">
            <img
              src="{scatter_png_data_uri}"
              alt="CNVkit scatter PNG plot"
              style="max-width:100%; height:auto; display:block;"
            />
          </div>

          <h3>CNVkit diagram</h3>
          <div style="border:1px solid #ccc; padding:10px; margin-top:10px; width:100%;">
            <img
              src="{diagram_png_data_uri}"
              alt="CNVkit diagram PNG plot"
              style="max-width:100%; height:auto; display:block;"
            />
          </div>
        </div>
      </details>

      {qc_plots_html}

      <h1>LOH / CNV Genes</h1>

      <div class="controls">
        <div class="control-group">
          <label>
            <input type="checkbox" id="hide-non-cnv">
            Show only LOH / CNV genes (loh_flag TRUE or CNV call ≠ NEUTRAL)
          </label>
          <label style="margin-left:12px;">
            <input type="checkbox" id="show-significant-pon">
            Also include genes with pon_call = "significant"
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

      {table_html}

      <script>
        (function() {{
          const table             = document.getElementById("report-table");
          const checkboxNonCnv    = document.getElementById("hide-non-cnv");
          const checkboxPonSig    = document.getElementById("show-significant-pon");
          const checkboxCancer    = document.getElementById("only-cancer-genes");
          const geneInput         = document.getElementById("gene-filter");
          const minTargetsInput   = document.getElementById("min-targets");

          if (!table) return;

          const colIndex = {col_idx_json};

          const geneColIdx        = colIndex["gene.symbol"] ?? -1;
          const lohColIdx         = colIndex["loh_flag"] ?? -1;
          const cnvkitColIdx      = colIndex["cnvkit_cnv_call"] ?? -1;
          const purecnColIdx      = colIndex["purecn_cnv_call"] ?? -1;
          const nTargetsColIdx    = colIndex["n.targets"] ?? -1;
          const ponCallColIdx     = colIndex["pon_call"] ?? -1;
          const cancerColIdx      = colIndex["is_cancer_gene"] ?? -1;

          function normalizeCellText(td) {{
            return (td && td.textContent ? td.textContent : "")
              .trim()
              .toLowerCase();
          }}

          function isCnvRow(row) {{
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

          function isCancerGeneRow(row) {{
            if (cancerColIdx === -1) return false;
            const cell = row.cells[cancerColIdx];
            const val = normalizeCellText(cell);
            return (val === "true" || val === "1" || val === "yes");
          }}

          function applyFilter() {{
            const hideNonCnv    = checkboxNonCnv && checkboxNonCnv.checked;
            const showPonSig    = checkboxPonSig && checkboxPonSig.checked;
            const onlyCancer    = checkboxCancer && checkboxCancer.checked;

            const geneList = (geneInput && geneInput.value ? geneInput.value : "")
              .split(",")
              .map(g => g.trim().toLowerCase())
              .filter(g => g.length > 0);

            const minTargets = minTargetsInput
              ? (parseInt(minTargetsInput.value, 10) || 0)
              : 0;

            const body = table.tBodies[0];
            if (!body) return;
            const rows = Array.from(body.rows);

            for (const row of rows) {{
              let hideRow = false;

              if (!hideRow && hideNonCnv) {{
                const cnv = isCnvRow(row);

                let isPonSignificant = false;
                if (ponCallColIdx !== -1) {{
                  const ponCell = row.cells[ponCallColIdx];
                  const ponVal = normalizeCellText(ponCell);
                  isPonSignificant = (ponVal === "significant");
                }}

                const treatedAsCnv = cnv || (showPonSig && isPonSignificant);
                if (!treatedAsCnv) {{
                  hideRow = true;
                }}
              }}

              if (!hideRow && onlyCancer && cancerColIdx !== -1) {{
                if (!isCancerGeneRow(row)) {{
                  hideRow = true;
                }}
              }}

              if (!hideRow && nTargetsColIdx !== -1 && minTargetsInput) {{
                const nCell = row.cells[nTargetsColIdx];
                const nText = normalizeCellText(nCell);
                if (nText.length > 0) {{
                  const nVal = parseFloat(nText);
                  if (!Number.isNaN(nVal) && nVal < minTargets) {{
                    hideRow = true;
                  }}
                }}
              }}

              if (!hideRow && geneColIdx !== -1 && geneList.length > 0) {{
                const geneCell = row.cells[geneColIdx];
                const geneText = normalizeCellText(geneCell);
                const matches = geneList.some(g => geneText === g);
                if (!matches) {{
                  hideRow = true;
                }}
              }}

              row.style.display = hideRow ? "none" : "";
            }}
          }}

          applyFilter();
          if (checkboxNonCnv)   checkboxNonCnv.addEventListener("change", applyFilter);
          if (checkboxPonSig)   checkboxPonSig.addEventListener("change", applyFilter);
          if (checkboxCancer)   checkboxCancer.addEventListener("change", applyFilter);
          if (geneInput)        geneInput.addEventListener("input", applyFilter);
          if (minTargetsInput)  minTargetsInput.addEventListener("input", applyFilter);
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
    required=True,
    help="cnvkit scatter PDF file.",
)
@click.option(
    "--cnvkit-diagram",
    type=click.Path(exists=True),
    required=True,
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

    pdf_first_page_to_png(cnvkit_scatter, scatter_png_path)
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
    # pon can be None (no PON file); build_gene_table handles that
    df_genes = build_gene_segment_table(
        cnr_path=cnr,
        cns_path=cns,
        cancer_genes=cancer_gene_set,
        refgene_path=refgene,
        transcript_selection="longest_tx",
        pon_path=pon,
        loh_path=loh_genes,
        cytoband_path=cytoband,
        sex=sex,
    )

    # ---------------
    # Read PureCN purity and ploidy estimation (optional)
    # ---------------
    purity_df = None
    if purity_csv:
        purity_df = pd.read_csv(purity_csv)

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
        gene_seg_df=df_genes,
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
        gene_seg_df=df_genes,
        outdir=chr_plots_dir,
        case_id=plot_case_id,
        pon_path=pon,
        include_y=True,
        neutral_target_factor=neutral_target_factor,
        highlight_only_cancer=highlight_only_cancer,
        focus_genes=["RB1", "DLEU1"],
        focus_padding_bp=100_000,
    )

    # ----------------------------
    # HTML report
    # ----------------------------
    out_html = str(out_prefix)
    csv_to_html_table(
        df=df_genes,
        out_html=out_html,
        scatter_png=str(scatter_png_path),
        diagram_png=str(diagram_png_path),
        chr_plots_dir=chr_plots_dir,
        summary_df=purity_df,
    )

    click.echo(f"[CNV QC] Finished report for {case_id}: {out_html}")


if __name__ == "__main__":
    sys.exit(main())
