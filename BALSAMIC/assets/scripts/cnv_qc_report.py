import sys
import click
import base64
import json
import pandas as pd
from pathlib import Path
from cnv_report_utils import plot_chromosomes, build_gene_segment_table, pdf_first_page_to_png, load_cancer_gene_set
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
    # Determine which chromosomes have CNV genes
    # CNV = (C != 2) OR (loh == TRUE) OR (type annotated and not NA/NAN/./empty)
    # ------------------------------------------------------------------
    df_chr = df.copy()
    if "chr" in df_chr.columns:
        chr_col = "chr"
    elif "chromosome" in df_chr.columns:
        chr_col = "chromosome"
    else:
        chr_col = None

    cnv_chr_set: set[str] = set()

    if chr_col is not None and {"C", "loh", "type"}.issubset(df_chr.columns):
        if "gene.symbol" in df_chr.columns:
            df_chr = df_chr[~df_chr["gene.symbol"].isin(["Antitarget", "-"])]

        # Normalize
        df_chr["loh_str"] = df_chr["loh"].astype(str).str.upper()
        df_chr["type_str"] = df_chr["type"].fillna("NA").astype(str).str.upper()

        # Missing / placeholder types
        bad_type = {"NA", "NAN", ".", "", "<NA>", "NONE"}
        is_type_real = df_chr["type_str"].notna() & ~df_chr["type_str"].isin(bad_type)

        cnv_mask = (df_chr["loh_str"] == "TRUE") | is_type_real

        cnv_chr_set = set(df_chr.loc[cnv_mask, chr_col].astype(str).tolist())

    # -------------------------------------------------------------
    # Per-chromosome QC PNGs (sorted, split by CNV / non-CNV)
    # -------------------------------------------------------------
    with_cnv_blocks: list[str] = []
    no_cnv_blocks: list[str] = []

    if chr_plots_dir is not None:
        qc_dir = Path(chr_plots_dir)
        if qc_dir.is_dir():
            png_files = list(qc_dir.glob("cnv_chr*segments.png"))

            def chr_sort_key(path: Path):
                stem = path.stem
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
                chr_label = stem
                if "chr" in stem:
                    chr_label = stem.split("chr", 1)[1]
                if "_" in chr_label:
                    chr_label = chr_label.split("_", 1)[0]

                block = f"""
                <div style="border:1px solid #ccc; padding:10px; margin-top:20px; width:100%;">
                  <h3 style="margin-top:0;">Chr {chr_label}</h3>
                  <img
                    src="{img_data_uri}"
                    alt="QC plot {stem}"
                    style="max-width:100%; height:auto; display:block;"
                  />
                </div>
                """

                if chr_col is not None and chr_label in cnv_chr_set:
                    with_cnv_blocks.append(block)
                else:
                    no_cnv_blocks.append(block)

    qc_plots_html = ""
    if with_cnv_blocks or no_cnv_blocks:
        qc_plots_html += """
        <h2>Per-chromosome CNV QC plots</h2>
        <p class="muted">
          These plots show log2 coverage, PON spread (if available), BAF and PureCN segments/genes
          per chromosome. CNV genes are highlighted in colour.
        </p>
        """
        if with_cnv_blocks:
            qc_plots_html += "<h3>Chromosomes with CNV genes</h3>\n" + "\n".join(
                with_cnv_blocks
            )
        if no_cnv_blocks:
            qc_plots_html += """
            <details style="margin-top:15px;">
              <summary style="cursor:pointer; font-weight:500;">
                Show chromosomes without CNV genes
              </summary>
            """
            qc_plots_html += "\n".join(no_cnv_blocks)
            qc_plots_html += "</details>"

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
      </style>
    </head>
    <body>
      <h1>CNV Report</h1>

      <h2>Sample summary</h2>
      {summary_html}

      <h2>CNVkit scatter</h2>
      <div style="border:1px solid #ccc; padding:10px; margin-top:20px; width:100%;">
        <img
          src="{scatter_png_data_uri}"
          alt="CNVkit scatter PNG plot"
          style="max-width:100%; height:auto; display:block;"
        />
      </div>

      <h2>CNVkit diagram</h2>
      <div style="border:1px solid #ccc; padding:10px; margin-top:20px; width:100%;">
        <img
          src="{diagram_png_data_uri}"
          alt="CNVkit diagram PNG plot"
          style="max-width:100%; height:auto; display:block;"
        />
      </div>

      {qc_plots_html}

      <h1>LOH / CNV Genes</h1>

      <div class="controls">
        <div class="control-group">
          <label>
            <input type="checkbox" id="hide-non-cnv">
            Show only LOH / CNV genes (loh = TRUE or type set)
          </label>
          <label style="margin-left:12px;">
            <input type="checkbox" id="show-beyond-pon-noise">
            Also include genes with pon_spread_flag = "beyond_pon_noise"
          </label>
        </div>
        <div class="control-group">
          <label>
            Min number.targets:
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
          const table           = document.getElementById("report-table");
          const checkboxNonCnv  = document.getElementById("hide-non-cnv");
          const checkboxBeyond  = document.getElementById("show-beyond-pon-noise");
          const geneInput       = document.getElementById("gene-filter");
          const minTargetsInput = document.getElementById("min-targets");

          if (!table) return;

          const colIndex = {col_idx_json};

          const geneColIdx      = colIndex["gene.symbol"] ?? -1;
          const lohColIdx       = colIndex["loh"] ?? -1;
          const typeColIdx      = colIndex["type"] ?? -1;
          const nTargetsColIdx  = colIndex["number.targets"] ?? -1;
          const ponSpreadColIdx = colIndex["pon_spread_flag"] ?? -1;

          function normalizeCellText(td) {{
            return (td && td.textContent ? td.textContent : "")
              .trim()
              .toLowerCase();
          }}

          function isCnvRow(row) {{
            let isCnv = false;
            let evaluated = false;

            if (lohColIdx !== -1) {{
              const lohCell = row.cells[lohColIdx];
              const lohVal = normalizeCellText(lohCell);
              if (lohVal === "true") {{
                evaluated = true;
                isCnv = true;
              }}
            }}

            if (!isCnv && typeColIdx !== -1) {{
              const typeCell = row.cells[typeColIdx];
              const typeVal = normalizeCellText(typeCell);
              const bad = ["na", "nan", ".", "", "<na>", "none"];
              if (typeVal && !bad.includes(typeVal)) {{
                evaluated = true;
                isCnv = true;
              }}
            }}

            return evaluated ? isCnv : false;
          }}

          function applyFilter() {{
            const hideNonCnv    = checkboxNonCnv && checkboxNonCnv.checked;
            const showBeyondPon = checkboxBeyond && checkboxBeyond.checked;

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

              // CNV / non-CNV filter with beyond_pon_noise override
              if (!hideRow && hideNonCnv) {{
                const isCnv = isCnvRow(row);

                let isBeyond = false;
                if (ponSpreadColIdx !== -1) {{
                  const ponCell = row.cells[ponSpreadColIdx];
                  const ponVal = normalizeCellText(ponCell);
                  isBeyond = (ponVal === "beyond_pon_noise");
                }}

                const treatedAsCnv = isCnv || (showBeyondPon && isBeyond);
                if (!treatedAsCnv) {{
                  hideRow = true;
                }}
              }}

              // min number.targets filter
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

              // gene.symbol filter
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
          if (checkboxBeyond)   checkboxBeyond.addEventListener("change", applyFilter);
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
    "--loh-regions",
    type=click.Path(exists=True),
    required=True,  # PureCN regions CSV (can be empty but must exist)
    help="PureCN LOH regions CSV (…_loh.csv).",
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
def main(
    loh_genes,
    loh_regions,
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

    # pon may be None; plot_chromosomes is written to handle that
    plot_chromosomes(
        pon_path=pon,
        cnr_path=cnr,
        cns_path=cns,
        vcf_path=vcf,
        genes_path=loh_genes,
        segs_path=loh_regions,
        outdir=chr_plots_dir,
        case_id=plot_case_id,
        cancer_genes=cancer_gene_set,
        neutral_target_factor=neutral_target_factor,
        is_exome=is_exome,
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
