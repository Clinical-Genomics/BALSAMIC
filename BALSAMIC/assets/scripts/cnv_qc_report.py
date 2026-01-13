import sys
import click
import base64
import json
import pandas as pd
from pathlib import Path
from cnv_report_utils import plot_chromosomes, build_gene_table, pdf_first_page_to_png


def csv_to_html_table(
    df: pd.DataFrame,
    out_html: str,
    scatter_png: str,
    diagram_png: str,
    chr_plots_dir: str | Path | None = None,
) -> None:
    out_path = Path(out_html)

    # --- PNG → Base64 (original PureCN plot) ---
    png_bytes = Path(scatter_png).read_bytes()
    png_base64 = base64.b64encode(png_bytes).decode("ascii")
    scatter_png_data_uri = f"data:image/png;base64,{png_base64}"

    # --- PNG → Base64 (original PureCN plot) ---
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

        df_chr["loh_str"] = df_chr["loh"].astype(str).str.upper()
        df_chr["type_str"] = df_chr["type"].astype(str).str.upper()

        is_type_real = ~df_chr["type_str"].isin(["NA", "NAN", ".", ""])

        cnv_mask = (
            ((df_chr["C"].notna()) & (df_chr["C"] != 2))
            | (df_chr["loh_str"] == "TRUE")
            | is_type_real
        )

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
          These plots show log2 coverage, PON spread (noise band), BAF and PureCN segments/genes
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

    # ---------------- Tables ----------------
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
    }}
    .control-group {{
      margin-bottom: 6px;
    }}
    table.dataframe {{
      border-collapse: collapse;
      width: 100%;
      font-size: 14px;
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
  </style>
</head>
<body>
  <h1>CNV Report</h1>

  <h2>PureCN plot (PNG embedded)</h2>
  <div style="border:1px solid #ccc; padding:10px; margin-top:20px; width:100%;">
    <img
      src="{scatter_png_data_uri}"
      alt="PureCN PNG plot"
      style="max-width:100%; height:auto; display:block;"
    />
  </div>

  <h2>PureCN plot (PNG embedded)</h2>
  <div style="border:1px solid #ccc; padding:10px; margin-top:20px; width:100%;">
    <img
      src="{diagram_png_data_uri}"
      alt="PureCN PNG plot"
      style="max-width:100%; height:auto; display:block;"
    />
  </div>

  {qc_plots_html}

  <h1>LOH Genes</h1>

  <div class="controls">
    <div class="control-group">
      <label>
        <input type="checkbox" id="hide-flagged" checked>
        Hide rows where <code>M.flagged</code> is True
      </label>
    </div>
    <div class="control-group">
      <label>
        <input type="checkbox" id="hide-non-cnv" checked>
        Hide genes not affected by CNV
        <span class="muted">
          (C &ne; 2, or loh = TRUE, or type annotated)
        </span>
      </label>
    </div>
    <div class="control-group">
      <label>
        Min <code>number.targets</code>:
        <input type="number" id="min-targets" value="2" min="0" step="1">
      </label>
    </div>
    <div class="control-group">
      <label>
        Filter <code>gene.symbol</code> (comma-separated):
        <input type="text" id="gene-filter" placeholder="e.g. FLT3, NPM1, DNMT3A">
      </label>
    </div>
    <div class="control-group">
      <button id="copy-table" type="button">
        Copy visible table to clipboard
      </button>
    </div>
    <div class="muted">
      Use the checkboxes to hide flagged calls and non-CNV genes,
      limit by min number.targets, and type one or more genes separated by commas
      to focus the list. Use the button to copy the currently visible table
      to the clipboard and paste into Excel.
    </div>
  </div>

  {table_html}

  <script>
    (function() {{
      const table           = document.getElementById("report-table");
      const checkboxFlagged = document.getElementById("hide-flagged");
      const checkboxNonCnv  = document.getElementById("hide-non-cnv");
      const geneInput       = document.getElementById("gene-filter");
      const minTargetsInput = document.getElementById("min-targets");
      const copyBtn         = document.getElementById("copy-table");

      if (!table) return;

      const colIndex = {col_idx_json};

      const flaggedColIdx  = colIndex["M.flagged"] ?? -1;
      const geneColIdx     = colIndex["gene.symbol"] ?? -1;
      const cColIdx        = colIndex["C"] ?? -1;
      const lohColIdx      = colIndex["loh"] ?? -1;
      const typeColIdx     = colIndex["type"] ?? -1;
      const nTargetsColIdx = colIndex["number.targets"] ?? -1;

      function normalizeCellText(td) {{
        return (td && td.textContent ? td.textContent : "")
          .trim()
          .toLowerCase();
      }}

      function isCnvRow(row) {{
        // CNV = C != 2 OR loh == TRUE OR type annotated (not NA/NAN/./empty)
        let isCnv = false;
        let evaluated = false;

        // C != 2
        if (cColIdx !== -1) {{
          const cCell = row.cells[cColIdx];
          const cText = normalizeCellText(cCell);
          if (cText.length > 0) {{
            const cVal = parseFloat(cText);
            if (!Number.isNaN(cVal)) {{
              evaluated = true;
              if (cVal !== 2) {{
                isCnv = true;
              }}
            }}
          }}
        }}

        // loh == TRUE
        if (!isCnv && lohColIdx !== -1) {{
          const lohCell = row.cells[lohColIdx];
          const lohVal = normalizeCellText(lohCell);
          if (lohVal === "true") {{
            evaluated = true;
            isCnv = true;
          }}
        }}

        // type annotated
        if (!isCnv && typeColIdx !== -1) {{
          const typeCell = row.cells[typeColIdx];
          const typeVal = normalizeCellText(typeCell);
          const bad = ["na", "nan", ".", ""];
          if (typeVal && !bad.includes(typeVal)) {{
            evaluated = true;
            isCnv = true;
          }}
        }}

        return evaluated ? isCnv : false;
      }}

      function applyFilter() {{
        const hideTrue   = checkboxFlagged && checkboxFlagged.checked;
        const hideNonCnv = checkboxNonCnv && checkboxNonCnv.checked;

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

          // M.flagged filter
          if (flaggedColIdx !== -1 && checkboxFlagged) {{
            const flaggedCell = row.cells[flaggedColIdx];
            const flaggedVal = normalizeCellText(flaggedCell);
            const isTrue = (flaggedVal === "true");
            if (hideTrue && isTrue) {{
              hideRow = true;
            }}
          }}

          // CNV filter
          if (!hideRow && hideNonCnv && checkboxNonCnv) {{
            if (!isCnvRow(row)) {{
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
      if (checkboxFlagged)  checkboxFlagged.addEventListener("change", applyFilter);
      if (checkboxNonCnv)   checkboxNonCnv.addEventListener("change", applyFilter);
      if (geneInput)        geneInput.addEventListener("input", applyFilter);
      if (minTargetsInput)  minTargetsInput.addEventListener("input", applyFilter);

      // --------- Copy visible table to clipboard ---------
      function copyVisibleTableToClipboard() {{
        const body = table.tBodies[0];
        if (!body) return;

        const lines = [];

        // Header
        let headerCells = [];
        if (table.tHead && table.tHead.rows.length > 0) {{
          headerCells = Array.from(table.tHead.rows[0].cells);
        }} else if (table.rows.length > 0) {{
          headerCells = Array.from(table.rows[0].cells);
        }}
        if (headerCells.length > 0) {{
          const headerLine = headerCells
            .map(th => th.textContent.trim())
            .join("\\t");
          lines.push(headerLine);
        }}

        // Visible rows only
        const rows = Array.from(body.rows);
        for (const row of rows) {{
          if (row.style.display === "none") continue;
          const cells = Array.from(row.cells);
          const line = cells
            .map(td => td.textContent.trim())
            .join("\\t");
          lines.push(line);
        }}

        const text = lines.join("\\n");

        if (navigator.clipboard && navigator.clipboard.writeText) {{
          navigator.clipboard.writeText(text).then(
            () => console.log("Table copied to clipboard"),
            err => console.warn("Clipboard copy failed:", err)
          );
        }} else {{
          const textarea = document.createElement("textarea");
          textarea.value = text;
          textarea.style.position = "fixed";
          textarea.style.left = "-9999px";
          document.body.appendChild(textarea);
          textarea.select();
          try {{
            document.execCommand("copy");
            console.log("Table copied to clipboard (fallback).");
          }} catch (e) {{
            console.warn("Clipboard copy failed (fallback):", e);
          }}
          document.body.removeChild(textarea);
        }}
      }}

      if (copyBtn) {{
        copyBtn.addEventListener("click", copyVisibleTableToClipboard);
      }}
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
    required=True,
    help="PureCN LOH genes CSV (…_genes.csv).",
)
@click.option(
    "--loh-regions",
    type=click.Path(exists=True),
    required=True,
    help="PureCN LOH regions CSV (…_loh.csv).",
)
@click.option(
    "--cnr", type=click.Path(exists=True), required=True, help="CNVkit tumor .cnr file."
)
@click.option(
    "--pon", type=click.Path(exists=True), required=True, help="CNVkit PON .cnn file."
)
@click.option(
    "--vcf",
    type=click.Path(exists=True),
    required=True,
    help="VCF with germline variants for BAF.",
)
@click.option(
    "--refgene", type=click.Path(exists=True), required=True, help="refgene.flat file."
)
@click.option(
    "--cytoband",
    type=click.Path(exists=True),
    required=True,
    help="cytoBand file for genome build.",
)
@click.option(
    "--case-id", type=str, required=True, help="Case ID for labels / outputs."
)
@click.option(
    "--cust-case-id", type=str, required=False, help="Cust Case ID for labels / outputs."
)
@click.option(
    "--purecn-scatter",
    type=click.Path(exists=True),
    required=True,
    help="purecn scatter PDF file.",
)
@click.option(
    "--purecn-diagram",
    type=click.Path(exists=True),
    required=True,
    help="purecn diagram PDF file.",
)
@click.option(
    "--out-prefix",
    type=click.Path(),
    required=True,
    help="Output prefix (e.g. /path/to/sample_cnv_qc).",
)
def main(
    loh_genes,
    loh_regions,
    cnr,
    pon,
    vcf,
    refgene,
    cytoband,
    case_id,
    cust_case_id,,
    purecn_scatter,
    purecn_diagram,
    out_prefix,
):
    """
    Build CNV QC plots + HTML report from PureCN + CNVkit outputs.
    """

    out_prefix = Path(out_prefix)
    outdir = out_prefix.parent
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------------
    # Convert PDF to PNG
    # ---------------
    pdf_first_page_to_png(purecn_diagram, f"{outdir}/purecn_diagram_{case_id}.png")
    pdf_first_page_to_png(purecn_scatter, f"{outdir}/purecn_scatter_{case_id}.png")

    # ----------------------------
    # Generate per-chromosome PNG plots in outdir
    # ----------------------------
    df_genes = build_gene_table(loh_regions, loh_genes, cnr, pon, refgene, cytoband)

    # ----------------------------
    # Generate per-chromosome PNG plots in outdir
    # ----------------------------
    chr_plots_dir = outdir / f"{case_id}_chr_plots"
    chr_plots_dir.mkdir(exist_ok=True, parents=True)
    if cust_case_id:
        plot_case_id = cust_case_id
    else:
        plot_case_id = case_id
    plot_chromosomes(pon, cnr, vcf, loh_genes, loh_regions, chr_plots_dir, plot_case_id)

    # ----------------------------
    # HTML report
    # ----------------------------
    out_html = str(out_prefix)
    csv_to_html_table(
        df=df_genes,
        out_html=out_html,
        scatter_png=f"{outdir}/purecn_scatter_{case_id}.png",
        diagram_png=f"{outdir}/purecn_diagram_{case_id}.png",
        chr_plots_dir=chr_plots_dir,
    )

    click.echo(f"[CNV QC] Finished report for {case_id}: {out_html}")


if __name__ == "__main__":
    # Use Click for argument parsing
    sys.exit(main())
