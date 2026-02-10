(function () {
  const geneTable = document.getElementById("report-table");
  const chunkTable = document.getElementById("chunk-table");

  // Gene-level controls
  const gHideNonCnv = document.getElementById("hide-non-cnv");
  const gOnlyCancer = document.getElementById("only-cancer-genes");
  const gGeneInput = document.getElementById("gene-filter");
  const gMinTargets = document.getElementById("min-targets");

  // Chunk-level controls
  const cHideNonCnv = document.getElementById("hide-non-cnv-chunk");
  const cIncludePonCnv = document.getElementById("show-pon-cnv-chunk");
  const cOnlyCancer = document.getElementById("only-cancer-genes-chunk");
  const cOnlySplit = document.getElementById("only-split-genes-chunk");
  const cGeneInput = document.getElementById("gene-filter-chunk");
  const cMinTargets = document.getElementById("min-targets-chunk");

  // ---------------------------------------------------------------------------
  // Column group toggles (Option 2)
  // Hidden by default; shown when group checkbox is checked.
  // ---------------------------------------------------------------------------

  // Edit these lists/prefixes to match your dataframes.
  // Keep the PON group as prefix-based so new pon_* columns auto-follow.
  const COL_GROUP_QC = [
    "depth_mean",
    "mean_weight",
    // add more QC-ish columns here if desired
  ];

  // PON group: any columns starting with these prefixes are treated as "PON details"
  const COL_GROUP_PON_PREFIXES = ["pon_"];

  // Optional extra PON cols not following prefix convention (leave empty if none)
  const COL_GROUP_PON_EXPLICIT = [
    // e.g. "ponChunkZ"
  ];

  // CNVkit group (optional)
  const COL_GROUP_CNVKIT_EXTRA = [
      "cnvkit_seg_log2", "cnvkit_seg_cn", "cnvkit_seg_cn1", "cnvkit_seg_cn2"
  ];

    // PureCN group (optional)
  const COL_GROUP_PURECN_EXTRA = [
      "purecn_C", "purecn_M", "purecn_M_flagged",
  ];

  // Default visibility per table
  const DEFAULT_GROUP_STATE = {
    "report-table": {
      qc: false,
      pon: false,
      cnvkit: false,
      purecn: false,
    },
    "chunk-table": {
      qc: false,
      pon: true,     // <-- PON columns shown by default here
      cnvkit: false,
      purecn: false,
    },
  };

  function getColIndexForTableId(tableId) {
    if (tableId === "report-table") {
      return (typeof colIndexGene !== "undefined" ? colIndexGene : {});
    }
    if (tableId === "chunk-table") {
      return (typeof colIndexChunk !== "undefined" ? colIndexChunk : {});
    }
    return {};
  }

  function getTableById(tableId) {
    return document.getElementById(tableId);
  }

  function setColumnVisibility(table, colIndex, colName, visible) {
    if (!table) return;

    const idx = (colIndex && (colName in colIndex)) ? colIndex[colName] : -1;
    if (idx === -1) return;

    const displayValue = visible ? "" : "none";

    // Header cells
    if (table.tHead && table.tHead.rows) {
      for (const hr of Array.from(table.tHead.rows)) {
        if (hr.cells && hr.cells[idx]) {
          hr.cells[idx].style.display = displayValue;
        }
      }
    }

    // Body cells
    if (table.tBodies) {
      for (const tb of Array.from(table.tBodies)) {
        for (const row of Array.from(tb.rows)) {
          if (row.cells && row.cells[idx]) {
            row.cells[idx].style.display = displayValue;
          }
        }
      }
    }
  }

  function columnsMatchingPrefixes(colIndex, prefixes) {
    const cols = Object.keys(colIndex || {});
    const out = [];
    for (const c of cols) {
      for (const p of prefixes) {
        if (c.startsWith(p)) {
          out.push(c);
          break;
        }
      }
    }
    return out;
  }

  function unique(arr) {
    return Array.from(new Set(arr));
  }

  function getGroupColumns(colIndex, groupName) {
    if (!colIndex) return [];

    if (groupName === "qc") {
      return COL_GROUP_QC.filter((c) => c in colIndex);
    }

    if (groupName === "pon") {
      const prefixCols = columnsMatchingPrefixes(colIndex, COL_GROUP_PON_PREFIXES);
      const explicitCols = COL_GROUP_PON_EXPLICIT.filter((c) => c in colIndex);
      return unique([...prefixCols, ...explicitCols]);
    }

    if (groupName === "cnvkit") {
      return COL_GROUP_CNVKIT_EXTRA.filter((c) => c in colIndex);
    }

    if (groupName === "purecn") {
      return COL_GROUP_PURECN_EXTRA.filter((c) => c in colIndex);
    }

    return [];
  }

  function applyColumnGroupsForTable(tableId, toggles) {
    const table = getTableById(tableId);
    if (!table) return;

    const colIndex = getColIndexForTableId(tableId);

    const qcCols = getGroupColumns(colIndex, "qc");
    const ponCols = getGroupColumns(colIndex, "pon");
    const cnvkitCols = getGroupColumns(colIndex, "cnvkit");
    const purecnCols = getGroupColumns(colIndex, "purecn");

    // Hidden by default: show only if the toggle is checked
    for (const c of qcCols) setColumnVisibility(table, colIndex, c, !!toggles.qc);
    for (const c of ponCols) setColumnVisibility(table, colIndex, c, !!toggles.pon);
    for (const c of cnvkitCols) setColumnVisibility(table, colIndex, c, !!toggles.cnvkit);
    for (const c of purecnCols) setColumnVisibility(table, colIndex, c, !!toggles.purecn);
  }

  function hookColumnGroupToggles() {
    // Gene toggles
    const gQc = document.getElementById("cols-qc-gene");
    const gPon = document.getElementById("cols-pon-gene");
    const gCnvkit = document.getElementById("cols-cnvkit-extra-gene");
    const gPurecn = document.getElementById("cols-purecn-extra-gene");

    // Chunk toggles
    const cQc = document.getElementById("cols-qc-chunk");
    const cPon = document.getElementById("cols-pon-chunk");
    const cCnvkit = document.getElementById("cols-cnvkit-extra-chunk");
    const cPurecn = document.getElementById("cols-purecn-extra-chunk");

    function applyAll() {
      applyColumnGroupsForTable("report-table", {
        qc: gQc && gQc.checked,
        pon: gPon && gPon.checked,
        cnvkit: gCnvkit && gCnvkit.checked,
        purecn: gPurecn && gPurecn.checked,
      });

      applyColumnGroupsForTable("chunk-table", {
        qc: cQc && cQc.checked,
        pon: cPon && cPon.checked,
        cnvkit: cCnvkit && cCnvkit.checked,
        purecn: cPurecn && cPurecn.checked,
      });
    }

    // --- Apply default states
    const geneDefaults = DEFAULT_GROUP_STATE["report-table"];
    if (gQc && geneDefaults) gQc.checked = geneDefaults.qc;
    if (gPon && geneDefaults) gPon.checked = geneDefaults.pon;
    if (gCnvkit && geneDefaults) gCnvkit.checked = geneDefaults.cnvkit;
    if (gPurecn && geneDefaults) gPurecn.checked = geneDefaults.purecn;

    const chunkDefaults = DEFAULT_GROUP_STATE["chunk-table"];
    if (cQc && chunkDefaults) cQc.checked = chunkDefaults.qc;
    if (cPon && chunkDefaults) cPon.checked = chunkDefaults.pon;
    if (cCnvkit && chunkDefaults) cCnvkit.checked = chunkDefaults.cnvkit;
    if (cPurecn && chunkDefaults) cPurecn.checked = chunkDefaults.purecn;

    // Apply initial state
    applyAll();

    // Events
    if (gQc) gQc.addEventListener("change", applyAll);
    if (gPon) gPon.addEventListener("change", applyAll);
    if (gCnvkit) gCnvkit.addEventListener("change", applyAll);
    if (gPurecn) gPurecn.addEventListener("change", applyAll);

    if (cQc) cQc.addEventListener("change", applyAll);
    if (cPon) cPon.addEventListener("change", applyAll);
    if (cCnvkit) cCnvkit.addEventListener("change", applyAll);
    if (cPurecn) cPurecn.addEventListener("change", applyAll);
  }

  // ---------------------------------------------------------------------------
  // Existing row-filter logic
  // ---------------------------------------------------------------------------

  function normalizeCellText(td) {
    return (td && td.textContent ? td.textContent : "").trim().toLowerCase();
  }

  function isCnvRow(row, colIndex) {
  const lohColIdx = (colIndex["purecn_loh_flag"] ?? colIndex["loh_flag"] ?? -1);
  const cnvkitColIdx = (colIndex["cnvkit_cnv_call"] ?? -1);
  const purecnColIdx = (colIndex["purecn_cnv_call"] ?? -1);

  // LOH: keep if true/1/yes
  if (lohColIdx !== -1) {
    const lohVal = normalizeCellText(row.cells[lohColIdx]);
    if (lohVal === "true" || lohVal === "1" || lohVal === "yes") return true;
  }

  // CNV calls: keep if present and not NEUTRAL/empty/na
  function isNonNeutralCall(cell) {
    const v = normalizeCellText(cell);
    if (!v) return false;
    if (v === "neutral") return false;
    if (v === "na" || v === "nan" || v === ".") return false;
    return true;
  }

  if (cnvkitColIdx !== -1 && isNonNeutralCall(row.cells[cnvkitColIdx])) return true;
  if (purecnColIdx !== -1 && isNonNeutralCall(row.cells[purecnColIdx])) return true;

  return false;
}

  function isCancerGeneRow(row, colIndex) {
    const idx = colIndex["is_cancer_gene"] ?? -1;
    if (idx === -1) return false;
    const val = normalizeCellText(row.cells[idx]);
    return val === "true" || val === "1" || val === "yes";
  }

  function isSplitGeneRow(row, colIndex) {
    const idx = colIndex["is_gene_split"] ?? -1;
    if (idx === -1) return false;
    const val = normalizeCellText(row.cells[idx]);
    return val === "true" || val === "1" || val === "yes";
  }

  function parseGeneList(value) {
    return (value || "")
      .split(",")
      .map((g) => g.trim().toLowerCase())
      .filter((g) => g.length > 0);
  }

  function filterTable(table, colIndex, cfg) {
    if (!table) return;
    const body = table.tBodies[0];
    if (!body) return;

    const geneColIdx = colIndex["gene.symbol"] ?? -1;
    const nTargetsColIdx = colIndex["n.targets"] ?? -1;
    const ponCnvColIdx = (colIndex["pon_cnv_call"] ?? -1);

    const rows = Array.from(body.rows);

    for (const row of rows) {
      let hideRow = false;

      if (!hideRow && cfg.hideNonCnv) {
        const cnv = isCnvRow(row, colIndex);

        let isPonCnv = false;
        if (ponCnvColIdx !== -1) {
          const ponVal = normalizeCellText(row.cells[ponCnvColIdx]);
          if (ponVal === "amplification" || ponVal === "deletion") isPonCnv = true;
        }

        const treatedAsCnv = cnv || (cfg.includePonCnv && isPonCnv);
        if (!treatedAsCnv) hideRow = true;
      }

      if (!hideRow && cfg.onlyCancer) {
        if (!isCancerGeneRow(row, colIndex)) hideRow = true;
      }

      if (!hideRow && cfg.onlySplit) {
        if (!isSplitGeneRow(row, colIndex)) hideRow = true;
      }

      if (!hideRow && nTargetsColIdx !== -1 && cfg.minTargets > 0) {
        const nText = normalizeCellText(row.cells[nTargetsColIdx]);
        if (nText.length > 0) {
          const nVal = parseFloat(nText);
          if (!Number.isNaN(nVal) && nVal < cfg.minTargets) hideRow = true;
        }
      }

      if (!hideRow && geneColIdx !== -1 && cfg.geneList.length > 0) {
        const geneText = normalizeCellText(row.cells[geneColIdx]);
        const matches = cfg.geneList.some((g) => geneText === g);
        if (!matches) hideRow = true;
      }

      row.style.display = hideRow ? "none" : "";
    }
  }

  function filterGeneTable() {
    if (!geneTable) return;
    const cfg = {
      hideNonCnv: !!(gHideNonCnv && gHideNonCnv.checked),
      includePonCnv: false,
      onlyCancer: !!(gOnlyCancer && gOnlyCancer.checked),
      onlySplit: false,
      geneList: parseGeneList(gGeneInput && gGeneInput.value),
      minTargets: gMinTargets ? parseInt(gMinTargets.value, 10) || 0 : 0,
    };
    filterTable(geneTable, (typeof colIndexGene !== "undefined" ? colIndexGene : {}), cfg);
  }

  function filterChunkTable() {
    if (!chunkTable) return;
    const cfg = {
      hideNonCnv: !!(cHideNonCnv && cHideNonCnv.checked),
      includePonCnv: !!(cIncludePonCnv && cIncludePonCnv.checked),
      onlyCancer: !!(cOnlyCancer && cOnlyCancer.checked),
      onlySplit: !!(cOnlySplit && cOnlySplit.checked),
      geneList: parseGeneList(cGeneInput && cGeneInput.value),
      minTargets: cMinTargets ? parseInt(cMinTargets.value, 10) || 0 : 0,
    };
    filterTable(chunkTable, (typeof colIndexChunk !== "undefined" ? colIndexChunk : {}), cfg);
  }

  function applyFilter() {
    filterGeneTable();
    filterChunkTable();
  }

  if (gHideNonCnv) gHideNonCnv.addEventListener("change", applyFilter);
  if (gOnlyCancer) gOnlyCancer.addEventListener("change", applyFilter);
  if (gGeneInput) gGeneInput.addEventListener("input", applyFilter);
  if (gMinTargets) gMinTargets.addEventListener("input", applyFilter);

  if (cHideNonCnv) cHideNonCnv.addEventListener("change", applyFilter);
  if (cIncludePonCnv) cIncludePonCnv.addEventListener("change", applyFilter);
  if (cOnlyCancer) cOnlyCancer.addEventListener("change", applyFilter);
  if (cOnlySplit) cOnlySplit.addEventListener("change", applyFilter);
  if (cGeneInput) cGeneInput.addEventListener("input", applyFilter);
  if (cMinTargets) cMinTargets.addEventListener("input", applyFilter);

  // ---------------------------------------------------------------------------
  // Table show/hide toggles
  // ---------------------------------------------------------------------------

  const geneToggleBtn = document.getElementById("toggle-gene-table");
  const geneContainer = document.getElementById("gene-table-container");
  if (geneToggleBtn && geneContainer) {
    let geneVisible = true;
    geneToggleBtn.addEventListener("click", () => {
      geneVisible = !geneVisible;
      geneContainer.style.display = geneVisible ? "block" : "none";
      geneToggleBtn.textContent = geneVisible ? "Hide gene-level table" : "Show gene-level table";
    });
  }

  const chunkToggleBtn = document.getElementById("toggle-chunk-table");
  const chunkContainer = document.getElementById("chunk-table-container");
  if (chunkToggleBtn && chunkContainer) {
    let chunkVisible = false;
    chunkToggleBtn.addEventListener("click", () => {
      chunkVisible = !chunkVisible;
      chunkContainer.style.display = chunkVisible ? "block" : "none";
      chunkToggleBtn.textContent = chunkVisible ? "Hide chunk-level table" : "Show chunk-level table";
    });
  }

  // ---------------------------------------------------------------------------
  // Plot modal viewer
  // ---------------------------------------------------------------------------

  const modal = document.getElementById("plot-modal");
  const modalImg = document.getElementById("plot-modal-image");
  const modalClose = document.getElementById("plot-modal-close");

  function openModal(src) {
    if (!modal || !modalImg) return;
    modalImg.src = src;
    modal.style.display = "block";
  }

  function closeModal() {
    if (!modal) return;
    modal.style.display = "none";
    if (modalImg) modalImg.src = "";
  }

  function hookClickablePlots() {
    const imgs = document.querySelectorAll(".clickable-plot");
    imgs.forEach((img) => {
      img.addEventListener("click", () => openModal(img.src));
    });
  }

  hookClickablePlots();

  if (modalClose) modalClose.addEventListener("click", closeModal);

  if (modal) {
    modal.addEventListener("click", (event) => {
      if (event.target === modal) closeModal();
    });
  }

  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape") closeModal();
  });

  // ---------------------------------------------------------------------------
  // Init order
  // ---------------------------------------------------------------------------

  // 1) Hook and apply column group toggles (unchecked => hidden)
  hookColumnGroupToggles();

  // 2) Apply row filters once
  applyFilter();
})();
