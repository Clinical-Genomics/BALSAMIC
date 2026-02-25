(function () {
  // ---------------------------------------------------------------------------
  // Small DOM / utility helpers
  // ---------------------------------------------------------------------------

  const $ = (id) => document.getElementById(id);

  const textOf = (cell) => (cell?.textContent ?? "").trim().toLowerCase();

  const isTruthyText = (v) => v === "true" || v === "1" || v === "yes";

  const parseIntSafe = (v) => {
    const n = parseInt(v, 10);
    return Number.isFinite(n) ? n : 0;
  };

  const unique = (arr) => Array.from(new Set(arr));

  // Prefer globalThis to avoid typeof checks + Sonar negated condition warning
  function getColIndex(tableId) {
    const map = {
      "report-table": globalThis.colIndexGene,
      "chunk-table": globalThis.colIndexChunk,
    };
    return map[tableId] ?? {};
  }

  function getTable(tableId) {
    return $(tableId);
  }

  function bindEvents(pairs, handler) {
    for (const [el, ev] of pairs) {
      if (el) el.addEventListener(ev, handler);
    }
  }

  const isPresentText = (v) => {
  v = (v ?? "").trim().toLowerCase();
  if (!v) return false;
  return !["na", "nan", "<na>", ".", "none", "null"].includes(v);
  };

  // ---------------------------------------------------------------------------
  // Tables + controls
  // ---------------------------------------------------------------------------

  const geneTable = getTable("report-table");
  const chunkTable = getTable("chunk-table");

  const controls = {
    gene: {
      hideNonCnv: $("hide-non-cnv"),
      onlyCancer: $("only-cancer-genes"),
      geneInput: $("gene-filter"),
      minTargets: $("min-targets"),
    },
    chunk: {
      hideNonCnv: $("hide-non-cnv-chunk"),
      includePonCnv: $("show-pon-cnv-chunk"),
      onlyCancer: $("only-cancer-genes-chunk"),
      onlySplit: $("only-split-genes-chunk"),
      geneInput: $("gene-filter-chunk"),
      minTargets: $("min-targets-chunk"),
    },
  };

  // ---------------------------------------------------------------------------
  // Column group toggles
  // ---------------------------------------------------------------------------

  const COL_GROUP_QC = ["depth_mean"];

  const COL_GROUP_PON_PREFIXES = ["pon_"];
  const COL_GROUP_PON_EXPLICIT = [];

  const COL_GROUP_CNVKIT_EXTRA = [
    "cnvkit_seg_log2",
    "cnvkit_seg_cn",
    "cnvkit_seg_cn1",
    "cnvkit_seg_cn2",
  ];

  const COL_GROUP_PURECN_EXTRA = ["purecn_C", "purecn_M", "purecn_M_flagged"];

  const DEFAULT_GROUP_STATE = {
    "report-table": { qc: false, pon: false, cnvkit: false, purecn: false },
    "chunk-table": { qc: false, pon: true, cnvkit: false, purecn: false },
  };

  function columnsMatchingPrefixes(colIndex, prefixes) {
    const cols = Object.keys(colIndex ?? {});
    return cols.filter((c) => prefixes.some((p) => c.startsWith(p)));
  }

  function getGroupColumns(colIndex, groupName) {
    if (!colIndex) return [];

    const groupDefs = {
      qc: () => COL_GROUP_QC.filter((c) => c in colIndex),
      pon: () => {
        const prefixCols = columnsMatchingPrefixes(colIndex, COL_GROUP_PON_PREFIXES);
        const explicitCols = COL_GROUP_PON_EXPLICIT.filter((c) => c in colIndex);
        return unique([...prefixCols, ...explicitCols]);
      },
      cnvkit: () => COL_GROUP_CNVKIT_EXTRA.filter((c) => c in colIndex),
      purecn: () => COL_GROUP_PURECN_EXTRA.filter((c) => c in colIndex),
    };

    return (groupDefs[groupName]?.() ?? []);
  }

  function setColumnVisibility(table, colIndex, colName, visible) {
    if (!table) return;
    const idx = colIndex?.[colName];
    if (!Number.isInteger(idx)) return;

    const displayValue = visible ? "" : "none";

    // header
    for (const hr of Array.from(table.tHead?.rows ?? [])) {
      const cell = hr.cells?.[idx];
      if (cell) cell.style.display = displayValue;
    }

    // body
    for (const tb of Array.from(table.tBodies ?? [])) {
      for (const row of Array.from(tb.rows ?? [])) {
        const cell = row.cells?.[idx];
        if (cell) cell.style.display = displayValue;
      }
    }
  }

  function applyColumnGroupsForTable(tableId, toggles) {
    const table = getTable(tableId);
    if (!table) return;

    const colIndex = getColIndex(tableId);

    const groups = [
      ["qc", toggles.qc],
      ["pon", toggles.pon],
      ["cnvkit", toggles.cnvkit],
      ["purecn", toggles.purecn],
    ];

    for (const [groupName, enabled] of groups) {
      const cols = getGroupColumns(colIndex, groupName);
      for (const c of cols) setColumnVisibility(table, colIndex, c, Boolean(enabled));
    }
  }

  function hookColumnGroupToggles() {
    const toggleEls = {
      "report-table": {
        qc: $("cols-qc-gene"),
        pon: $("cols-pon-gene"),
        cnvkit: $("cols-cnvkit-extra-gene"),
        purecn: $("cols-purecn-extra-gene"),
      },
      "chunk-table": {
        qc: $("cols-qc-chunk"),
        pon: $("cols-pon-chunk"),
        cnvkit: $("cols-cnvkit-extra-chunk"),
        purecn: $("cols-purecn-extra-chunk"),
      },
    };

    // apply defaults
    for (const [tableId, toggles] of Object.entries(toggleEls)) {
      const defaults = DEFAULT_GROUP_STATE[tableId];
      if (!defaults) continue;

      for (const [k, el] of Object.entries(toggles)) {
        if (el) el.checked = Boolean(defaults[k]);
      }
    }

    function applyAll() {
      for (const [tableId, toggles] of Object.entries(toggleEls)) {
        applyColumnGroupsForTable(tableId, {
          qc: toggles.qc?.checked,
          pon: toggles.pon?.checked,
          cnvkit: toggles.cnvkit?.checked,
          purecn: toggles.purecn?.checked,
        });
      }
    }

    applyAll();

    // bind events
    const events = [];
    for (const toggles of Object.values(toggleEls)) {
      for (const el of Object.values(toggles)) {
        if (el) events.push([el, "change"]);
      }
    }
    bindEvents(events, applyAll);
  }

  // ---------------------------------------------------------------------------
  // Row filtering (split into small predicates)
  // ---------------------------------------------------------------------------

  function parseGeneList(value) {
    return (value ?? "")
      .split(",")
      .map((g) => g.trim().toLowerCase())
      .filter(Boolean);
  }

  function getIdx(colIndex, ...names) {
    for (const n of names) {
      const idx = colIndex?.[n];
      if (Number.isInteger(idx)) return idx;
    }
    return -1;
  }

  function isNonNeutralCall(cell) {
    const v = textOf(cell);
    if (!v) return false;
    if (v === "neutral") return false;
    if (v === "na" || v === "nan" || v === ".") return false;
    return true;
  }

  function rowHasCnvOrLoh(row, colIndex) {
  const lohIdx = getIdx(colIndex, "purecn_loh_flag");
  const cnvkitIdx = getIdx(colIndex, "cnvkit_cnv_call");
  const purecnIdx = getIdx(colIndex, "purecn_cnv_call");

  if (lohIdx !== -1 && isPresentText(textOf(row.cells[lohIdx]))) return true;

  if (cnvkitIdx !== -1 && isNonNeutralCall(row.cells[cnvkitIdx])) return true;
  if (purecnIdx !== -1 && isNonNeutralCall(row.cells[purecnIdx])) return true;

  return false;
  }

  function rowIsCancerGene(row, colIndex) {
    const idx = getIdx(colIndex, "is_cancer_gene");
    return idx !== -1 && isTruthyText(textOf(row.cells[idx]));
  }

  function rowIsSplitGene(row, colIndex) {
    const idx = getIdx(colIndex, "is_gene_split");
    return idx !== -1 && isTruthyText(textOf(row.cells[idx]));
  }

  function rowHasMinTargets(row, colIndex, minTargets) {
    if (!minTargets) return true;
    const idx = getIdx(colIndex, "n.targets");
    if (idx === -1) return true;

    const nText = textOf(row.cells[idx]);
    if (!nText) return true;

    const nVal = parseFloat(nText);
    return Number.isNaN(nVal) ? true : nVal >= minTargets;
  }

  function rowMatchesGeneList(row, colIndex, geneList) {
    if (!geneList?.length) return true;
    const idx = getIdx(colIndex, "gene.symbol");
    if (idx === -1) return true;
    const gene = textOf(row.cells[idx]);
    return geneList.includes(gene);
  }

  function rowHasPonCnvCall(row, colIndex) {
    const idx = getIdx(colIndex, "pon_cnv_call");
    if (idx === -1) return false;
    const v = textOf(row.cells[idx]);
    return v === "amplification" || v === "deletion";
  }

  function shouldHideRow(row, colIndex, cfg) {
    // CNV filter
    if (cfg.hideNonCnv) {
      const cnv = rowHasCnvOrLoh(row, colIndex);
      const pon = cfg.includePonCnv ? rowHasPonCnvCall(row, colIndex) : false;
      if (!cnv && !pon) return true;
    }

    if (cfg.onlyCancer && !rowIsCancerGene(row, colIndex)) return true;
    if (cfg.onlySplit && !rowIsSplitGene(row, colIndex)) return true;
    if (!rowHasMinTargets(row, colIndex, cfg.minTargets)) return true;
    if (!rowMatchesGeneList(row, colIndex, cfg.geneList)) return true;

    return false;
  }

  function filterTable(tableId, table, cfg) {
    if (!table) return;
    const body = table.tBodies?.[0];
    if (!body) return;

    const colIndex = getColIndex(tableId);

    for (const row of Array.from(body.rows)) {
      row.style.display = shouldHideRow(row, colIndex, cfg) ? "none" : "";
    }
  }

  function buildGeneCfg() {
    const c = controls.gene;
    return {
      hideNonCnv: Boolean(c.hideNonCnv?.checked),
      includePonCnv: false,
      onlyCancer: Boolean(c.onlyCancer?.checked),
      onlySplit: false,
      geneList: parseGeneList(c.geneInput?.value),
      minTargets: parseIntSafe(c.minTargets?.value),
    };
  }

  function buildChunkCfg() {
    const c = controls.chunk;
    return {
      hideNonCnv: Boolean(c.hideNonCnv?.checked),
      includePonCnv: Boolean(c.includePonCnv?.checked),
      onlyCancer: Boolean(c.onlyCancer?.checked),
      onlySplit: Boolean(c.onlySplit?.checked),
      geneList: parseGeneList(c.geneInput?.value),
      minTargets: parseIntSafe(c.minTargets?.value),
    };
  }

  function applyFilter() {
    filterTable("report-table", geneTable, buildGeneCfg());
    filterTable("chunk-table", chunkTable, buildChunkCfg());
  }

  // bind filter events
  bindEvents(
    [
      [controls.gene.hideNonCnv, "change"],
      [controls.gene.onlyCancer, "change"],
      [controls.gene.geneInput, "input"],
      [controls.gene.minTargets, "input"],

      [controls.chunk.hideNonCnv, "change"],
      [controls.chunk.includePonCnv, "change"],
      [controls.chunk.onlyCancer, "change"],
      [controls.chunk.onlySplit, "change"],
      [controls.chunk.geneInput, "input"],
      [controls.chunk.minTargets, "input"],
    ],
    applyFilter
  );

  // ---------------------------------------------------------------------------
  // Table show/hide toggles
  // ---------------------------------------------------------------------------

  function hookShowHide(buttonId, containerId, initialVisible, textWhenVisible, textWhenHidden) {
    const btn = $(buttonId);
    const container = $(containerId);
    if (!btn || !container) return;

    let visible = Boolean(initialVisible);
    container.style.display = visible ? "block" : "none";
    btn.textContent = visible ? textWhenVisible : textWhenHidden;

    btn.addEventListener("click", () => {
      visible = !visible;
      container.style.display = visible ? "block" : "none";
      btn.textContent = visible ? textWhenVisible : textWhenHidden;
    });
  }

  hookShowHide(
    "toggle-gene-table",
    "gene-table-container",
    true,
    "Hide gene-level table",
    "Show gene-level table"
  );

  hookShowHide(
    "toggle-chunk-table",
    "chunk-table-container",
    false,
    "Hide chunk-level table",
    "Show chunk-level table"
  );

  // ---------------------------------------------------------------------------
  // Plot modal viewer
  // ---------------------------------------------------------------------------

  const modal = $("plot-modal");
  const modalImg = $("plot-modal-image");
  const modalClose = $("plot-modal-close");

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
    for (const img of document.querySelectorAll(".clickable-plot")) {
      img.addEventListener("click", () => openModal(img.src));
    }
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
  // Init
  // ---------------------------------------------------------------------------

  hookColumnGroupToggles();
  applyFilter();
})();
