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
      "report-table": globalThis.colIndexSegments,
      "region-table": globalThis.colIndexRegions,
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

  function getVisibleColumnIndexes(table) {
    if (!table?.tHead?.rows?.length) return [];
    const headerCells = Array.from(table.tHead.rows[0].cells);
    const idxs = [];
    for (let i = 0; i < headerCells.length; i++) {
      // If a column group toggle hid it, display will be "none"
      if (headerCells[i].style.display !== "none") idxs.push(i);
    }
    return idxs;
  }

  function tableToTSVVisible(table) {
    if (!table) return "";

    const body = table.tBodies?.[0];
    if (!body) return "";

    const colIdxs = getVisibleColumnIndexes(table);
    if (!colIdxs.length) return "";

    const lines = [];

    // header row
    const headerCells = Array.from(table.tHead.rows[0].cells);
    lines.push(colIdxs.map((i) => (headerCells[i]?.textContent ?? "").trim()).join("\t"));

    // visible body rows
    for (const row of Array.from(body.rows)) {
      if (row.style.display === "none") continue;

      const cells = Array.from(row.cells);
      lines.push(colIdxs.map((i) => (cells[i]?.textContent ?? "").trim()).join("\t"));
    }

    return lines.join("\n");
  }

  async function copyTextToClipboard(text) {
    if (!text) return;

    // Prefer async clipboard API
    if (navigator.clipboard?.writeText) {
      await navigator.clipboard.writeText(text);
      return;
    }

    // Fallback for older browsers / stricter contexts
    const ta = document.createElement("textarea");
    ta.value = text;
    ta.style.position = "fixed";
    ta.style.left = "-9999px";
    ta.style.top = "0";
    document.body.appendChild(ta);
    ta.focus();
    ta.select();
    document.execCommand("copy");
    document.body.removeChild(ta);
  }

  // ---------------------------------------------------------------------------
  // Tables + controls
  // ---------------------------------------------------------------------------
  const segTable = getTable("report-table");
  const regionTable = getTable("region-table");

  const controls = {
    seg: {
      hideNonCnv: $("hide-non-cnv-seg"),
      geneInput: $("gene-filter-seg"),
      rangeChr: $("range-chr-seg"),
      rangeStart: $("range-start-seg"),
      rangeEnd: $("range-end-seg"),
    },
    region: {
      hideNonCnv: $("hide-non-cnv-region"),
      includePonCnv: $("show-pon-cnv-region"),
      onlyCancer: $("only-cancer-genes-region"),
      geneInput: $("gene-filter-region"),
      minTargets: $("min-targets-region"),
      rangeChr: $("range-chr-region"),
      rangeStart: $("range-start-region"),
      rangeEnd: $("range-end-region"),
      onlyPonIndication: $("only-pon-indication-region"),
    },
  };

  // ---------------------------------------------------------------------------
  // Column group toggles
  // ---------------------------------------------------------------------------

  const COL_GROUP_QC = ["depth_mean", "min_log2", "max_log2"];

  const COL_GROUP_PON_PREFIXES = ["pon_"];
  const COL_GROUP_PON_EXPLICIT = [];

  const COL_GROUP_CNVKIT_EXTRA = [
    "cnvkit_adjusted_log2",
    "cnvkit_seg_cn",
    "cnvkit_seg_cn1",
    "cnvkit_seg_cn2",
    "cnvkit_seg_depth",
  ];

  const COL_GROUP_PURECN_EXTRA = ["purecn_C", "purecn_M", "purecn_M_flagged", "purecn_num_snps"];

  const DEFAULT_GROUP_STATE = {
    "report-table": { qc: false, pon: false, cnvkit: false, purecn: false },
    "region-table": { qc: false, pon: true, cnvkit: false, purecn: false },
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
        cnvkit: $("cols-cnvkit-extra-gene"),
        purecn: $("cols-purecn-extra-gene"),
      },
      "region-table": {
        qc: $("cols-qc-region"),
        pon: $("cols-pon-region"),
        cnvkit: $("cols-cnvkit-extra-region"),
        purecn: $("cols-purecn-extra-region"),
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
  // LOH for both tables: PureCN type present (non-empty)
  const lohIdx = getIdx(colIndex, "purecn_type");
  if (lohIdx !== -1 && isPresentText(textOf(row.cells[lohIdx]))) return true;

  // Segments table: unified call column
  const cnvIdx = getIdx(colIndex, "cnv_call");
  if (cnvIdx !== -1 && isNonNeutralCall(row.cells[cnvIdx])) return true;

  // region table (or legacy): separate call columns
  const cnvkitIdx = getIdx(colIndex, "cnvkit_cnv_call");
  const purecnIdx = getIdx(colIndex, "purecn_cnv_call");
  if (cnvkitIdx !== -1 && isNonNeutralCall(row.cells[cnvkitIdx])) return true;
  if (purecnIdx !== -1 && isNonNeutralCall(row.cells[purecnIdx])) return true;

  return false;
}

  function rowIsCancerGene(row, colIndex) {
    const idx = getIdx(colIndex, "is_cancer_gene");
    return idx !== -1 && isTruthyText(textOf(row.cells[idx]));
  }

  function rowHasPonIndication(row, colIndex) {
    const idx = getIdx(colIndex, "pon_region_indication");
    if (idx === -1) return false;

    const v = textOf(row.cells[idx]);
    if (!v) return false;
    if (v === "neutral") return false;
    if (v === "na" || v === "nan" || v === "<na>" || v === ".") return false;

    return true; // any non-neutral indication counts
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

  // prefer segment column if present
  const idxList = getIdx(colIndex, "gene.symbol");
  if (idxList !== -1) {
    const s = (row.cells[idxList]?.textContent ?? "").trim().toLowerCase();
    if (!s) return false;

    // "TP53,EGFR,MYC" -> ["tp53","egfr","myc"]
    const genesInRow = s
      .split(",")
      .map((g) => g.trim().toLowerCase())
      .filter(Boolean);

    // any match is enough
    return genesInRow.some((g) => geneList.includes(g));
  }

  // fallback to old single-gene column
  const idxOne = getIdx(colIndex, "gene.symbol");
  if (idxOne === -1) return true;

  const gene = textOf(row.cells[idxOne]);
  return geneList.includes(gene);
}

  function rowHasPonCnvCall(row, colIndex) {
    const idx = getIdx(colIndex, "pon_cnv_call");
    if (idx === -1) return false;
    const v = textOf(row.cells[idx]);
    return v === "amplification" || v === "deletion";
  }

  function shouldHideRow(row, colIndex, cfg) {
    // apply genomic range filter
    if (cfg.range && !rowInRange(row, colIndex, cfg.range)) return true;

    // Additive "OR" gate for region table:
    if (cfg.isRegionTable) {
      const wantCnv = Boolean(cfg.hideNonCnv);
      const wantPonInd = Boolean(cfg.onlyPonIndication);

      if (wantCnv || wantPonInd) {
        const matchCnv = wantCnv
          ? (rowHasCnvOrLoh(row, colIndex) ||
             (cfg.includePonCnv ? rowHasPonCnvCall(row, colIndex) : false))
          : false;

        const matchPonInd = wantPonInd ? rowHasPonIndication(row, colIndex) : false;

        if (!matchCnv && !matchPonInd) return true;
      }
    } else {
      if (cfg.hideNonCnv) {
        const cnv = rowHasCnvOrLoh(row, colIndex);
        const pon = cfg.includePonCnv ? rowHasPonCnvCall(row, colIndex) : false;
        if (!cnv && !pon) return true;
      }
    }

    if (cfg.onlyCancer && !rowIsCancerGene(row, colIndex)) return true;
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

  function buildSegCfg() {
    const c = controls.seg;
    return {
      isRegionTable: false,
      hideNonCnv: Boolean(c.hideNonCnv?.checked),
      includePonCnv: false,
      onlyCancer: false,          // segment UI doesn’t have this right now
      geneList: parseGeneList(c.geneInput?.value),
      minTargets: 0,              // segment UI doesn’t have this right now
      range: {
        chr: c.rangeChr?.value ?? "",
        start: parseFloatSafe(c.rangeStart?.value),
        end: parseFloatSafe(c.rangeEnd?.value),
      },
    };
  }

  function buildRegionCfg() {
    const c = controls.region;
    return {
      isRegionTable: true,
      hideNonCnv: Boolean(c.hideNonCnv?.checked),
      includePonCnv: Boolean(c.includePonCnv?.checked),
      onlyPonIndication: Boolean(c.onlyPonIndication?.checked),
      onlyCancer: Boolean(c.onlyCancer?.checked),
      geneList: parseGeneList(c.geneInput?.value),
      minTargets: parseIntSafe(c.minTargets?.value),
      range: {
        chr: c.rangeChr?.value ?? "",
        start: parseFloatSafe(c.rangeStart?.value),
        end: parseFloatSafe(c.rangeEnd?.value),
      },
    };
  }

  function applyFilter() {
    filterTable("report-table", segTable, buildSegCfg());
    filterTable("region-table", regionTable, buildRegionCfg());
  }

  // bind filter events
  bindEvents(
    [
      [controls.seg.hideNonCnv, "change"],
      [controls.seg.geneInput, "input"],
      [controls.seg.rangeChr, "input"],
      [controls.seg.rangeStart, "input"],
      [controls.seg.rangeEnd, "input"],

      [controls.region.hideNonCnv, "change"],
      [controls.region.includePonCnv, "change"],
      [controls.region.onlyCancer, "change"],
      [controls.region.geneInput, "input"],
      [controls.region.minTargets, "input"],
      [controls.region.rangeChr, "input"],
      [controls.region.rangeStart, "input"],
      [controls.region.rangeEnd, "input"],
      [controls.region.onlyPonIndication, "change"],
    ],
    applyFilter
  );

  const copySegBtn = $("copy-segment-table");
  if (copySegBtn) {
    copySegBtn.addEventListener("click", async () => {
      const tsv = tableToTSVVisible(segTable);
      await copyTextToClipboard(tsv);
      copySegBtn.textContent = "Copied!";
      setTimeout(() => (copySegBtn.textContent = "Copy visible rows"), 1000);
    });
  }

  const copyRegionBtn = $("copy-region-table");
  if (copyRegionBtn) {
    copyRegionBtn.addEventListener("click", async () => {
      const tsv = tableToTSVVisible(regionTable);
      await copyTextToClipboard(tsv);
      copyRegionBtn.textContent = "Copied!";
      setTimeout(() => (copyRegionBtn.textContent = "Copy visible rows"), 1000);
    });
  }

  function normChr(s) {
  s = (s ?? "").trim().toLowerCase();
  if (s.startsWith("chr")) s = s.slice(3);
  return s;
  }

  function parseFloatSafe(v) {
    const n = parseFloat(v);
    return Number.isFinite(n) ? n : null;
  }

  function rowInRange(row, colIndex, range) {
    const wantChr = normChr(range.chr);
    const wantStart = range.start;
    const wantEnd = range.end;

  // no range filter set
  if (!wantChr && wantStart == null && wantEnd == null) return true;

  // find columns
  const chrIdx = getIdx(colIndex, "chr", "chromosome");
  const startIdx = getIdx(colIndex, "region_start", "start");
  const endIdx = getIdx(colIndex, "region_end", "end");

  if (chrIdx === -1 || startIdx === -1 || endIdx === -1) {
    // can't apply range if columns aren't available
    return true;
  }

  const rowChr = normChr(textOf(row.cells[chrIdx]));
  const rowStart = parseFloatSafe(textOf(row.cells[startIdx]));
  const rowEnd = parseFloatSafe(textOf(row.cells[endIdx]));

  if (!rowChr || rowStart == null || rowEnd == null) return true;

  // chr must match if specified
  if (wantChr && rowChr !== wantChr) return false;

  // overlap check:
  // keep row if it overlaps [wantStart, wantEnd]
  const qStart = wantStart != null ? wantStart : -Infinity;
  const qEnd = wantEnd != null ? wantEnd : Infinity;

  return rowEnd >= qStart && rowStart <= qEnd;
  }

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
    "toggle-segment-table",
    "segment-table-container",
    true,
    "Hide segment table",
    "Show segment table"
  );

  hookShowHide(
    "toggle-region-table",
    "region-table-container",
    false,
    "Hide region-level table",
    "Show region-level table"
  );

  hookShowHide(
    "toggle-seg-col-glossary",
    "seg-col-glossary-container",
    false,
    "Hide segment glossary",
    "Show segment glossary"
  );

  hookShowHide(
    "toggle-region-col-glossary",
    "region-col-glossary-container",
    false,
    "Hide region glossary",
    "Show region glossary"
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
