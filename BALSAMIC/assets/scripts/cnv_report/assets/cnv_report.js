(function () {
  const geneTable  = document.getElementById("report-table");
  const chunkTable = document.getElementById("chunk-table");

  // Gene-level controls
  const gHideNonCnv = document.getElementById("hide-non-cnv");
  const gOnlyCancer = document.getElementById("only-cancer-genes");
  const gGeneInput  = document.getElementById("gene-filter");
  const gMinTargets = document.getElementById("min-targets");

  // Chunk-level controls
  const cHideNonCnv    = document.getElementById("hide-non-cnv-chunk");
  const cIncludePonCnv = document.getElementById("show-pon-cnv-chunk");
  const cOnlyCancer    = document.getElementById("only-cancer-genes-chunk");
  const cOnlySplit     = document.getElementById("only-split-genes-chunk");
  const cGeneInput     = document.getElementById("gene-filter-chunk");
  const cMinTargets    = document.getElementById("min-targets-chunk");

  function normalizeCellText(td) {
    return (td && td.textContent ? td.textContent : "")
      .trim()
      .toLowerCase();
  }

  function isCnvRow(row, colIndex) {
    const lohColIdx    = (colIndex["loh_flag"] ?? -1);
    const cnvkitColIdx = (colIndex["cnvkit_cnv_call"] ?? -1);
    const purecnColIdx = (colIndex["purecn_cnv_call"] ?? -1);

    let isCnv = false;

    if (lohColIdx !== -1) {
      const lohVal = normalizeCellText(row.cells[lohColIdx]);
      if (lohVal === "true") isCnv = true;
    }

    if (!isCnv && cnvkitColIdx !== -1) {
      const cVal = normalizeCellText(row.cells[cnvkitColIdx]);
      if (cVal === "deletion" || cVal === "amplification") isCnv = true;
    }

    if (!isCnv && purecnColIdx !== -1) {
      const pVal = normalizeCellText(row.cells[purecnColIdx]);
      if (pVal === "deletion" || pVal === "amplification") isCnv = true;
    }

    return isCnv;
  }

  function isCancerGeneRow(row, colIndex) {
    const idx = (colIndex["is_cancer_gene"] ?? -1);
    if (idx === -1) return false;
    const val = normalizeCellText(row.cells[idx]);
    return (val === "true" || val === "1" || val === "yes");
  }

  function isSplitGeneRow(row, colIndex) {
    const idx = (colIndex["is_gene_split"] ?? -1);
    if (idx === -1) return false;
    const val = normalizeCellText(row.cells[idx]);
    return (val === "true" || val === "1" || val === "yes");
  }

  function parseGeneList(value) {
    return (value || "")
      .split(",")
      .map(g => g.trim().toLowerCase())
      .filter(g => g.length > 0);
  }

  function filterTable(table, colIndex, cfg) {
    if (!table) return;
    const body = table.tBodies[0];
    if (!body) return;

    const geneColIdx     = (colIndex["gene.symbol"] ?? -1);
    const nTargetsColIdx = (colIndex["n.targets"] ?? -1);
    const ponCnvColIdx   = (colIndex["pon_cnv_call"] ?? -1);

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
        const matches = cfg.geneList.some(g => geneText === g);
        if (!matches) hideRow = true;
      }

      row.style.display = hideRow ? "none" : "";
    }
  }

  function filterGeneTable() {
    if (!geneTable) return;
    const cfg = {
      hideNonCnv:   !!(gHideNonCnv && gHideNonCnv.checked),
      includePonCnv: false,
      onlyCancer:   !!(gOnlyCancer && gOnlyCancer.checked),
      onlySplit:    false,
      geneList:     parseGeneList(gGeneInput && gGeneInput.value),
      minTargets:   gMinTargets ? (parseInt(gMinTargets.value, 10) || 0) : 0,
    };
    filterTable(geneTable, (typeof colIndexGene !== "undefined" ? colIndexGene : {}), cfg);
  }

  function filterChunkTable() {
    if (!chunkTable) return;
    const cfg = {
      hideNonCnv:    !!(cHideNonCnv && cHideNonCnv.checked),
      includePonCnv: !!(cIncludePonCnv && cIncludePonCnv.checked),
      onlyCancer:    !!(cOnlyCancer && cOnlyCancer.checked),
      onlySplit:     !!(cOnlySplit && cOnlySplit.checked),
      geneList:      parseGeneList(cGeneInput && cGeneInput.value),
      minTargets:    cMinTargets ? (parseInt(cMinTargets.value, 10) || 0) : 0,
    };
    filterTable(chunkTable, (typeof colIndexChunk !== "undefined" ? colIndexChunk : {}), cfg);
  }

  function applyFilter() {
    filterGeneTable();
    filterChunkTable();
  }

  if (gHideNonCnv) gHideNonCnv.addEventListener("change", applyFilter);
  if (gOnlyCancer) gOnlyCancer.addEventListener("change", applyFilter);
  if (gGeneInput)  gGeneInput.addEventListener("input", applyFilter);
  if (gMinTargets) gMinTargets.addEventListener("input", applyFilter);

  if (cHideNonCnv)    cHideNonCnv.addEventListener("change", applyFilter);
  if (cIncludePonCnv) cIncludePonCnv.addEventListener("change", applyFilter);
  if (cOnlyCancer)    cOnlyCancer.addEventListener("change", applyFilter);
  if (cOnlySplit)     cOnlySplit.addEventListener("change", applyFilter);
  if (cGeneInput)     cGeneInput.addEventListener("input", applyFilter);
  if (cMinTargets)    cMinTargets.addEventListener("input", applyFilter);

  applyFilter();

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

  const modal      = document.getElementById("plot-modal");
  const modalImg   = document.getElementById("plot-modal-image");
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
    imgs.forEach(img => {
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
})();