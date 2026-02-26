from __future__ import annotations

# Standard library
from pathlib import Path
from typing import Dict, List, Tuple, Mapping, Iterable

# Third-party
import pandas as pd

from cnv_constants import CHR


def _load_bins_tsv(
    path: str | Path,
    *,
    chr_col: str = "chromosome",
    gene_col: str = "gene",
    rename: Mapping[str, str] | None = None,
    out_cols: Iterable[str],
) -> pd.DataFrame:
    """
    Load CNVkit-style bin table (.cnr/.cnn).

    - Rename chromosome column -> CHR ("chr")
    - Drop Antitarget bins and NA genes
    - Rename "-" -> "backbone"
    - Apply optional column renames
    - Select required out_cols (must all exist)
    """
    df = pd.read_csv(path, sep="\t")

    # normalize chromosome column name
    df = df.rename(columns={chr_col: CHR})

    # gene.symbol cleanup
    gene = df[gene_col].astype("string")
    mask_keep = gene.notna() & ~gene.isin(["Antitarget"])
    df = df.loc[mask_keep].copy()
    df["gene.symbol"] = gene.loc[mask_keep].replace("-", "backbone")

    # explode multi-gene bins
    df["gene.symbol"] = df["gene.symbol"].str.split(r"\s*,\s*", regex=True)
    df = df.explode("gene.symbol", ignore_index=True)
    df["gene.symbol"] = df["gene.symbol"].astype("string").str.strip()

    # optional renames (e.g. log2->pon_log2, spread->pon_spread)
    if rename is not None:
        df = df.rename(columns=dict(rename))

    # select required columns
    out_cols = list(out_cols)
    return df[out_cols].copy()


def load_cnr_bins(cnr_path: str | Path) -> pd.DataFrame:
    return _load_bins_tsv(
        cnr_path,
        out_cols=[CHR, "start", "end", "gene.symbol", "log2", "depth"],
    )


def load_pon_bins(pon_path: str | Path) -> pd.DataFrame:
    return _load_bins_tsv(
        pon_path,
        rename={"log2": "pon_log2", "spread": "pon_spread"},
        out_cols=[CHR, "start", "end", "pon_log2", "pon_spread", "gene.symbol", "gc"],
    )


def load_purecn_segments(pure_cn: str | Path) -> pd.DataFrame:
    """
    Load PureCN segment CSV, keep relevant columns, and assign stable ordering.

    Keeps the original `type` column as-is (no derived boolean LOH flag).
    """
    cns = pd.read_csv(pure_cn, sep=",").copy()

    keep_columns = [
        "chr",
        "start",
        "end",
        "seg.mean",
        "num.snps",
        "M",
        "M.flagged",
        "C",
        "maf.observed",
        "type",
    ]

    keep = [c for c in keep_columns if c in cns.columns]
    cns = cns[keep].copy()

    cns["loh_flag"] = cns["type"].astype("string")

    cns = cns.sort_values(["chr", "start"], kind="stable").reset_index(drop=True)
    return cns


def _load_cns_segments(cns_path: str | Path, cns_file_type: str) -> pd.DataFrame:
    """
    Load CNVkit CNS, normalize chr, keep available segment columns,
    and assign stable segment_id within each file.
    """
    cns = pd.read_csv(cns_path, sep="\t").copy()

    chr_col = "chromosome"

    if cns_file_type == "calls":
        base_cols = ["start", "end", "log2", "baf", "cn", "cn1", "cn2", "depth"]
    else:
        base_cols = ["start", "end", "log2"]
    keep = [c for c in base_cols if c in cns.columns]
    cns = cns[[chr_col] + keep].rename(columns={chr_col: "chr"}).copy()
    cns = cns.sort_values(["chr", "start"], kind="stable").reset_index(drop=True)
    cns["segment_id"] = cns.index.astype(str)
    return cns


def load_cnvkit_segments_with_raw(
    cns_path: str | Path,
    cns_init_path: str | Path,
) -> pd.DataFrame:
    """Load CNVkit called segments and attach raw_log2 from the init segments file."""

    segs = _load_cns_segments(cns_path, "calls").copy()
    segs_init = _load_cns_segments(cns_init_path, "raw").copy()

    segs_init = segs_init.rename(columns={"log2": "raw_log2"})
    segs = segs.merge(
        segs_init[["chr", "start", "end", "raw_log2"]],
        on=["chr", "start", "end"],
        how="left",
        validate="one_to_one",
    )
    return segs


def load_refgene_exons(refgene_path: str | Path) -> Dict[Tuple[str, str], dict]:
    """Load refgene and create an exon map from the longest transcript of a gene."""
    cols = [
        "gene_symbol",
        "transcript",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
    ]
    rg = pd.read_csv(refgene_path, sep="\t", header=None, names=cols)

    rg["txStart"] = rg["txStart"].astype(int)
    rg["txEnd"] = rg["txEnd"].astype(int)

    exon_map: Dict[Tuple[str, str], dict] = {}

    for (chrom, gene), sub in rg.groupby(["chrom", "gene_symbol"], sort=False):
        sub = sub.copy()
        sub["tx_len"] = sub["txEnd"] - sub["txStart"]
        best = sub.loc[sub["tx_len"].idxmax()]

        starts = [int(x) for x in str(best["exonStarts"]).split(",") if x != ""]
        ends = [int(x) for x in str(best["exonEnds"]).split(",") if x != ""]
        exons: List[Tuple[int, int]] = list(zip(starts, ends))
        if not exons:
            continue

        exons.sort(key=lambda x: x[0])
        if best["strand"] == "-":
            exons = exons[::-1]

        exon_map[(chrom, gene)] = {
            "transcript": best["transcript"],
            "strand": best["strand"],
            "txStart": int(best["txStart"]),
            "txEnd": int(best["txEnd"]),
            "exons": exons,
        }

    return exon_map


def load_cytobands(path: str | Path) -> pd.DataFrame:
    """Load cytoband UCSC file; return normalized df with integer coords."""
    cols = ["chr", "chromStart", "chromEnd", "name", "gieStain"]
    cyto = pd.read_csv(path, sep="\t", header=None, names=cols)
    cyto["chr"] = cyto["chr"].astype(str).str.replace("^chr", "", regex=True)
    cyto["start_int"] = cyto["chromStart"].astype(int)
    cyto["end_int"] = cyto["chromEnd"].astype(int)
    return cyto


def load_cancer_gene_set(
    path: str | Path,
    min_occurrence: int = 1,
    oncokb_source: str = "ONCOKB",
) -> set[str]:
    """
    Load cancer gene list TSV (GENE / OCCURRENCE / SOURCE) into a set of gene symbols.

    Rules
    -----
    - Always include non-empty genes from non-ONCOKB sources (e.g. cust000), regardless of OCCURRENCE.
    - For ONCOKB rows only, require OCCURRENCE >= min_occurrence (numeric; non-numeric treated as 0).

    Parameters
    ----------
    path
        Path to TSV with columns: GENE, OCCURRENCE, SOURCE
    min_occurrence
        Minimum occurrence threshold applied only to ONCOKB source rows.
    oncokb_source
        Source label used for OncoKB rows (default "ONCOKB").

    Returns
    -------
    set[str]
        Selected gene symbols.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    gene_col = "GENE"
    occ_col = "OCCURRENCE"
    src_col = "SOURCE"

    # Normalize
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df[src_col] = df[src_col].astype(str).str.strip().str.upper()

    # Parse occurrence (NA/non-numeric -> 0)
    occ = pd.to_numeric(df[occ_col], errors="coerce").fillna(0).astype(int)

    is_oncokb = df[src_col].eq(str(oncokb_source).strip().upper())

    # ONCOKB: apply threshold; non-ONCOKB: always keep
    keep = (~is_oncokb) | (occ >= int(min_occurrence))

    selected = df.loc[keep, gene_col]
    selected = selected[selected != ""]
    return set(selected)
