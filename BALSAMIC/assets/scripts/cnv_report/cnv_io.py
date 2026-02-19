from __future__ import annotations

# Standard library
from pathlib import Path
from typing import Dict, List, Tuple

# Third-party
import pandas as pd


from cnv_report_utils import strip_chr_prefix, explode_multigene_bins
from cnv_constants import CHR


def load_cnr_bins(cnr_path: str | Path) -> pd.DataFrame:
    """Load CNVkit CNR, normalize chr, validate columns, explode genes."""
    cnr = pd.read_csv(cnr_path, sep="\t").copy()

    chr_col = "chromosome"

    cnr = explode_multigene_bins(cnr)
    cnr = cnr.rename(columns={chr_col: CHR})
    return cnr


def load_pon_bins(pon_path: str | Path) -> pd.DataFrame:
    """
    Load CNVkit PON .cnn, normalize chr, validate required columns,
    rename log2/spread -> pon_log2/pon_spread.
    """
    pon = pd.read_csv(pon_path, sep="\t").copy()

    chr_col = "chromosome"

    pon = pon.rename(columns={chr_col: CHR, "log2": "pon_log2", "spread": "pon_spread"})
    return pon[["chr", "start", "end", "pon_log2", "pon_spread"]]


def load_purecn_segments(pure_cn: str | Path) -> pd.DataFrame:
    """
    Load CNVkit CNS, normalize chr, keep available segment columns,
    and assign stable segment_id within each file.
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

    rg["chrom"] = rg["chrom"].astype(str).str.replace("^chr", "", regex=True)
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
    cyto["chr"] = strip_chr_prefix(cyto["chr"])
    cyto["start_int"] = cyto["chromStart"].astype(int)
    cyto["end_int"] = cyto["chromEnd"].astype(int)
    return cyto


def load_cancer_gene_set(
    path: str | Path, min_occurrence: int = 1, only_annotated: bool = True
) -> set[str]:
    """
    Load and filter a cancer gene list TSV into a set of gene symbols.

    Applies an occurrence threshold, optional OncoKB annotation filter,
    and optional restriction to classic cancer gene types (ONCOGENE / TSG)
    when the relevant columns exist.

    Parameters
    ----------
    path : str | Path
        Path to cancer gene list TSV file.
    min_occurrence : int, optional
        Minimum occurrence count required to include a gene (default: 1).
    only_annotated : bool, optional
        If True, require OncoKB annotation = YES when the column is present
        (default: True).

    Returns
    -------
    set[str]
        Gene symbols passing filters.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = df.columns.str.strip()

    symbol_col = "Hugo Symbol"
    occ_col = "# of occurrence within resources (Column J-P)"
    onco_col = "OncoKB Annotated"
    type_col = "Gene Type"

    if symbol_col not in df.columns or occ_col not in df.columns:
        missing = [c for c in (symbol_col, occ_col) if c not in df.columns]
        raise ValueError(f"Cancer gene TSV missing required columns: {missing}")

    df[occ_col] = pd.to_numeric(df[occ_col], errors="coerce").fillna(0).astype(int)

    mask = df[occ_col] >= int(min_occurrence)

    if only_annotated and onco_col in df.columns:
        mask &= df[onco_col].astype(str).str.upper().eq("YES")

    if type_col in df.columns:
        mask &= df[type_col].isin(["ONCOGENE", "TSG"])

    selected = df.loc[mask, symbol_col].astype(str).str.strip()
    selected = selected[selected != ""]
    return set(selected)
