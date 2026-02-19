from __future__ import annotations

# Standard library
from pathlib import Path
from typing import Dict, List, Tuple, Mapping, Iterable

# Third-party
import pandas as pd


from cnv_report_utils import strip_chr_prefix
from cnv_constants import CHR



def _load_bins_tsv(
    path: str | Path,
    *,
    chr_col: str = "chromosome",
    gene_col: str = "gene",
    explode_genes: bool = True,
    rename: Mapping[str, str] | None = None,
    out_cols: Iterable[str] | None = None,
) -> pd.DataFrame:
    """
    Generic loader for CNVkit-style bin tables (.cnr, .cnn, etc).

    - Normalizes chromosome column to CHR ("chr").
    - Creates/normalizes "gene.symbol".
    - Drops Antitarget bins.
    - Renames "-" bins to "backbone".
    - Optionally explodes multi-gene bins.
    - Optional renaming + column selection.
    """
    df = pd.read_csv(path, sep="\t")

    # base renames (chromosome -> chr)
    df = df.rename(columns={chr_col: CHR})

    if gene_col in df.columns:
        gene = df[gene_col].astype("string")

        # Drop only Antitarget and NA genes
        mask_keep = gene.notna() & ~gene.isin(["Antitarget"])
        df = df.loc[mask_keep].copy()
        gene = gene.loc[mask_keep]

        # Rename "-" to "backbone"
        gene = gene.replace("-", "backbone")

        df["gene.symbol"] = gene
    else:
        df["gene.symbol"] = pd.Series(
            [pd.NA] * len(df), dtype="string", index=df.index
        )

    if explode_genes and gene_col in df.columns:
        df["gene.symbol"] = (
            df["gene.symbol"]
            .astype("string")
            .str.split(r"\s*,\s*", regex=True)
        )
        df = df.explode("gene.symbol", ignore_index=True)
        df["gene.symbol"] = df["gene.symbol"].astype("string").str.strip()

    # drop empty symbols from malformed inputs like "TP53,,RB1"
    df = df[df["gene.symbol"].notna() & df["gene.symbol"].ne("")].copy()

    # apply optional renames (e.g. log2->pon_log2, spread->pon_spread)
    if rename:
        df = df.rename(columns=dict(rename))

    if out_cols is not None:
        out_cols = list(out_cols)
        missing = [c for c in out_cols if c not in df.columns]
        if missing:
            raise KeyError(f"Missing required columns in {path}: {missing}")
        df = df[out_cols].copy()

    return df


def load_cnr_bins(cnr_path: str | Path) -> pd.DataFrame:
    return _load_bins_tsv(
        cnr_path,
        explode_genes=True,
        out_cols=[CHR, "start", "end", "gene.symbol", "log2", "depth", "weight"],
    )


def load_pon_bins(pon_path: str | Path) -> pd.DataFrame:
    return _load_bins_tsv(
        pon_path,
        explode_genes=True,  # keep behavior consistent with CNR (safe either way)
        rename={"log2": "pon_log2", "spread": "pon_spread"},
        out_cols=[CHR, "start", "end", "pon_log2", "pon_spread", "gene.symbol"],
    )


def load_purecn_segments(pure_cn: str | Path) -> pd.DataFrame:
    """
    Load PureCN segment file, keep relevant columns,
    assign stable ordering, and derive LOH flag.

    The loh_flag column:
        - True  -> if type contains "LOH"
        - False -> if type present but not LOH
        - <NA>  -> if type is missing
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

    # ------------------------------------------------------------------
    # Derive LOH flag from the type column
    # ------------------------------------------------------------------
    loh_flag_col = "loh_flag"
    cns[loh_flag_col] = pd.array([pd.NA] * len(cns), dtype="boolean")

    if "type" in cns.columns:
        loh_mask = (
            cns["type"]
            .astype("string")          # preserves <NA>
            .str.upper()
            .str.contains("LOH", na=pd.NA)
        )
        cns[loh_flag_col] = loh_mask.astype("boolean")

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
