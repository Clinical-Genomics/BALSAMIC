from __future__ import annotations

# Standard library
from pathlib import Path
from typing import Dict, List, Tuple, Mapping, Iterable, Sequence, Literal
import math
from dataclasses import dataclass

# Third-party
import pandas as pd
from cyvcf2 import VCF


from cnv_constants import CHR


def coerce_columns(
    df: pd.DataFrame,
    *,
    ints: Sequence[str] = (),
    floats: Sequence[str] = (),
    strings: Sequence[str] = (),
    booleans: Sequence[str] = (),
    bool_map: Mapping[str, Mapping[object, object]] | None = None,
    required: Sequence[str] = (),
) -> pd.DataFrame:
    """
    Coerce dataframe columns to predictable dtypes if present.

    - ints: nullable Int64
    - floats: numeric float (NaN allowed)
    - strings: pandas StringDtype
    - booleans: pandas BooleanDtype (nullable)
    - bool_map: per-column replacement map before boolean casting
    - required: raise ValueError if missing
    """
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df.copy()

    for c in strings:
        if c in df.columns:
            df[c] = df[c].astype("string")

    for c in ints:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")

    for c in floats:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Float64")

    for c in booleans:
        if c in df.columns:
            if bool_map and c in bool_map:
                df[c] = df[c].replace(bool_map[c])
            df[c] = df[c].astype("boolean")

    return df


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
    df = pd.read_csv(
        path,
        sep="\t",
        dtype={
            chr_col: "string",
            gene_col: "string",
        },
        low_memory=False,
    )

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

    if rename is not None:
        df = df.rename(columns=dict(rename))

    df = coerce_columns(
        df,
        ints=("start", "end"),
        floats=("log2", "gc", "depth", "pon_log2", "pon_spread"),
        strings=(CHR, "gene.symbol"),
    )

    out_cols = list(out_cols)
    df = coerce_columns(df, required=out_cols)  # nice error message
    return df[out_cols].copy()


def load_cnr_bins(cnr_path: str) -> pd.DataFrame:
    return _load_bins_tsv(
        cnr_path,
        out_cols=[CHR, "start", "end", "gene.symbol", "log2", "depth"],
    )


def load_pon_bins(pon_path: str) -> pd.DataFrame:
    return _load_bins_tsv(
        pon_path,
        rename={"log2": "pon_log2", "spread": "pon_spread"},
        out_cols=[CHR, "start", "end", "pon_log2", "pon_spread", "gene.symbol", "gc"],
    )


def load_purecn_segments(
    pure_cn: str,
    *,
    prefix: str = "purecn_",
) -> pd.DataFrame:
    """
    Load PureCN LOHregions/segments CSV and return a table with prefixed columns
    suitable for direct annotation merges.

    Outputs (when present in input):
      - chr, purecn_seg_start, purecn_seg_end
      - purecn_seg_mean, purecn_num_snps, purecn_M, purecn_M_flagged, purecn_C,
        purecn_maf_observed, purecn_type
      - purecn_type (nullable boolean): type contains 'LOH'
    """
    cns = pd.read_csv(pure_cn, sep=",", low_memory=False)

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

    cns = coerce_columns(
        cns,
        strings=("chr", "type"),
        ints=("start", "end", "num.snps"),
        floats=("seg.mean", "maf.observed", "M", "C"),
        booleans=("M.flagged",),
        bool_map={"M.flagged": {"TRUE": True, "FALSE": False, "NA": pd.NA}},
    )

    # Rename columns into your prefixed schema
    rename_map = {
        "start": f"{prefix}seg_start",
        "end": f"{prefix}seg_end",
        "seg.mean": f"{prefix}seg_mean_log2",
        "num.snps": f"{prefix}num_snps",
        "M": f"{prefix}M",
        "M.flagged": f"{prefix}M_flagged",
        "C": f"{prefix}C",
        "maf.observed": f"{prefix}maf_observed",
        "type": f"{prefix}type",
    }
    cns = cns.rename(columns={k: v for k, v in rename_map.items() if k in cns.columns})

    # Append "(unreliable)" if PureCN flagged the M estimate
    type_col = f"{prefix}type"
    flag_col = f"{prefix}M_flagged"

    flagged = cns[flag_col].fillna(False)
    cns[type_col] = cns[type_col].astype("string")
    cns.loc[flagged & cns[type_col].notna(), type_col] = (
        cns.loc[flagged & cns[type_col].notna(), type_col] + " (unreliable)"
    )

    return cns.sort_values(["chr", f"{prefix}seg_start"], kind="stable").reset_index(
        drop=True
    )


CnsType = Literal["calls", "raw"]


def _load_cns_segments(
    cns_path: str | Path,
    cns_file_type: CnsType,
    *,
    prefix: str = "cnvkit_",
) -> pd.DataFrame:
    """
    Load CNVkit CNS and return prefixed segment columns.

    Output columns (depending on file type / availability):
      - chr
      - {prefix}seg_start, {prefix}seg_end
      - {prefix}seg_log2
      - (calls only, if present) {prefix}seg_baf, {prefix}seg_cn, {prefix}seg_cn1, {prefix}seg_cn2, {prefix}seg_depth
    """
    df = pd.read_csv(cns_path, sep="\t", low_memory=False)

    chr_col = "chromosome"

    base_cols = (
        ["start", "end", "log2", "baf", "cn", "cn1", "cn2", "depth"]
        if cns_file_type == "calls"
        else ["start", "end", "log2"]
    )
    keep = [c for c in base_cols if c in df.columns]

    df = df[[chr_col] + keep].rename(columns={chr_col: "chr"}).copy()

    # dtype normalization (unprefixed, while easy)
    df = coerce_columns(
        df,
        strings=("chr",),
        ints=("start", "end"),
        floats=("log2", "baf", "cn", "cn1", "cn2", "depth"),
        required=("chr", "start", "end", "log2"),
    )

    df = df.sort_values(["chr", "start"], kind="stable").reset_index(drop=True)

    # Rename into prefixed schema
    rename_map = {
        "start": f"{prefix}seg_start",
        "end": f"{prefix}seg_end",
        "log2": f"{prefix}seg_log2",
        "baf": f"{prefix}seg_baf",
        "cn": f"{prefix}seg_cn",
        "cn1": f"{prefix}seg_cn1",
        "cn2": f"{prefix}seg_cn2",
        "depth": f"{prefix}seg_depth",
    }
    df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

    return df


def load_cnvkit_segments_with_raw(
    cns_path: str,
    cns_init_path: str,
    *,
    prefix: str = "cnvkit_",
) -> pd.DataFrame:
    """
    Load CNVkit called segments and attach raw log2 from init segments.

    Output includes (if present):
      - chr
      - {prefix}seg_start, {prefix}seg_end
      - {prefix}seg_log2 (called)
      - {prefix}raw_log2 (from init)
    """
    segs = _load_cns_segments(cns_path, "calls", prefix=prefix)
    init = _load_cns_segments(cns_init_path, "raw", prefix=prefix)

    # init has {prefix}seg_log2; rename that to {prefix}raw_log2 for clarity
    seg_log2_col = f"{prefix}seg_log2"
    raw_log2_col = f"{prefix}seg_raw_log2"
    init = init.rename(columns={seg_log2_col: raw_log2_col})

    # Merge on chr + prefixed coords
    on_cols = ["chr", f"{prefix}seg_start", f"{prefix}seg_end"]
    init_cols = on_cols + ([raw_log2_col] if raw_log2_col in init.columns else [])

    out = segs.merge(
        init[init_cols],
        on=on_cols,
        how="left",
        validate="one_to_one",
    )
    return out


def load_refgene_exons(refgene_path: str) -> Dict[Tuple[str, str], dict]:
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


def load_cytobands(path: str) -> pd.DataFrame:
    """Load cytoband UCSC file; return normalized df with integer coords."""
    cols = ["chr", "chromStart", "chromEnd", "name", "gieStain"]
    cyto = pd.read_csv(path, sep="\t", header=None, names=cols)
    cyto["chr"] = cyto["chr"].astype(str).str.replace("^chr", "", regex=True)
    cyto["start_int"] = cyto["chromStart"].astype(int)
    cyto["end_int"] = cyto["chromEnd"].astype(int)
    return cyto


def load_cancer_gene_set(
    path: str,
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


def load_vcf_with_vaf(
    vcf_path: str | Path,
    chr_order: list[str],
) -> pd.DataFrame:
    """
    Load a single-sample VCF using cyvcf2 and compute VAF from AD (first ALT only).

    Returns DataFrame columns:
      CHROM, POS, VAF
    """
    chr_set = set(map(str, chr_order))
    rows: list[dict[str, float | int | str]] = []

    vcf = VCF(str(vcf_path))

    for rec in vcf:
        chrom = str(rec.CHROM)
        if chrom.startswith("chr"):
            chrom = chrom[3:]

        if chrom not in chr_set:
            continue

        pos = int(rec.POS)

        vaf = math.nan
        ad = rec.format("AD")  # typically shape (n_samples, n_alleles) or None
        if ad is not None and len(ad) > 0:
            try:
                ad0 = ad[0]  # first/only sample
                # ad0 should be like [ref, alt1, alt2, ...]
                if ad0 is not None and len(ad0) >= 2:
                    ref_c = int(ad0[0])
                    alt_c = int(ad0[1])
                    total = ref_c + alt_c
                    vaf = (alt_c / total) if total > 0 else math.nan
            except (TypeError, ValueError, IndexError):
                vaf = math.nan

        rows.append({"CHROM": chrom, "POS": pos, "VAF": float(vaf)})

    return pd.DataFrame(rows, columns=["CHROM", "POS", "VAF"])
