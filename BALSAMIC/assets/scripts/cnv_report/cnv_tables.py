from __future__ import annotations

# Standard library
from typing import Dict, Mapping, Tuple

# Third-party
import numpy as np
import pandas as pd

# Local
from BALSAMIC.constants.analysis import Gender
from cnv_constants import TableSpec, GENE_TABLE_SPEC, SEGMENT_TABLE_SPEC


def finalize_table(df: pd.DataFrame, spec: TableSpec) -> pd.DataFrame:
    """
    Finalize a table based on a TableSpec: rename, reorder, round floats, stable-sort by interval.
    """
    out = df.copy()

    # 1) Rename (only if present)
    if spec.renames:
        present = {
            k: v
            for k, v in spec.renames.items()
            if k in out.columns and v not in out.columns
        }
        if present:
            out = out.rename(columns=present)

    # 2) Reorder columns (keep extras at end)
    out = _reorder_columns(out, list(spec.column_order))

    # 3) Round floats
    existing_floats = [c for c in spec.float_columns if c in out.columns]
    if existing_floats:
        out[existing_floats] = (
            out[existing_floats]
            .apply(pd.to_numeric, errors="coerce")
            .round(spec.decimals)
        )

    # 4) Sort
    chr_col, start_col, end_col = spec.sort_keys
    if all(c in out.columns for c in (chr_col, start_col, end_col)):
        out = _stable_sort_by_chr_interval(out, chr_col, start_col, end_col)

    return out


def _chrom_sort_key(chrom: str) -> tuple[int, int | str]:
    """Stable sort key: autosomes numeric first, then X/Y, then other contigs."""
    try:
        return (0, int(chrom))
    except ValueError:
        c = str(chrom)
        if c in ("X", "x"):
            return (1, 23)
        if c in ("Y", "y"):
            return (1, 24)
        return (2, c)


def _stable_sort_by_chr_interval(
    df: pd.DataFrame,
    chr_col: str,
    start_col: str,
    end_col: str,
    *,
    tmp_col: str = "chr_sort",
) -> pd.DataFrame:
    """Stable-sort by (chr, start, end) using `_chrom_sort_key`."""
    df = df.copy()
    df[tmp_col] = df[chr_col].map(_chrom_sort_key)
    df = df.sort_values(by=[tmp_col, start_col, end_col], kind="stable").drop(
        columns=[tmp_col]
    )
    return df


def _reorder_columns(df: pd.DataFrame, preferred: list[str]) -> pd.DataFrame:
    """Move preferred columns first (when present), preserving the rest."""
    preferred_present = [c for c in preferred if c in df.columns]
    return df[preferred_present + [c for c in df.columns if c not in preferred_present]]


def _flatten_agg_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Flatten MultiIndex columns produced by pandas .agg()."""
    df = df.copy()
    df.columns = [
        "_".join(col).strip("_") if isinstance(col, tuple) else col
        for col in df.columns
    ]
    return df


def _pon_abs_z(effect: float, spread: float, n_targets: float, *, min_n: int) -> float:
    """Compute abs(effect) / (spread / sqrt(n)) with defensive checks."""
    if pd.isna(effect) or pd.isna(spread) or spread <= 0:
        return np.nan
    if pd.isna(n_targets) or n_targets < min_n:
        return np.nan
    sigma_eff = spread / np.sqrt(float(n_targets))
    if sigma_eff <= 0:
        return np.nan
    return abs(effect) / sigma_eff


def _pon_direction(effect: float) -> str:
    """Map effect sign to 'gain'/'loss'/'neutral' (or '' when missing)."""
    if pd.isna(effect):
        return ""
    if effect > 0:
        return "gain"
    if effect < 0:
        return "loss"
    return "neutral"


def _pon_significance(z: float, *, noise_lt: float, borderline_lt: float) -> str:
    """Map z-score to 'noise'/'borderline'/'significant' (or '' when missing)."""
    if pd.isna(z):
        return ""
    if z < noise_lt:
        return "noise"
    if z < borderline_lt:
        return "borderline"
    return "significant"


def _pon_cnv_call_from_log2(
    *,
    is_significant: bool,
    mean_log2: float,
    gain_gt: float,
    loss_lt: float,
    non_significant_value: str,
    neutral_value: str,
) -> str:
    """Return AMPLIFICATION/DELETION/neutral/nonsignificant depending on thresholds."""
    if not is_significant or pd.isna(mean_log2):
        return non_significant_value
    if mean_log2 > gain_gt:
        return "GAIN"
    if mean_log2 < loss_lt:
        return "LOSS"
    return neutral_value


def _compute_gene_level_pon_from_bins(bins: pd.DataFrame) -> pd.DataFrame | None:
    """
    Gene-level PON stats computed from per-bin data (independent of segments/chunks).
    """
    if not {"pon_log2", "pon_spread"}.issubset(bins.columns):
        return None

    has_pon = bins["pon_spread"].notna()
    if not has_pon.any():
        return None

    gene_agg = (
        bins.loc[has_pon]
        .groupby(["chr", "gene.symbol"], as_index=False)
        .agg(
            gene_n_targets=("log2", "count"),
            gene_mean_log2=("log2", "mean"),
            pon_gene_mean_log2=("pon_log2", "mean"),
            pon_gene_mean_spread=("pon_spread", "mean"),
        )
    )

    MIN_N_FOR_PON = 2
    gene_agg["pon_gene_effect"] = (
        gene_agg["gene_mean_log2"] - gene_agg["pon_gene_mean_log2"]
    )

    gene_agg["pon_gene_z"] = gene_agg.apply(
        lambda r: _pon_abs_z(
            r["pon_gene_effect"],
            r["pon_gene_mean_spread"],
            r["gene_n_targets"],
            min_n=MIN_N_FOR_PON,
        ),
        axis=1,
    )

    gene_agg["pon_gene_direction"] = gene_agg["pon_gene_effect"].apply(_pon_direction)
    gene_agg["pon_gene_significance"] = gene_agg["pon_gene_z"].apply(
        lambda z: _pon_significance(z, noise_lt=1.5, borderline_lt=3.0)
    )

    gene_agg["pon_gene_indication"] = gene_agg.apply(
        lambda r: _pon_cnv_call_from_log2(
            is_significant=str(r.get("pon_gene_significance", "")).strip().lower()
            == "significant",
            mean_log2=r.get("gene_mean_log2", np.nan),
            gain_gt=0.08,
            loss_lt=-0.08,
            non_significant_value="",
            neutral_value="",
        ),
        axis=1,
    )

    return gene_agg[
        [
            "chr",
            "gene.symbol",
            "pon_gene_mean_log2",
            "pon_gene_mean_spread",
            "pon_gene_effect",
            "pon_gene_z",
            "pon_gene_direction",
            "pon_gene_significance",
            "pon_gene_indication",
        ]
    ]


def classify_cnv_from_total_cn_sex_aware(
    cn: float | int | None,
    chrom: str | int,
    sex: Gender,
) -> str:
    """
    Sex-aware classifier based on absolute total copy number.

    - Autosomes (1–22): baseline = 2
    - Chr X: baseline = 2 (female), 1 (male)
    - Chr Y: baseline = 1 (male), 0 (female)

    Returns: 'DELETION', 'AMPLIFICATION', or 'NEUTRAL'
    """
    if pd.isna(cn):
        return ""

    if chrom.isdigit():
        expected = 2
    elif chrom.upper() == "X":
        expected = 2 if sex == Gender.FEMALE else 1
    elif chrom.upper() == "Y":
        expected = 1 if sex == Gender.MALE else 0
    else:
        expected = 2

    cn_int = int(round(float(cn)))
    exp_int = int(round(expected))

    if cn_int < exp_int:
        return "DELETION"
    if cn_int > exp_int:
        return "AMPLIFICATION"
    return "NEUTRAL"


def add_cnv_calls_wide(
    df: pd.DataFrame,
    *,
    sex: Gender,
    chr_col: str = "chr",
    cnvkit_total_cn_col: str = "cnvkit_seg_cn",
    purecn_total_cn_col: str = "purecn_C",
    out_cnvkit_col: str = "cnvkit_cnv_call",
    out_purecn_col: str = "purecn_cnv_call",
) -> pd.DataFrame:
    out = df.copy()

    if cnvkit_total_cn_col in out.columns:
        out[out_cnvkit_col] = out.apply(
            lambda r: classify_cnv_from_total_cn_sex_aware(
                r[cnvkit_total_cn_col], r[chr_col], sex
            ),
            axis=1,
        )

    if purecn_total_cn_col in out.columns:
        out[out_purecn_col] = out.apply(
            lambda r: classify_cnv_from_total_cn_sex_aware(
                r[purecn_total_cn_col], r[chr_col], sex
            ),
            axis=1,
        )

    return out


'''
def _compute_exons_hit_for_region(
    exon_map: dict[tuple[str, str], dict[str, Any]],
    chrom: str,
    gene: str,
    region_start: float | int | None,
    region_end: float | int | None,
) -> str:
    """
    Return which exons overlap [region_start, region_end) for (chrom, gene).

    Output:
      - "" if no overlap / no exon model
      - "whole_gene" if the region fully covers the gene span
      - otherwise: "exon N" or "exon N-M" (span of hit exons)
    """
    if region_start is None or region_end is None:
        return ""
    if pd.isna(region_start) or pd.isna(region_end):
        return ""

    key = (str(chrom), str(gene))
    info = exon_map.get(key)
    if not info:
        return ""

    exons = info.get("exons") or []
    if not exons:
        return ""

    # Gene span
    gene_start = min(s for s, _e in exons)
    gene_end = max(e for _s, e in exons)

    # Whole-gene coverage
    rs = float(region_start)
    re = float(region_end)
    if rs <= float(gene_start) and re >= float(gene_end):
        return "whole_gene"

    # Overlapping exon indices (1-based)
    hit = [
        i for i, (s, e) in enumerate(exons, start=1) if float(e) > rs and float(s) < re
    ]
    if not hit:
        return ""

    lo = min(hit)
    hi = max(hit)
    return f"{lo}" if lo == hi else f"{lo}-{hi}"


def _add_exons_hit_column(
    df: pd.DataFrame,
    exon_map: dict[tuple[str, str], dict[str, Any]],
    *,
    chr_col: str = "chr",
    gene_col: str = "gene.symbol",
    start_col: str,
    end_col: str,
    out_col: str,
) -> pd.DataFrame:
    """Add/overwrite df[out_col] by computing exon overlaps vs [start_col, end_col]."""

    def _one(row: pd.Series) -> str:
        return _compute_exons_hit_for_region(
            exon_map=exon_map,
            chrom=str(row.get(chr_col, "")),
            gene=str(row.get(gene_col, "")),
            region_start=row.get(start_col),
            region_end=row.get(end_col),
        )

    df[out_col] = df.apply(_one, axis=1)
    return df
'''

############################
# CHUNK LEVEL
############################


def annotate_regions_with_overlapping_segments(
    regions_df: pd.DataFrame,
    segs_df: pd.DataFrame,
    field_map: Mapping[str, str] | None = None,
    *,
    region_chr_col: str = "chr",
    region_start_col: str = "region_start",
    region_end_col: str = "region_end",
    seg_chr_col: str = "chr",
    seg_start_col: str = "start",
    seg_end_col: str = "end",
) -> pd.DataFrame:
    """
    Generic annotator: for each region, pick the segment with max bp overlap and copy fields.

    field_map maps segs_df column -> regions_df output column name.
    seg_*_col allow segment coordinate columns to be prefixed (e.g. cnvkit_seg_start).
    """
    out = regions_df.copy()

    # Keep only fields present in segs_df
    field_map = {src: dst for src, dst in field_map.items() if src in segs_df.columns}

    # Per-chrom map for speed
    segs_by_chr: dict[str, pd.DataFrame] = {
        str(ch): df.reset_index(drop=True)
        for ch, df in segs_df.groupby(seg_chr_col, sort=False)
    }

    # Ensure destination columns exist
    for src, dst in field_map.items():
        if dst not in out.columns:
            out[dst] = pd.Series(pd.NA, index=out.index, dtype=segs_df[src].dtype)

    # Main loop
    for i, row in out.iterrows():
        ch = str(row.get(region_chr_col, ""))
        segs_chr = segs_by_chr.get(ch)
        if segs_chr is None or segs_chr.empty:
            continue

        picked = _pick_best_overlapping_segment(
            segs_chr,
            row.get(region_start_col, np.nan),
            row.get(region_end_col, np.nan),
            seg_start_col=seg_start_col,
            seg_end_col=seg_end_col,
        )
        if picked is None:
            continue

        for src, dst in field_map.items():
            val = picked.get(src, pd.NA)
            out.at[i, dst] = val

    return out


def _pick_best_overlapping_segment(
    segs_chr: pd.DataFrame,
    region_start: float,
    region_end: float,
    *,
    seg_start_col: str = "start",
    seg_end_col: str = "end",
) -> pd.Series | None:
    """Return the segment row (Series) with max bp overlap vs region, or None."""
    if segs_chr is None or segs_chr.empty:
        return None
    if pd.isna(region_start) or pd.isna(region_end):
        return None

    mask = (segs_chr[seg_end_col] > region_start) & (
        segs_chr[seg_start_col] < region_end
    )
    if not mask.any():
        return None

    ssub = segs_chr.loc[mask]
    ov0 = np.maximum(ssub[seg_start_col].to_numpy(), region_start)
    ov1 = np.minimum(ssub[seg_end_col].to_numpy(), region_end)
    best = int((ov1 - ov0).argmax())
    return ssub.iloc[best]


def annotate_regions_with_cnvkit_segments(
    regions_df: pd.DataFrame,
    cns_df: pd.DataFrame,
    prefix: str = "cnvkit_",
) -> pd.DataFrame:
    seg_start = f"{prefix}seg_start"
    seg_end = f"{prefix}seg_end"

    field_map = {c: c for c in cns_df.columns if c.startswith(prefix)}

    return annotate_regions_with_overlapping_segments(
        regions_df,
        cns_df,
        field_map=field_map,
        seg_start_col=seg_start,
        seg_end_col=seg_end,
    )


def annotate_regions_with_purecn_lohregions(
    regions_df: pd.DataFrame,
    lohregions_df: pd.DataFrame,
    prefix: str = "purecn_",
) -> pd.DataFrame:
    seg_start = f"{prefix}seg_start"
    seg_end = f"{prefix}seg_end"

    field_map = {c: c for c in lohregions_df.columns if c.startswith(prefix)}

    return annotate_regions_with_overlapping_segments(
        regions_df,
        lohregions_df,
        field_map=field_map,
        seg_start_col=seg_start,
        seg_end_col=seg_end,
    )


def merge_cnr_with_pon(
    cnr_df: pd.DataFrame,
    pon_df: pd.DataFrame | None,
    *,
    key_cols: tuple[str, str, str, str] = ("chr", "start", "end", "gene.symbol"),
    pon_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Merge exploded CNR bins with exploded PON bins on (chr,start,end,gene.symbol).

    - Drops backbone rows from both.
    - De-duplicates PON on the merge key to avoid many-to-many inflation.
    - Keeps CNR rows even if PON is missing (left join).

    If pon_cols is given, only those PON columns (plus the key) are merged in.
    """
    cnr = cnr_df.copy()
    pon = pon_df.copy() if pon_df is not None else None

    g_cnr = cnr["gene.symbol"].astype("string").str.strip()
    cnr = cnr.loc[g_cnr.ne("backbone")].copy()

    if pon is None:
        return cnr

    g_pon = pon["gene.symbol"].astype("string").str.strip()
    pon = pon.loc[g_pon.ne("backbone")].copy()

    # Choose which PON columns to bring (prevents collisions like log2/depth/etc.)
    key = list(key_cols)
    if pon_cols is not None:
        keep = [c for c in (key + pon_cols) if c in pon.columns]
        pon = pon[keep].copy()

    # Ensure PON is unique per key (critical for stable merge)
    pon = pon.drop_duplicates(subset=key)

    return cnr.merge(
        pon,
        how="left",
        on=key,
        validate="many_to_one",  # cnr may repeat; pon must be unique per key
    )


def create_gene_chunks(cnr_df: pd.DataFrame, pon_df: pd.DataFrame | None = None):
    """
    Always returns a gene/chunk table.

    If pon_df is provided:
      - run PON-based chunking (your existing logic)
    If pon_df is None/empty:
      - return one row per gene ("genelevel") without PON-derived chunking
    """

    has_pon = pon_df is not None and not pon_df.empty

    if has_pon:
        # --- your existing path ---
        bins = merge_cnr_with_pon(
            cnr_df=cnr_df,
            pon_df=pon_df,
            key_cols=("chr", "start", "end", "gene.symbol"),
            pon_cols=["pon_log2", "pon_spread"],
        )
    else:
        # --- no PON: bins are just the sample bins ---
        bins = cnr_df.copy()
        # Ensure columns exist so later code can be uniform
        bins["pon_log2"] = np.nan
        bins["pon_spread"] = np.nan

    bins = bins.sort_values(["chr", "gene.symbol", "start"], kind="stable").reset_index(
        drop=True
    )

    # ------------------------------------------------------------------
    # If NO PON: do NOT assign chunks. Just one "chunk" per gene.
    # ------------------------------------------------------------------
    if not has_pon:
        bins["chunk_id"] = "genelevel"

        agg_dict: dict[str, list[str]] = {
            "start": ["min"],
            "end": ["max"],
            "log2": ["count", "mean", "min", "max"],
            "pon_log2": ["mean"],
            "pon_spread": ["mean"],
        }

        chunks_df = bins.groupby(
            ["chr", "gene.symbol", "chunk_id"], as_index=False
        ).agg(agg_dict)
        chunks_df = _flatten_agg_columns(chunks_df).rename(
            columns={
                "start_min": "region_start",
                "end_max": "region_end",
                "log2_count": "n_targets",
                "log2_mean": "mean_log2",
                "log2_min": "min_log2",
                "log2_max": "max_log2",
                "pon_log2_mean": "pon_mean_log2",
                "pon_spread_mean": "pon_mean_spread",
            }
        )
        chunks_df["n.targets"] = chunks_df["n_targets"]

        # PON-derived outputs: keep columns but make them empty/neutral
        chunks_df["pon_chunk_effect"] = np.nan
        chunks_df["pon_chunk_z"] = np.nan
        chunks_df["pon_chunk_significance"] = ""
        chunks_df["pon_chunk_indication"] = "NEUTRAL"

        return chunks_df

    # CREATE INITIAL CHUNKS
    MIN_GENE_TARGETS = 8
    MIN_RUN_BINS = 4
    Z_BIN_THRESH = 1.5
    Z_RUN_THRESH = 3.0
    SMOOTH_WINDOW = 3

    bins["chunk_id"] = "no_chunk"
    next_id = 0

    def _run_z_for_indices(ix: list[int]) -> float:
        sub = bins.loc[ix]
        eff = sub["log2"].to_numpy() - sub["pon_log2"].fillna(0.0).to_numpy()
        if eff.size == 0:
            return np.nan
        mean_eff = float(np.nanmean(eff))
        if not np.isfinite(mean_eff):
            return np.nan

        sigma = sub["pon_spread"].fillna(0.0).to_numpy()
        sigma_safe = np.where(sigma <= 0, 1e-3, sigma)
        sigma_mean = float(np.nanmean(sigma_safe)) if sigma_safe.size else 1e-3

        n = len(eff)
        sigma_eff = sigma_mean / np.sqrt(float(n))
        if sigma_eff <= 0:
            return np.nan
        return abs(mean_eff) / sigma_eff

    def _assign_gene_chunks(df_gene: pd.DataFrame) -> None:
        nonlocal next_id

        if df_gene.shape[0] < MIN_GENE_TARGETS:
            return
        if df_gene["pon_spread"].isna().all():
            return

        idxs = df_gene.index.to_numpy()
        eff = df_gene["log2"].to_numpy() - df_gene["pon_log2"].fillna(0.0).to_numpy()
        sigma = df_gene["pon_spread"].fillna(0.0).to_numpy()
        sigma_safe = np.where(sigma <= 0, 1e-3, sigma)
        z_raw = eff / sigma_safe

        z_smooth = (
            pd.Series(z_raw, index=df_gene.index)
            .rolling(window=SMOOTH_WINDOW, center=True, min_periods=1)
            .median()
            .to_numpy()
        )

        run: list[int] = []
        run_sign: int | None = None

        def _finalize(ix: list[int]) -> None:
            nonlocal next_id
            if len(ix) < MIN_RUN_BINS:
                return
            rz = _run_z_for_indices(ix)
            if not np.isfinite(rz) or rz < Z_RUN_THRESH:
                return
            label = f"genechunk_{next_id}"
            next_id += 1
            bins.loc[ix, "chunk_id"] = label

        for i, z in zip(idxs, z_smooth):
            if not np.isfinite(z) or abs(z) < Z_BIN_THRESH:
                if run:
                    _finalize(run)
                    run = []
                    run_sign = None
                continue

            sgn = 1 if z > 0 else -1
            if run_sign is None or sgn == run_sign:
                run.append(i)
                run_sign = sgn
            else:
                _finalize(run)
                run = [i]
                run_sign = sgn

        if run:
            _finalize(run)

    for (_ch, _g), df_gene in bins.groupby(["chr", "gene.symbol"], sort=False):
        _assign_gene_chunks(df_gene)

    # MERGE ADJACENT CHUNKS

    MAX_BRIDGE_BINS = 4
    BRIDGE_DELTA = 0.12
    SMALL_SEG_N = 3
    MERGE_DELTA = 0.10

    def _gene_runs(df_gene: pd.DataFrame) -> list[dict]:
        eff = df_gene["log2"].to_numpy() - df_gene["pon_log2"].fillna(0.0).to_numpy()
        labels = df_gene["chunk_id"].tolist()

        runs: list[dict] = []
        pos = [0]
        cur = labels[0]
        for i in range(1, len(labels)):
            if labels[i] == cur:
                pos.append(i)
            else:
                ix = df_gene.index[pos].tolist()
                mean_eff = float(np.nanmean(eff[pos])) if pos else np.nan
                runs.append({"positions": pos, "indices": ix, "mean_eff": mean_eff})
                cur = labels[i]
                pos = [i]
        if pos:
            ix = df_gene.index[pos].tolist()
            mean_eff = float(np.nanmean(eff[pos])) if pos else np.nan
            runs.append({"positions": pos, "indices": ix, "mean_eff": mean_eff})
        return runs, eff

    def _cleanup_gene_chunks(df_gene: pd.DataFrame) -> None:
        if df_gene.empty or df_gene.shape[0] == 1:
            return

        runs, eff_all = _gene_runs(df_gene)
        if len(runs) <= 1:
            return

        new_runs: list[dict] = []
        i = 0
        while i < len(runs):
            if i <= len(runs) - 3:
                A, B, C = runs[i], runs[i + 1], runs[i + 2]
                if (
                    len(B["indices"]) <= MAX_BRIDGE_BINS
                    and np.isfinite(A["mean_eff"])
                    and np.isfinite(C["mean_eff"])
                    and abs(A["mean_eff"] - C["mean_eff"]) <= BRIDGE_DELTA
                ):
                    merged_pos = A["positions"] + B["positions"] + C["positions"]
                    merged_ix = A["indices"] + B["indices"] + C["indices"]
                    merged_mean = (
                        float(np.nanmean(eff_all[merged_pos])) if merged_pos else np.nan
                    )
                    new_runs.append(
                        {
                            "positions": merged_pos,
                            "indices": merged_ix,
                            "mean_eff": merged_mean,
                        }
                    )
                    i += 3
                    continue

            r = runs[i]
            if new_runs:
                last = new_runs[-1]
                if (
                    np.isfinite(r["mean_eff"])
                    and np.isfinite(last["mean_eff"])
                    and abs(r["mean_eff"] - last["mean_eff"]) <= MERGE_DELTA
                    and (
                        len(r["indices"]) <= SMALL_SEG_N
                        or len(last["indices"]) <= SMALL_SEG_N
                    )
                ):
                    merged_pos = last["positions"] + r["positions"]
                    merged_ix = last["indices"] + r["indices"]
                    merged_mean = (
                        float(np.nanmean(eff_all[merged_pos])) if merged_pos else np.nan
                    )
                    new_runs[-1] = {
                        "positions": merged_pos,
                        "indices": merged_ix,
                        "mean_eff": merged_mean,
                    }
                else:
                    new_runs.append(r)
            else:
                new_runs.append(r)

            i += 1

        for run_idx, run in enumerate(new_runs):
            bins.loc[run["indices"], "chunk_id"] = f"genechunk_clean_{run_idx}"

    for (_ch, _g), df_gene in bins.groupby(["chr", "gene.symbol"], sort=False):
        _cleanup_gene_chunks(df_gene)

    bins = bins.sort_values(["chr", "gene.symbol", "start"], kind="stable").reset_index(
        drop=True
    )

    # Collapse to gene × chunk
    agg_dict: dict[str, list[str]] = {
        "start": ["min"],
        "end": ["max"],
        "log2": ["count", "mean", "min", "max"],
        "pon_log2": ["mean"],
        "pon_spread": ["mean"],
    }

    chunks_df = bins.groupby(["chr", "gene.symbol", "chunk_id"], as_index=False).agg(
        agg_dict
    )
    chunks_df = _flatten_agg_columns(chunks_df).rename(
        columns={
            "start_min": "region_start",
            "end_max": "region_end",
            "log2_count": "n_targets",
            "log2_mean": "mean_log2",
            "log2_min": "min_log2",
            "log2_max": "max_log2",
            "pon_log2_mean": "pon_mean_log2",
            "pon_spread_mean": "pon_mean_spread",
        }
    )
    chunks_df["n.targets"] = chunks_df["n_targets"]

    # PON deviation per chunk
    MIN_CHUNK_TARGETS_FOR_CALL = 5

    chunks_df["pon_chunk_effect"] = chunks_df["mean_log2"] - chunks_df["pon_mean_log2"]

    # Compute z with min_n=5 so small chunks naturally become NaN/0 depending on _pon_abs_z
    chunks_df["pon_chunk_z"] = chunks_df.apply(
        lambda r: _pon_abs_z(
            r["pon_chunk_effect"],
            r["pon_mean_spread"],
            r.get("n_targets", r.get("n.targets", np.nan)),
            min_n=MIN_CHUNK_TARGETS_FOR_CALL,
        ),
        axis=1,
    )

    chunks_df["pon_chunk_direction"] = chunks_df["pon_chunk_effect"].apply(
        _pon_direction
    )

    chunks_df["pon_chunk_significance"] = chunks_df["pon_chunk_z"].apply(
        lambda z: _pon_significance(z, noise_lt=2.0, borderline_lt=5.0)
    )

    # Hard gate: never allow calls for chunks with too few targets
    too_small = (
        chunks_df["n_targets"].fillna(0).astype(int) < MIN_CHUNK_TARGETS_FOR_CALL
    )
    chunks_df.loc[
        too_small, "pon_chunk_significance"
    ] = "noise"  # or "not_significant" if that's what you use
    chunks_df.loc[too_small, "pon_chunk_indication"] = "NEUTRAL"

    # Only compute indication for eligible chunks
    eligible = ~too_small
    chunks_df.loc[eligible, "pon_chunk_indication"] = chunks_df.loc[eligible].apply(
        lambda r: _pon_cnv_call_from_log2(
            is_significant=str(r.get("pon_chunk_significance", "")).strip().lower()
            == "significant",
            mean_log2=r.get("mean_log2", np.nan),
            gain_gt=0.07,
            loss_lt=-0.07,
            non_significant_value="NEUTRAL",
            neutral_value="NEUTRAL",
        ),
        axis=1,
    )

    return chunks_df


def build_gene_chunk_table(
    cnr_df: pd.DataFrame,
    cns_df: pd.DataFrame,
    sex: Gender,
    pon_df: pd.DataFrame | None = None,
    cancer_genes: set[str] | None = None,
    loh_regions_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Build per-gene table always.
    If pon_df exists -> chunked
    Else -> one-row-per-gene 'genelevel' chunks.
    """

    cancer_genes = cancer_genes or set()

    # Always create something
    chunks_df = create_gene_chunks(cnr_df=cnr_df, pon_df=pon_df)

    # Annotate with CNVkit segments
    chunks_df = annotate_regions_with_cnvkit_segments(chunks_df, cns_df)

    # Attach PureCN LOHregions (optional)
    if loh_regions_df is not None and not loh_regions_df.empty:
        chunks_df = annotate_regions_with_purecn_lohregions(chunks_df, loh_regions_df)

    chunks_df["is_cancer_gene"] = chunks_df["gene.symbol"].isin(cancer_genes)

    # CNV calls
    chunks_df = add_cnv_calls_wide(chunks_df, sex=sex)

    # Drop columns if present
    drop_cols = [
        c
        for c in ["chunk_id", "n_targets", "pon_chunk_direction"]
        if c in chunks_df.columns
    ]
    if drop_cols:
        chunks_df = chunks_df.drop(columns=drop_cols)

    return finalize_table(chunks_df, GENE_TABLE_SPEC)


############################
# SEGMENT LEVEL
############################
def add_overlapping_genes_from_bins(
    seg_df: pd.DataFrame,
    cnr_df: pd.DataFrame,
    *,
    genes_col: str = "gene.symbol",
    out_targets_col: str = "n.targets",
    min_targets: int = 2,
    cancer_genes: set[str] | None = None,
    drop_genes: set[str] | None = None,
    # which bins count as "targets" for n.targets
    target_drop_genes: set[str] | None = None,
) -> pd.DataFrame:
    """
    For each segment row (chr/start/end), find overlapping CNVkit bins in cnr_df and:
      1) add comma-separated gene list in `genes_col` (genes w/ >= min_targets bins)
      2) add `out_targets_col` = number of UNIQUE overlapping *target* bins

    Notes:
      - cnr_df is typically exploded by gene.symbol, so `out_targets_col` counts unique bins
        by (chr,start,end) to avoid overcounting.
      - `drop_genes` affects the gene list.
      - `target_drop_genes` affects counting targets (defaults to drop_genes if not provided).
    """
    out = seg_df.copy()

    # Normalize segment types
    out["chr"] = out["chr"].astype("string")
    out["start"] = pd.to_numeric(out["start"], errors="coerce").astype("Int64")
    out["end"] = pd.to_numeric(out["end"], errors="coerce").astype("Int64")

    bins = cnr_df.copy()
    bins["chr"] = bins["chr"].astype("string")
    bins["start"] = pd.to_numeric(bins["start"], errors="coerce").astype("Int64")
    bins["end"] = pd.to_numeric(bins["end"], errors="coerce").astype("Int64")
    bins["gene.symbol"] = bins["gene.symbol"].astype("string")

    # Drop unusable rows early
    out = out.dropna(subset=["chr", "start", "end"]).copy()
    bins = bins.dropna(subset=["chr", "start", "end", "gene.symbol"]).copy()

    # normalize drop sets (lowercase)
    drop_set: set[str] = set()
    if drop_genes:
        drop_set = {str(g).strip().lower() for g in drop_genes if str(g).strip()}

    target_drop_set: set[str] = set(drop_set)
    if target_drop_genes is not None:
        target_drop_set = {
            str(g).strip().lower() for g in target_drop_genes if str(g).strip()
        }

    # Prepare output columns
    out[genes_col] = ""
    out[out_targets_col] = 0

    # Work per chromosome
    for chrom, seg_g in out.groupby("chr", sort=False):
        bins_g = bins[bins["chr"] == chrom].copy()
        if bins_g.empty:
            continue

        bins_g = bins_g.sort_values("start", kind="stable").reset_index(drop=True)

        b_start = bins_g["start"].to_numpy(dtype=np.int64, copy=False)
        b_end = bins_g["end"].to_numpy(dtype=np.int64, copy=False)
        b_gene = bins_g["gene.symbol"].to_numpy(dtype=object, copy=False)
        # for unique-bin counting (avoid exploded overcount)
        b_key_start = b_start
        b_key_end = b_end

        seg_idx = seg_g.index.to_numpy()
        s_start = seg_g["start"].to_numpy(dtype=np.int64, copy=False)
        s_end = seg_g["end"].to_numpy(dtype=np.int64, copy=False)

        for i, row_idx in enumerate(seg_idx):
            lo = s_start[i]
            hi = s_end[i]
            if hi <= lo:
                continue

            cut = np.searchsorted(b_start, hi, side="left")
            if cut == 0:
                continue

            cand_ends = b_end[:cut]
            mask = cand_ends > lo
            if not np.any(mask):
                continue

            # ----------------------------
            # (A) n.targets = unique bins overlapping, excluding backbone etc
            # ----------------------------
            genes_overlap = (
                pd.Series(b_gene[:cut][mask], dtype="string").fillna("").astype(str)
            )
            if target_drop_set:
                keep_targets = (
                    ~genes_overlap.str.strip().str.lower().isin(target_drop_set)
                )
            else:
                keep_targets = genes_overlap.ne("")

            if keep_targets.any():
                # (A) n.targets = number of UNIQUE bins overlapping this segment (no exclusions)
                starts = b_key_start[:cut][mask]
                ends = b_key_end[:cut][mask]

                # unique (start,end) pairs => unique bins
                uniq = np.unique(np.stack([starts, ends], axis=1), axis=0)
                out.at[row_idx, out_targets_col] = int(uniq.shape[0])

            # ----------------------------
            # (B) gene list (>= min_targets bins per gene) with optional cancer restriction
            # ----------------------------
            genes_for_list = genes_overlap

            if drop_set:
                genes_for_list = genes_for_list[
                    ~genes_for_list.str.strip().str.lower().isin(drop_set)
                ]

            if cancer_genes:
                cancer_set = {str(g).strip() for g in cancer_genes if str(g).strip()}
                if cancer_set:
                    genes_for_list = genes_for_list[genes_for_list.isin(cancer_set)]

            if genes_for_list.empty:
                continue

            vc = genes_for_list.value_counts(dropna=True)
            vc = vc[vc >= min_targets]
            if vc.empty:
                continue

            gene_list = sorted(vc.index.astype(str).tolist())
            out.at[row_idx, genes_col] = ",".join(gene_list)

    return out


def annotate_segments_with_cytoband(
    df_segments: pd.DataFrame,
    cyto: pd.DataFrame,
    *,
    seg_chr_col: str = "chr",
    seg_start_col: str = "start",
    seg_end_col: str = "end",
    cyto_chr_col: str = "chr",
    cyto_start_col: str = "start_int",
    cyto_end_col: str = "end_int",
    cyto_name_col: str = "name",
) -> pd.DataFrame:
    """
    Add 'cytoband' based on overlap between segment (chr,start,end) and cytobands.
    """
    df = df_segments.copy()
    df["cytoband"] = pd.Series(pd.NA, index=df.index, dtype="string")

    # make sure types are sane
    df[seg_chr_col] = df[seg_chr_col].astype("string")
    df[seg_start_col] = pd.to_numeric(df[seg_start_col], errors="coerce").astype(
        "Int64"
    )
    df[seg_end_col] = pd.to_numeric(df[seg_end_col], errors="coerce").astype("Int64")

    cy = cyto.copy()
    cy[cyto_chr_col] = cy[cyto_chr_col].astype("string")
    cy[cyto_start_col] = pd.to_numeric(cy[cyto_start_col], errors="coerce").astype(
        "Int64"
    )
    cy[cyto_end_col] = pd.to_numeric(cy[cyto_end_col], errors="coerce").astype("Int64")

    for chrom, seg_chr in df.groupby(seg_chr_col, sort=False):
        cyto_chr = cy[cy[cyto_chr_col] == chrom]
        if cyto_chr.empty:
            continue

        cyto_chr = cyto_chr.sort_values(cyto_start_col, kind="stable")
        band_starts = cyto_chr[cyto_start_col].to_numpy()
        band_ends = cyto_chr[cyto_end_col].to_numpy()
        band_names = cyto_chr[cyto_name_col].to_numpy()

        # iterate ONLY rows on this chromosome
        for idx, row in seg_chr.iterrows():
            if pd.isna(row[seg_start_col]) or pd.isna(row[seg_end_col]):
                continue

            s_start = int(row[seg_start_col])
            s_end = int(row[seg_end_col])

            mask = (band_starts <= s_end) & (band_ends >= s_start)
            if not mask.any():
                continue

            names = band_names[mask]
            if len(names) == 1:
                label = f"{chrom}{names[0]}"
            else:
                label = f"{chrom}{names[0]}-{chrom}{names[-1]}"

            df.at[idx, "cytoband"] = label

    return df


def build_segment_table(
    cnr_df: pd.DataFrame,
    cns_df: pd.DataFrame,
    cytoband_df: pd.DataFrame,
    sex: Gender,
    cancer_genes: set[str],
    is_exome: bool,
    loh_regions_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """
    Build a unified segment table by stacking CNVkit and PureCN segments.

    Output core columns:
      - chr, start, end, caller

    CNVkit-specific (if present in cns_df):
      - cnvkit_seg_cn, cnvkit_seg_cn1, cnvkit_seg_cn2, cnvkit_seg_depth
      - cnvkit_seg_log2, cnvkit_seg_raw_log2

    PureCN-specific (if present in loh_regions_df):
      - purecn_C, purecn_M, purecn_M_flagged, purecn_seg_mean_log2, purecn_num_snps
      - purecn_maf_observed, purecn_type
      - purecn_is_loh (derived: purecn_type contains 'LOH')

    Notes:
      - Duplicated chr/start/end rows are allowed (CNVkit and PureCN can overlap).
      - This function *stacks* callers (long format). It does not attempt interval intersection.
    """
    out_parts: list[pd.DataFrame] = []

    # ---- CNVkit ----
    cnv = cns_df.copy()
    # normalize coords -> shared names
    cnv = cnv.rename(columns={"cnvkit_seg_start": "start"})
    cnv = cnv.rename(columns={"cnvkit_seg_end": "end"})
    cnv = cnv.rename(columns={"cnvkit_seg_log2": "cnvkit_adjusted_log2"})
    cnv = cnv.rename(columns={"cnvkit_seg_raw_log2": "log2"})
    cnv = cnv.rename(columns={"cnvkit_seg_baf": "baf_maf"})

    cnv["caller"] = "CNVkit"
    out_parts.append(cnv)

    # ---- PureCN ----
    if loh_regions_df is not None and not loh_regions_df.empty:
        pc = loh_regions_df.copy()

        pc = pc.rename(columns={"purecn_seg_start": "start"})
        pc = pc.rename(columns={"purecn_seg_end": "end"})
        pc = pc.rename(columns={"purecn_seg_mean_log2": "log2"})
        pc = pc.rename(columns={"purecn_maf_observed": "baf_maf"})

        pc["caller"] = "PureCN"

        out_parts.append(pc)

    segments = pd.concat(out_parts, axis=0, ignore_index=True, sort=False)
    segments = segments.sort_values(
        ["chr", "start", "end", "caller"], kind="stable"
    ).reset_index(drop=True)

    # Add overlapping genes (limit to cancer genes only if provided)
    segments = add_overlapping_genes_from_bins(
        segments,
        cnr_df,
        genes_col="gene.symbol",
        min_targets=2,
        cancer_genes=cancer_genes if is_exome else None,
        drop_genes={"backbone"},
    )
    segments = add_cnv_calls_wide(
        segments, sex=sex
    )  # makes cnvkit_cnv_call and purecn_cnv_call

    # Then unify
    segments["cnv_call"] = pd.NA
    mask_cnv = segments["caller"].astype("string").str.upper().eq("CNVKIT")
    mask_pc = segments["caller"].astype("string").str.upper().eq("PURECN")

    segments.loc[mask_cnv, "cnv_call"] = segments.loc[mask_cnv, "cnvkit_cnv_call"]
    segments.loc[mask_pc, "cnv_call"] = segments.loc[mask_pc, "purecn_cnv_call"]

    segments = annotate_segments_with_cytoband(segments, cytoband_df)
    segments["segment_size"] = segments["end"] - segments["start"]

    segments = segments.drop(columns=["cnvkit_cnv_call", "purecn_cnv_call"])
    return finalize_table(segments, SEGMENT_TABLE_SPEC)
