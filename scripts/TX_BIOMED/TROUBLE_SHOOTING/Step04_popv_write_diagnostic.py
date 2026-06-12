#!/usr/bin/env python3
"""
Step04_popv_write_diagnostic.py

Targeted, fast (100-cell) reproduction of the h5ad write failure that hits
right after the popV merge in Step04_Cell_Type_Annotation.py:

    TypeError: Can't implicitly convert non-string objects to strings

Goal: decide between candidate fixes with evidence, not guesswork.

ABSOLUTE PATHS
--------------
This script can be launched from anywhere (e.g. the TROUBLE_SHOOTING dir).
All paths are resolved to absolute up front and printed before anything runs.
By default PROJECT_ROOT and SCRIPTS_DIR are auto-detected from this file's
location; if that guesses wrong, hardcode the two constants in the CONFIG
block below. Every path can also be overridden on the command line.

What it does
------------
1. Load the merged QC h5ad (raw counts).
2. Subsample to ~100 query cells, stratified by puck.
3. Run popV inference exactly as Step04 does (same repo, batch key, gene key).
4. Build the merged query obs (same merge pattern as Step04).
5. Run THREE write tiers and report PASS/FAIL + round-trip check for each:

   Tier 0  RAW         : write all popV columns (reproduce the bug)
   Tier 1  WHITELIST   : keep only popv_* + _predict_cells + over_clustering
                         (Jake's minimal chunk, verbatim)
   Tier 2  WHITELIST+   : Tier 1, plus a NARROW coercion that only touches
                         object / numpy.str_ columns and EMPTY categoricals.
                         Floats and int-categoricals are left untouched.

It does NOT modify Step04 or write any production outputs. Everything goes
to a scratch dir (default: <this script's dir>/write_diagnostic).

Run (on Titan, sc_pre or NETWORK env):
    python Step04_popv_write_diagnostic.py
    python Step04_popv_write_diagnostic.py --merged /abs/path/all_pucks_merged_QC.h5ad
    python Step04_popv_write_diagnostic.py --n-cells 100
    python Step04_popv_write_diagnostic.py --scratch /abs/path/scratch

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import argparse
import traceback
import numpy as np
import pandas as pd
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# CONFIG  -- set these to absolute paths to skip auto-detection (None = auto)
# =============================================================================
PROJECT_ROOT = None   # e.g. "/master/jlehle/WORKING/<spatial_project>"
SCRIPTS_DIR  = None   # dir containing spatial_config.py / Step04_*.py

# popV settings, identical to Step04.run_popv_annotation
HF_REPO = "popV/tabula_sapiens_All_Cells"
QUERY_BATCH_KEY = "puck_id"
GENE_SYMBOL_KEY = "feature_name"

# Jake's whitelist, verbatim
WHITELIST_EXTRA = ['_predict_cells', 'over_clustering']

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


# -----------------------------------------------------------------------------
# absolute-path resolution
# -----------------------------------------------------------------------------
def _abspath(p):
    return os.path.abspath(os.path.expanduser(p)) if p else p


def _autodetect_project_root():
    """Walk up from this file looking for data/outputs/03_qc."""
    d = THIS_DIR
    for _ in range(6):
        if os.path.isdir(os.path.join(d, "data", "outputs", "03_qc")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    # fall back to the parent of this script's dir (TROUBLE_SHOOTING -> root)
    return os.path.dirname(THIS_DIR)


def _autodetect_scripts_dir(root):
    """Find the dir that actually contains spatial_config.py."""
    candidates = [
        os.path.join(root, "scripts"),
        THIS_DIR,
        root,
        os.path.join(os.path.dirname(THIS_DIR), "scripts"),
    ]
    for c in candidates:
        if os.path.isfile(os.path.join(c, "spatial_config.py")):
            return c
    return candidates[0]


def resolve_config(args):
    """Resolve every path to an absolute path and return a dict."""
    root = _abspath(PROJECT_ROOT) if PROJECT_ROOT else _autodetect_project_root()
    scripts = _abspath(SCRIPTS_DIR) if SCRIPTS_DIR else _autodetect_scripts_dir(root)

    # import spatial_config from the resolved (absolute) scripts dir
    dir_03_qc = dir_04_annotation = None
    if os.path.isfile(os.path.join(scripts, "spatial_config.py")):
        sys.path.insert(0, scripts)
        try:
            import spatial_config as scfg
            d3 = getattr(scfg, "DIR_03_QC", None)
            d4 = getattr(scfg, "DIR_04_ANNOTATION", None)
            # spatial_config dirs may be relative; anchor them to root
            dir_03_qc = _abspath(d3 if os.path.isabs(str(d3)) else os.path.join(root, str(d3))) if d3 else None
            dir_04_annotation = _abspath(d4 if os.path.isabs(str(d4)) else os.path.join(root, str(d4))) if d4 else None
        except Exception as e:
            print(f"[WARN] could not import spatial_config from {scripts}: {e}")

    if dir_03_qc is None:
        dir_03_qc = _abspath(os.path.join(root, "data", "outputs", "03_qc"))
    if dir_04_annotation is None:
        dir_04_annotation = _abspath(os.path.join(root, "data", "outputs", "04_annotation"))

    merged = _abspath(args.merged) if args.merged else os.path.join(dir_03_qc, "all_pucks_merged_QC.h5ad")
    cache = _abspath(args.cache) if args.cache else os.path.join(dir_04_annotation, "popv_cache")
    scratch = _abspath(args.scratch) if args.scratch else os.path.join(THIS_DIR, "write_diagnostic")

    return {
        "PROJECT_ROOT": root,
        "SCRIPTS_DIR": scripts,
        "DIR_03_QC": dir_03_qc,
        "DIR_04_ANNOTATION": dir_04_annotation,
        "MERGED_H5AD": merged,
        "CACHE_DIR": cache,
        "SCRATCH_DIR": scratch,
    }


def banner(msg):
    print("\n" + "=" * 70)
    print(msg)
    print("=" * 70)


# -----------------------------------------------------------------------------
# obs audit: show exactly what each column looks like to the h5ad writer
# -----------------------------------------------------------------------------
def audit_obs(obs, label):
    print(f"\n[AUDIT] {label}: {obs.shape[1]} columns")
    flagged = []
    for col in obs.columns:
        s = obs[col]
        dt = str(s.dtype)
        if isinstance(s.dtype, pd.CategoricalDtype):
            cats = list(s.cat.categories)
            cat_types = {type(c).__name__ for c in cats}
            note = ""
            if len(cats) == 0:
                note = "  <-- EMPTY categorical"
                flagged.append(col)
            elif not cat_types <= {"str"}:
                note = f"  <-- non-str categories {cat_types}"
                flagged.append(col)
            print(f"    {col:<48} category  n_cat={len(cats):<4} types={cat_types}{note}")
        elif dt == "object":
            nonnull = s.dropna()
            py_types = {type(v).__name__ for v in nonnull.head(50)}
            note = ""
            if len(nonnull) == 0:
                note = "  <-- all-NaN object"
                flagged.append(col)
            elif not py_types <= {"str"}:
                note = f"  <-- non-str objects {py_types}"
                flagged.append(col)
            print(f"    {col:<48} object    py_types={py_types or set()}{note}")
        elif dt == "boolean" and s.isna().all():
            print(f"    {col:<48} boolean   <-- all-NA nullable bool")
            flagged.append(col)
    if flagged:
        print(f"  [AUDIT] potentially-unwritable columns: {flagged}")
    else:
        print(f"  [AUDIT] no obviously unwritable columns detected")
    return flagged


# -----------------------------------------------------------------------------
# write + round-trip test
# -----------------------------------------------------------------------------
def try_write_and_read(adata, path, check_col="popv_prediction"):
    """Returns (ok: bool, message: str)."""
    try:
        adata.write_h5ad(path)
    except Exception as e:
        return False, f"WRITE FAILED: {type(e).__name__}: {e}"

    try:
        rb = sc.read_h5ad(path)
    except Exception as e:
        return False, f"WRITE ok but READ FAILED: {type(e).__name__}: {e}"

    if rb.n_obs != adata.n_obs:
        return False, f"round-trip cell count mismatch {adata.n_obs} -> {rb.n_obs}"

    if check_col in adata.obs.columns and check_col in rb.obs.columns:
        before = adata.obs[check_col].astype(str).values
        after = rb.obs[check_col].astype(str).values
        if not np.array_equal(before, after):
            n_diff = int((before != after).sum())
            return False, f"round-trip '{check_col}' values changed ({n_diff} cells)"
    return True, f"PASS (write + read-back ok, n_obs={rb.n_obs})"


# -----------------------------------------------------------------------------
# the candidate fixes
# -----------------------------------------------------------------------------
def make_whitelist_obs(annotated_obs, original_columns):
    """Tier 1: Jake's verbatim whitelist."""
    popv_columns = [c for c in annotated_obs.columns
                    if c.startswith('popv_') or c in WHITELIST_EXTRA]
    return annotated_obs[popv_columns].copy(), popv_columns


def narrow_coerce(obs):
    """
    Tier 2: minimal coercion. Only fixes the shapes the writer chokes on,
    and leaves everything else exactly as-is.
      - object / numpy.str_ columns  -> fillna('') then python str
      - EMPTY categoricals           -> str
      - non-str categoricals (int)   -> LEFT ALONE (kept numeric/categorical)
      - floats, ints, bools          -> LEFT ALONE
    """
    obs = obs.copy()
    touched = []
    for col in obs.columns:
        s = obs[col]
        if s.dtype == object:
            obs[col] = s.where(s.notna(), '').astype(str)
            touched.append((col, "object->str"))
        elif isinstance(s.dtype, pd.CategoricalDtype) and len(s.cat.categories) == 0:
            obs[col] = s.astype(str)
            touched.append((col, "empty-cat->str"))
    return obs, touched


# -----------------------------------------------------------------------------
# popV (mirror of Step04.run_popv_annotation, trimmed)
# -----------------------------------------------------------------------------
def run_popv(adata, cache_dir):
    import popv
    if GENE_SYMBOL_KEY not in adata.var.columns:
        adata.var[GENE_SYMBOL_KEY] = adata.var_names
    original_columns = list(adata.obs.columns)

    os.makedirs(cache_dir, exist_ok=True)
    print(f"  pulling {HF_REPO} (cache: {cache_dir}) ...")
    hmo = popv.hub.HubModel.pull_from_huggingface_hub(HF_REPO, cache_dir=cache_dir)

    print(f"  annotating {adata.n_obs} query cells ...")
    annotated = hmo.annotate_data(
        adata,
        query_batch_key=QUERY_BATCH_KEY,
        prediction_mode="inference",
        gene_symbols=GENE_SYMBOL_KEY,
    )
    if '_dataset' in annotated.obs.columns:
        annotated = annotated[annotated.obs['_dataset'] == 'query'].copy()
    return annotated, original_columns


def merge_back(orig_adata, popv_obs):
    """Same index-merge as Step04."""
    merged = pd.merge(orig_adata.obs, popv_obs,
                      left_index=True, right_index=True, how='inner')
    out = orig_adata[orig_adata.obs_names.isin(merged.index)].copy()
    merged = merged.loc[out.obs_names]
    out.obs = merged
    return out


# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--merged", default=None, help="absolute path to all_pucks_merged_QC.h5ad")
    ap.add_argument("--cache", default=None, help="absolute path to popv_cache dir")
    ap.add_argument("--scratch", default=None, help="absolute path for diagnostic outputs")
    ap.add_argument("--n-cells", type=int, default=100)
    args = ap.parse_args()

    banner("STEP04 popV WRITE DIAGNOSTIC (100-cell subset)")
    cfg = resolve_config(args)

    print("\n[RESOLVED CONFIG] (all absolute)")
    for k, v in cfg.items():
        exists = os.path.exists(v)
        tag = "ok " if exists else "MISSING"
        print(f"  {k:<18} [{tag}] {v}")

    if not os.path.exists(cfg["MERGED_H5AD"]):
        print(f"\nFATAL: merged h5ad not found at the resolved path above.")
        print(f"       Pass --merged /abs/path, or set PROJECT_ROOT in the CONFIG block.")
        sys.exit(1)

    scratch = cfg["SCRATCH_DIR"]
    cache_dir = cfg["CACHE_DIR"]
    os.makedirs(scratch, exist_ok=True)

    # --- load + subsample --------------------------------------------------
    adata = sc.read_h5ad(cfg["MERGED_H5AD"])
    print(f"\n  loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    np.random.seed(getattr(sys.modules.get("spatial_config", object()), "RANDOM_SEED", 0))
    if QUERY_BATCH_KEY in adata.obs.columns:
        pucks = list(adata.obs[QUERY_BATCH_KEY].unique())
        per = max(1, args.n_cells // max(1, len(pucks)))
        idx = []
        for p in pucks:
            pool = adata.obs_names[adata.obs[QUERY_BATCH_KEY] == p].tolist()
            take = min(per, len(pool))
            idx += list(np.random.choice(pool, size=take, replace=False))
        idx = idx[:args.n_cells]
    else:
        idx = list(np.random.choice(adata.obs_names,
                                    size=min(args.n_cells, adata.n_obs), replace=False))
    sub = adata[idx].copy()
    pv = (sub.obs[QUERY_BATCH_KEY].value_counts().to_dict()
          if QUERY_BATCH_KEY in sub.obs else 'n/a')
    print(f"  subsampled to {sub.n_obs} cells (pucks: {pv})")

    # --- popV --------------------------------------------------------------
    banner("RUNNING popV ON SUBSET")
    annotated, original_columns = run_popv(sub, cache_dir)
    n_new = len([c for c in annotated.obs.columns if c not in original_columns])
    print(f"  annotated obs columns: {annotated.obs.shape[1]} ({n_new} new)")

    # ===================================================================
    # TIER 0 : RAW (reproduce the bug)
    # ===================================================================
    banner("TIER 0 : RAW (all popV columns, current Step04 behavior)")
    popv_obs_all = annotated.obs.drop(columns=original_columns, errors='ignore')
    raw_merged = merge_back(sub, popv_obs_all)
    audit_obs(raw_merged.obs, "Tier 0 obs (raw)")
    ok0, msg0 = try_write_and_read(raw_merged, os.path.join(scratch, "tier0_raw.h5ad"))
    print(f"\n  >>> TIER 0 RESULT: {msg0}")

    # ===================================================================
    # TIER 1 : WHITELIST (Jake's verbatim chunk)
    # ===================================================================
    banner("TIER 1 : WHITELIST (popv_* + _predict_cells + over_clustering)")
    wl_obs, kept = make_whitelist_obs(annotated.obs, original_columns)
    print(f"  kept {len(kept)} columns")
    wl_merged = merge_back(sub, wl_obs)
    flagged1 = audit_obs(wl_merged.obs, "Tier 1 obs (whitelist)")
    ok1, msg1 = try_write_and_read(wl_merged, os.path.join(scratch, "tier1_whitelist.h5ad"))
    print(f"\n  >>> TIER 1 RESULT: {msg1}")
    if flagged1:
        print(f"  >>> NOTE: whitelist still carries flagged columns: {flagged1}")
        print(f"           (if Tier 1 failed, these are the prime suspects)")

    # ===================================================================
    # TIER 2 : WHITELIST + NARROW COERCION
    # ===================================================================
    banner("TIER 2 : WHITELIST + narrow coercion (object/str_ + empty-cat only)")
    wl_obs2, _ = make_whitelist_obs(annotated.obs, original_columns)
    coerced_obs, touched = narrow_coerce(wl_obs2)
    print(f"  coerced {len(touched)} columns: {touched}")
    c_merged = merge_back(sub, coerced_obs)
    audit_obs(c_merged.obs, "Tier 2 obs (coerced)")
    ok2, msg2 = try_write_and_read(c_merged, os.path.join(scratch, "tier2_coerced.h5ad"))
    print(f"\n  >>> TIER 2 RESULT: {msg2}")

    # ===================================================================
    banner("SUMMARY")
    print(f"  Tier 0 RAW         : {'PASS' if ok0 else 'FAIL'}  ({msg0})")
    print(f"  Tier 1 WHITELIST   : {'PASS' if ok1 else 'FAIL'}  ({msg1})")
    print(f"  Tier 2 WHITELIST+  : {'PASS' if ok2 else 'FAIL'}  ({msg2})")
    print()
    if ok1:
        print("  -> Minimal whitelist is sufficient. Ship Jake's chunk in Step04.")
    elif ok2:
        print("  -> Whitelist alone insufficient; the narrow coercion closes it.")
        print("     Check the Tier 1 flagged columns above to see which shape did it.")
    else:
        print("  -> Both failed. Re-check the Tier 2 audit for an unhandled dtype.")
    print(f"\n  scratch outputs: {scratch}")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        sys.exit(1)
