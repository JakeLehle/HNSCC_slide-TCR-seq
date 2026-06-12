#!/usr/bin/env python3
"""
diag_popv_column_audit.py

Targeted diagnostic to find the ACTUAL problematic column(s) from popV.

The first diagnostic showed predictions.csv columns are fine. But popV
adds extra columns to adata.obs in memory (_predict_cells, popv_onclass_seen,
etc.) that never make it to predictions.csv. One of those is the likely
culprit.

Strategy: Run popV on a tiny subset (200 cells) to capture the full
in-memory column inventory, then test each column for h5ad write
compatibility. This takes ~5-10 minutes, not 5.5 hours.

Run in sc_pre env on login node (interactive).
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import popv
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# PATHS
# =============================================================================
QC_H5AD = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/03_qc/all_pucks_merged_QC.h5ad"
TEST_DIR = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/TROUBLESHOOTING/diag_popv_audit"
CACHE_DIR = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/04_annotation/popv_cache"
os.makedirs(TEST_DIR, exist_ok=True)

print("=" * 70)
print("POPV COLUMN AUDIT: Running popV on tiny subset")
print("=" * 70)

# =============================================================================
# STEP 1: Load a tiny subset (200 cells, one puck)
# =============================================================================
print("\n[STEP 1] Loading QC h5ad and subsetting to 200 cells...")
adata_full = sc.read_h5ad(QC_H5AD)
print(f"  Full: {adata_full.shape}")

# Take 200 cells from the first puck
mask = adata_full.obs['puck_id'] == adata_full.obs['puck_id'].cat.categories[0]
adata_sub = adata_full[mask].copy()[:200]
print(f"  Subset: {adata_sub.shape}")

# Record original columns
original_obs_columns = list(adata_sub.obs.columns)
print(f"  Original obs columns ({len(original_obs_columns)}): {original_obs_columns}")

# Ensure gene symbol column
if 'feature_name' not in adata_sub.var.columns:
    adata_sub.var['feature_name'] = adata_sub.var_names

# =============================================================================
# STEP 2: Run popV annotation on tiny subset
# =============================================================================
print(f"\n[STEP 2] Running popV annotation on 200 cells (should take ~5 min)...")

huggingface_repo = "popV/tabula_sapiens_All_Cells"
hmo = popv.hub.HubModel.pull_from_huggingface_hub(huggingface_repo, cache_dir=CACHE_DIR)

adata_annotated = hmo.annotate_data(
    adata_sub,
    query_batch_key="puck_id",
    prediction_mode="inference",
    gene_symbols="feature_name",
)

print(f"\n  Annotated adata shape: {adata_annotated.shape}")

# =============================================================================
# STEP 3: Inspect ALL columns popV added
# =============================================================================
print(f"\n[STEP 3] Full column inventory of annotated adata")
print(f"  Total obs columns: {len(adata_annotated.obs.columns)}")

new_columns = [c for c in adata_annotated.obs.columns if c not in original_obs_columns]
print(f"  New columns added by popV ({len(new_columns)}):")

for col in new_columns:
    series = adata_annotated.obs[col]
    dtype = series.dtype

    # Detailed type inspection
    info = f"    {col:50s} dtype={str(dtype):25s}"

    if hasattr(dtype, 'categories'):
        cat_dtype = series.cat.categories.dtype
        cat_values = list(series.cat.categories[:5])
        cat_py_types = set(type(v).__name__ for v in series.cat.categories)
        info += f" cat.dtype={str(cat_dtype):15s} cat_py_types={cat_py_types}"
        info += f" sample_cats={cat_values}"
    else:
        sample_values = series.dropna().head(3).tolist()
        py_types = set(type(v).__name__ for v in series.dropna().head(20))
        info += f" py_types={py_types} sample={sample_values}"

    print(info)

# =============================================================================
# STEP 4: Filter to query cells (mimicking Step04)
# =============================================================================
print(f"\n[STEP 4] Filter to query cells...")
if '_dataset' in adata_annotated.obs.columns:
    n_before = adata_annotated.n_obs
    adata_annotated = adata_annotated[
        adata_annotated.obs['_dataset'] == 'query'
    ].copy()
    print(f"  Filtered: {n_before} -> {adata_annotated.n_obs} cells")

# =============================================================================
# STEP 5: Simulate the merge (exactly as Step04 does it)
# =============================================================================
print(f"\n[STEP 5] Simulating Step04 merge...")

popv_obs = adata_annotated.obs.drop(columns=original_obs_columns, errors='ignore')
print(f"  popv_obs columns ({len(popv_obs.columns)}): {list(popv_obs.columns)}")

merged_df = pd.merge(
    adata_sub.obs, popv_obs,
    left_index=True, right_index=True, how='inner'
)
adata_merged = adata_sub[adata_sub.obs_names.isin(merged_df.index)].copy()
merged_df = merged_df.loc[adata_merged.obs_names]
adata_merged.obs = merged_df

print(f"  Merged adata: {adata_merged.shape}")
print(f"  Merged obs columns ({len(adata_merged.obs.columns)}):")

# =============================================================================
# STEP 6: Column-by-column write test on MERGED adata
# =============================================================================
print(f"\n[STEP 6] Column-by-column h5ad write test (merged adata)")

# Start with clean adata (no popV columns)
adata_test = adata_sub[adata_sub.obs_names.isin(adata_merged.obs_names)].copy()
test_path = os.path.join(TEST_DIR, "test_col.h5ad")

# Test baseline
try:
    adata_test.write_h5ad(test_path)
    print(f"  BASELINE (no popV columns): PASS")
    os.remove(test_path)
except Exception as e:
    print(f"  BASELINE: FAIL ({type(e).__name__}: {e})")

# Add popV columns one at a time
for col in popv_obs.columns:
    adata_test.obs[col] = merged_df.loc[adata_test.obs_names, col]
    try:
        adata_test.write_h5ad(test_path)
        print(f"    + {col:50s} PASS")
        os.remove(test_path)
    except Exception as e:
        dtype = adata_test.obs[col].dtype
        print(f"    + {col:50s} FAIL ({type(e).__name__})")
        print(f"      dtype: {dtype}")
        if hasattr(dtype, 'categories'):
            cats = adata_test.obs[col].cat.categories
            print(f"      categories.dtype: {cats.dtype}")
            print(f"      categories: {list(cats)}")
            print(f"      category python types: {set(type(v).__name__ for v in cats)}")
        else:
            vals = adata_test.obs[col].dropna().head(5)
            print(f"      sample values: {list(vals)}")
            print(f"      python types: {set(type(v).__name__ for v in vals)}")

        # Test if converting to string fixes it
        adata_test.obs[col] = adata_test.obs[col].astype(str).astype('category')
        try:
            adata_test.write_h5ad(test_path)
            print(f"      -> .astype(str).astype('category') FIXES IT")
            os.remove(test_path)
        except Exception as e2:
            print(f"      -> .astype(str).astype('category') still fails: {e2}")
            adata_test.obs.drop(columns=[col], inplace=True)

# =============================================================================
# STEP 7: Test full merged write
# =============================================================================
print(f"\n[STEP 7] Full merged adata write test...")
try:
    adata_merged.write_h5ad(os.path.join(TEST_DIR, "test_full_merged.h5ad"))
    print(f"  PASS: Full merged adata writes OK")
except Exception as e:
    print(f"  FAIL: {type(e).__name__}: {e}")
    print(f"\n  Applying targeted fixes and retrying...")

    # Fix: convert any non-string categoricals
    for col in adata_merged.obs.columns:
        if hasattr(adata_merged.obs[col].dtype, 'categories'):
            cat_types = set(type(v).__name__ for v in adata_merged.obs[col].cat.categories)
            if cat_types != {'str'}:
                print(f"    Fixing {col}: categories are {cat_types}, converting to str")
                adata_merged.obs[col] = adata_merged.obs[col].astype(str).astype('category')

    try:
        adata_merged.write_h5ad(os.path.join(TEST_DIR, "test_full_fixed.h5ad"))
        print(f"  PASS after fixing non-string categoricals")
    except Exception as e2:
        print(f"  STILL FAILS: {type(e2).__name__}: {e2}")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("COLUMN AUDIT COMPLETE")
print("=" * 70)
print(f"\n  Columns NOT in predictions.csv (popV-internal):")
predictions_cols = [
    'popv_celltypist_prediction', 'popv_knn_bbknn_prediction',
    'popv_knn_harmony_prediction', 'popv_knn_on_scvi_prediction',
    'popv_onclass_prediction', 'popv_scanvi_prediction',
    'popv_svm_prediction', 'popv_xgboost_prediction',
    'popv_prediction', 'popv_prediction_score',
    'popv_majority_vote_prediction', 'popv_majority_vote_score',
    'popv_parent'
]
extra_cols = [c for c in popv_obs.columns if c not in predictions_cols]
for col in extra_cols:
    dtype = popv_obs[col].dtype
    if hasattr(dtype, 'categories'):
        cat_types = set(type(v).__name__ for v in popv_obs[col].cat.categories)
        print(f"    {col:50s} dtype={dtype} cat_types={cat_types}")
    else:
        py_types = set(type(v).__name__ for v in popv_obs[col].dropna().head(20))
        print(f"    {col:50s} dtype={dtype} py_types={py_types}")
