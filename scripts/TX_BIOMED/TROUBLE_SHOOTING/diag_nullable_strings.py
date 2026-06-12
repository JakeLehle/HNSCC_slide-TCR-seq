#!/usr/bin/env python3
"""
diag_nullable_strings.py

Diagnostic to reproduce and isolate the h5ad write failure from Step04.
Loads the QC h5ad + popV predictions.csv, merges them, then tests writing
column-by-column to identify exactly which column(s) carry the problematic
dtype that h5py can't serialize.

Run in sc_pre env on login node (no SLURM needed, just memory).
"""

import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# PATHS
# =============================================================================
QC_H5AD = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/03_qc/all_pucks_merged_QC.h5ad"
PREDICTIONS_CSV = "/master/jlehle/WORKING/slide-TCR-seq-working/scripts/tmp/popv_output/predictions.csv"
TEST_DIR = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/TROUBLESHOOTING/diag_nullable"
os.makedirs(TEST_DIR, exist_ok=True)

# =============================================================================
# STEP 1: Report package versions
# =============================================================================
print("=" * 70)
print("STEP 1: Package versions")
print("=" * 70)
print(f"  Python:   {sys.version}")
print(f"  pandas:   {pd.__version__}")
print(f"  numpy:    {np.__version__}")
print(f"  anndata:  {ad.__version__}")
print(f"  scanpy:   {sc.__version__}")

try:
    import popv
    print(f"  popv:     {popv.__version__}")
except Exception as e:
    print(f"  popv:     import failed ({e})")

try:
    import squidpy
    print(f"  squidpy:  {squidpy.__version__}")
except Exception:
    print(f"  squidpy:  not installed")

try:
    import h5py
    print(f"  h5py:     {h5py.__version__}")
except Exception:
    print(f"  h5py:     not installed")

# Check pandas string settings
print(f"\n  pd.options.future.infer_string = {pd.options.future.infer_string}")
print(f"  pd.api.types.infer_dtype test:")
test_series = pd.Series(["a", "b", "c"])
print(f"    pd.Series(['a','b','c']).dtype = {test_series.dtype}")
test_cat = pd.Categorical(["a", "b", "c"])
print(f"    pd.Categorical(['a','b','c']).categories.dtype = {test_cat.categories.dtype}")

# =============================================================================
# STEP 2: Load and inspect predictions.csv
# =============================================================================
print("\n" + "=" * 70)
print("STEP 2: Inspect predictions.csv")
print("=" * 70)

preds = pd.read_csv(PREDICTIONS_CSV, index_col=0)
print(f"  Shape: {preds.shape}")
print(f"  Columns: {list(preds.columns)}")
print(f"\n  Column dtypes from CSV load:")
for col in preds.columns:
    dtype = preds[col].dtype
    nunique = preds[col].nunique()
    n_null = preds[col].isna().sum()
    sample = preds[col].dropna().iloc[0] if len(preds[col].dropna()) > 0 else "ALL_NULL"
    print(f"    {col:45s} dtype={str(dtype):20s} unique={nunique:5d} null={n_null:5d} sample='{sample}'")

# =============================================================================
# STEP 3: Load QC h5ad and check baseline
# =============================================================================
print("\n" + "=" * 70)
print("STEP 3: Load QC h5ad (baseline)")
print("=" * 70)

adata = sc.read_h5ad(QC_H5AD)
print(f"  Shape: {adata.shape}")
print(f"  obs columns: {list(adata.obs.columns)}")
print(f"\n  Baseline obs dtypes:")
for col in adata.obs.columns:
    dtype = adata.obs[col].dtype
    if hasattr(dtype, 'categories'):
        cat_dtype = adata.obs[col].cat.categories.dtype
        print(f"    {col:45s} dtype={str(dtype):20s} categories.dtype={cat_dtype}")
    else:
        print(f"    {col:45s} dtype={str(dtype):20s}")

# Test baseline write
print(f"\n  Testing baseline h5ad write...")
test_path = os.path.join(TEST_DIR, "test_baseline.h5ad")
try:
    adata.write_h5ad(test_path)
    print(f"  PASS: Baseline writes OK ({os.path.getsize(test_path) / 1e6:.1f} MB)")
    os.remove(test_path)
except Exception as e:
    print(f"  FAIL: Baseline write error: {type(e).__name__}: {e}")

# =============================================================================
# STEP 4: Simulate the merge (what popV does internally)
# =============================================================================
print("\n" + "=" * 70)
print("STEP 4: Merge predictions into adata (simulating popV)")
print("=" * 70)

# Align indices
common = adata.obs_names.intersection(preds.index)
print(f"  Common barcodes: {len(common)} / {adata.n_obs} adata, {len(preds)} preds")

# Subset adata to common
adata_merged = adata[common].copy()

# Add prediction columns
for col in preds.columns:
    vals = preds.loc[adata_merged.obs_names, col]
    # popV stores these as categoricals
    adata_merged.obs[col] = pd.Categorical(vals)

print(f"  Merged shape: {adata_merged.shape}")
print(f"  Merged obs columns: {list(adata_merged.obs.columns)}")

# =============================================================================
# STEP 5: Inspect merged dtypes (the culprit hunt)
# =============================================================================
print("\n" + "=" * 70)
print("STEP 5: Inspect merged obs dtypes")
print("=" * 70)

for col in adata_merged.obs.columns:
    dtype = adata_merged.obs[col].dtype
    if hasattr(dtype, 'categories'):
        cat_dtype = adata_merged.obs[col].cat.categories.dtype
        is_nullable = "string" in str(cat_dtype).lower() or "arrow" in str(cat_dtype).lower()
        flag = " *** SUSPECT ***" if is_nullable else ""
        print(f"    {col:45s} dtype={str(dtype):20s} categories.dtype={str(cat_dtype):20s}{flag}")
    else:
        is_nullable = "string" in str(dtype).lower() or "arrow" in str(dtype).lower()
        flag = " *** SUSPECT ***" if is_nullable else ""
        print(f"    {col:45s} dtype={str(dtype):20s}{flag}")

# Also check for non-string objects hiding in string columns
print(f"\n  Checking for mixed-type columns:")
for col in adata_merged.obs.columns:
    series = adata_merged.obs[col]
    if series.dtype == "object":
        types_found = set(type(x).__name__ for x in series.dropna().head(100))
        if len(types_found) > 1:
            print(f"    {col:45s} MIXED TYPES: {types_found}")

# =============================================================================
# STEP 6: Column-by-column write test
# =============================================================================
print("\n" + "=" * 70)
print("STEP 6: Column-by-column h5ad write test")
print("=" * 70)

# Start with baseline adata (no predictions), add one column at a time
adata_test = adata[common].copy()
test_path = os.path.join(TEST_DIR, "test_colwise.h5ad")

for col in preds.columns:
    vals = preds.loc[adata_test.obs_names, col]
    adata_test.obs[col] = pd.Categorical(vals)
    try:
        adata_test.write_h5ad(test_path)
        print(f"    + {col:45s} PASS")
        os.remove(test_path)
    except Exception as e:
        print(f"    + {col:45s} FAIL ({type(e).__name__}: {e})")
        # Diagnostic: what is in this column?
        print(f"      dtype:            {adata_test.obs[col].dtype}")
        if hasattr(adata_test.obs[col].dtype, 'categories'):
            print(f"      categories.dtype: {adata_test.obs[col].cat.categories.dtype}")
            print(f"      categories[:5]:   {list(adata_test.obs[col].cat.categories[:5])}")
        print(f"      sample values:    {list(adata_test.obs[col].dropna().head(5))}")
        print(f"      null count:       {adata_test.obs[col].isna().sum()}")
        # Remove the failing column so we can continue testing the rest
        adata_test.obs.drop(columns=[col], inplace=True)

# =============================================================================
# STEP 7: Test the revert_nullable_strings fix
# =============================================================================
print("\n" + "=" * 70)
print("STEP 7: Test revert_nullable_strings fix")
print("=" * 70)

def revert_nullable_strings(adata):
    """
    Revert nullable/arrow string dtypes to object dtype for h5ad compatibility.
    From: https://github.com/scverse/anndata/issues/2221
    """
    pd.options.future.infer_string = False

    adata.obs.index = adata.obs.index.astype(object)
    adata.var.index = adata.var.index.astype(object)
    if adata.raw is not None:
        adata.raw.var.index = adata.raw.var.index.astype(object)

    obs_cols = adata.obs.select_dtypes(include=['category', 'string[python]']).columns
    for col in obs_cols:
        adata.obs[col] = pd.Series(adata.obs[col], dtype="object")
        adata.obs[col] = adata.obs[col].astype("category")

    var_cols = adata.var.select_dtypes(include=['category', 'string[python]']).columns
    for col in var_cols:
        adata.var[col] = pd.Series(adata.var[col], dtype="object")
        adata.var[col] = adata.var[col].astype("category")

    return adata


# Re-create the full merge
adata_fix = adata[common].copy()
for col in preds.columns:
    vals = preds.loc[adata_fix.obs_names, col]
    adata_fix.obs[col] = pd.Categorical(vals)

print(f"  Before fix:")
for col in preds.columns:
    if col in adata_fix.obs.columns and hasattr(adata_fix.obs[col].dtype, 'categories'):
        print(f"    {col:45s} categories.dtype = {adata_fix.obs[col].cat.categories.dtype}")

adata_fix = revert_nullable_strings(adata_fix)

print(f"\n  After fix:")
for col in preds.columns:
    if col in adata_fix.obs.columns and hasattr(adata_fix.obs[col].dtype, 'categories'):
        print(f"    {col:45s} categories.dtype = {adata_fix.obs[col].cat.categories.dtype}")

fix_path = os.path.join(TEST_DIR, "test_fixed.h5ad")
try:
    adata_fix.write_h5ad(fix_path)
    print(f"\n  PASS: Fixed adata writes OK ({os.path.getsize(fix_path) / 1e6:.1f} MB)")
    os.remove(fix_path)
except Exception as e:
    print(f"\n  FAIL: Fix did not resolve: {type(e).__name__}: {e}")

# =============================================================================
# SUMMARY
# =============================================================================
print("\n" + "=" * 70)
print("DIAGNOSTIC COMPLETE")
print("=" * 70)
print(f"  predictions.csv path: {PREDICTIONS_CSV}")
print(f"  predictions.csv rows: {len(preds)}")
print(f"  QC h5ad cells:        {adata.n_obs}")
print(f"  Common after merge:   {len(common)}")
print(f"  Troubleshooting dir:  {TEST_DIR}")
