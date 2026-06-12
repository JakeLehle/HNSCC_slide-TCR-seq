#!/usr/bin/env python3
"""
diag_h5py_fix_verify.py

Quick verification that h5py downgrade fixes the write issue.
Creates a tiny adata with the exact problematic column types from
the popV audit, tries to write. Takes ~5 seconds.

Run after: conda install h5py=3.15.1
"""

import numpy as np
import pandas as pd
import anndata as ad
import h5py
import tempfile, os

print(f"h5py version: {h5py.__version__}")
print(f"anndata version: {ad.__version__}")
print(f"pandas version: {pd.__version__}")

# Create minimal adata (10 cells, 5 genes)
adata = ad.AnnData(
    X=np.random.rand(10, 5).astype(np.float32),
    obs=pd.DataFrame(index=[f"cell_{i}" for i in range(10)]),
    var=pd.DataFrame(index=[f"gene_{i}" for i in range(5)]),
)

# Reproduce the exact problematic columns from the audit
# 1. observation_joinid: object dtype, ALL NaN (the confirmed crash culprit)
adata.obs["observation_joinid"] = pd.Series([np.nan] * 10, dtype=object)

# 2. Integer categoricals (popv_prediction_score, popv_majority_vote_score)
adata.obs["popv_prediction_score"] = pd.Categorical([4, 5, 3, 4, 5, 6, 7, 3, 4, 5])
adata.obs["popv_majority_vote_score"] = pd.Categorical([1, 2, 3, 4, 5, 1, 2, 3, 4, 5])

# 3. Empty categoricals (reference metadata columns, all NaN for query cells)
adata.obs["donor_id"] = pd.Categorical([np.nan] * 10)

# 4. numpy str_ categoricals (popv_onclass_seen)
adata.obs["popv_onclass_seen"] = pd.Series(
    [np.str_("basal cell")] * 10, dtype=object
)

# 5. Normal popV prediction column (should always work)
adata.obs["popv_prediction"] = pd.Categorical(
    ["basal cell", "T cell", "basal cell", "macrophage", "T cell",
     "basal cell", "T cell", "basal cell", "macrophage", "T cell"]
)

print(f"\nColumn types:")
for col in adata.obs.columns:
    dtype = adata.obs[col].dtype
    if hasattr(dtype, "categories"):
        print(f"  {col:40s} categorical, cat.dtype={dtype.categories.dtype}")
    else:
        print(f"  {col:40s} {dtype}")

# Try to write
path = os.path.join(tempfile.gettempdir(), "h5py_fix_test.h5ad")
try:
    adata.write_h5ad(path)
    size = os.path.getsize(path) / 1e3
    print(f"\nPASS: h5ad write succeeded ({size:.1f} KB)")
    print(f"h5py {h5py.__version__} handles all problematic dtypes correctly.")
    print(f"Safe to re-run Step04.")
    os.remove(path)
except TypeError as e:
    print(f"\nFAIL: {e}")
    print(f"h5py {h5py.__version__} still cannot write these dtypes.")
    print(f"Do NOT re-run Step04 yet.")
