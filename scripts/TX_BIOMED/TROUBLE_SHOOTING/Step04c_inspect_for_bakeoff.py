#!/usr/bin/env python3
"""
Step04c_inspect_for_bakeoff.py  (READ-ONLY diagnostic)

Inspects the annotated spatial object and the cell2location single-cell
reference so we know exactly what is on disk before writing the marker
bake-off. Prints a structured report; writes nothing.

Env: sc_pre, on theia (TX Biomed).
NOTE: paths are hardcoded to the TX Biomed tree on purpose. spatial_config.py
still carries the UTSA PROJECT_ROOT and its startup guard would abort here, so
this script deliberately does NOT import it.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as ssp

PROJECT_ROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
ANNOTATED = os.path.join(PROJECT_ROOT, "data/outputs/04_annotation/all_pucks_annotated.h5ad")
REFERENCE = os.path.join("/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/reference_for_cell2location.h5ad")


def mat_stats(M):
    """min/max/dtype + integer/negative check on a sample of the first rows."""
    try:
        sub = M[:1000]
        data = sub.data if ssp.issparse(sub) else np.asarray(sub).ravel()
        if data.size == 0:
            return f"dtype={getattr(M, 'dtype', '?')}, empty sample"
        mn, mx = float(np.nanmin(data)), float(np.nanmax(data))
        return (f"dtype={getattr(M, 'dtype', '?')}, min={mn:.4g}, max={mx:.4g}, "
                f"integer_like={np.allclose(data, np.round(data))}, "
                f"has_negative={mn < 0}  (first 1000-row sample)")
    except Exception as e:
        return f"could not sample ({e})"


def hr(t):
    print("\n" + "=" * 78 + f"\n  {t}\n" + "=" * 78)


def inspect(path, label, is_reference):
    hr(label)
    if not os.path.exists(path):
        print(f"  MISSING: {path}")
        return
    print(f"  path: {path}")
    print(f"  size: {os.path.getsize(path) / 1e6:.1f} MB")
    ad = sc.read_h5ad(path, backed='r')
    print(f"  shape: {ad.n_obs:,} obs x {ad.n_vars:,} vars")
    print(f"  X: {mat_stats(ad.X)}")
    print(f"  layers: {list(ad.layers.keys())}")
    for k in ad.layers:
        print(f"    layer['{k}']: {mat_stats(ad.layers[k])}")
    print(f"  var_names[:10]: {list(ad.var_names[:10])}")
    print(f"  ENSG-looking var_names in first 2000: "
          f"{sum(str(v).startswith('ENSG') for v in ad.var_names[:2000])}")
    print(f"  obsm keys: {list(ad.obsm.keys())}")
    print(f"  uns keys:  {list(ad.uns.keys())}")
    print(f"  obs columns ({ad.obs.shape[1]}): {list(ad.obs.columns)}")

    if is_reference:
        for col in ['final_annotation', 'subject id', 'subject_id', 'cell_type']:
            if col in ad.obs.columns:
                print(f"\n  obs['{col}'] value_counts:")
                for k, v in ad.obs[col].astype(str).value_counts().items():
                    print(f"    {k}: {v:,}")
    else:
        for col in ['puck_id', 'popv_prediction', 'c2l_argmax',
                    'final_annotation', 'consensus_annotation', 'consensus_agreement']:
            if col in ad.obs.columns:
                vc = ad.obs[col].astype(str).value_counts()
                print(f"\n  obs['{col}'] ({vc.shape[0]} levels, top 15):")
                for k, v in vc.head(15).items():
                    print(f"    {k}: {v:,}")
        if 'clusters' in ad.obs.columns:
            cl = ad.obs['clusters'].astype(str).value_counts()
            print(f"\n  obs['clusters']: {cl.shape[0]} clusters, sizes {cl.min():,}-{cl.max():,}")
        if 'popv_prediction_score' in ad.obs.columns:
            s = pd.to_numeric(ad.obs['popv_prediction_score'], errors='coerce')
            print(f"  popv_prediction_score: min={s.min():.3g}, median={s.median():.3g}, "
                  f"max={s.max():.3g}, n_nan={int(s.isna().sum())}")
        c2l_cols = [c for c in ad.obs.columns if c.startswith('c2l_') and c != 'c2l_argmax']
        print(f"  c2l abundance obs columns ({len(c2l_cols)}): {c2l_cols}")
        for col in ['x_coord', 'y_coord']:
            print(f"  obs has '{col}': {col in ad.obs.columns}")


if __name__ == "__main__":
    inspect(ANNOTATED, "ANNOTATED SPATIAL OBJECT (all_pucks_annotated.h5ad)", is_reference=False)
    inspect(REFERENCE, "CELL2LOCATION SINGLE-CELL REFERENCE", is_reference=True)
    print("\nDONE. Paste this whole report back.")
