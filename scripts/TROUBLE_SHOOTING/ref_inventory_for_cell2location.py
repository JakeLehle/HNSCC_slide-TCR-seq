#!/usr/bin/env python3
"""
ref_inventory_for_cell2location.py

Small diagnostic: is the HNSCC single-cell object a usable cell2location
reference for the Slide-seq spatial data? Checks the things cell2location
actually depends on:

  1. Cell-type granularity + per-type counts (does it have the fine T / B /
     APC subtypes we need, and enough cells per type to learn a signature?)
  2. Raw counts (cell2location's RegressionModel needs raw counts; finds
     whether they live in X, a layer, or .raw)
  3. Gene space (symbols, and overlap with the spatial gene set)

Also peeks at a second object (adata_v_pp.h5ad) in case adata_final was
subset to basal cells for the NMF paper and the fuller reference is there.

Read-only. Saves a cell-type count table + bar plot.

Run (base or sc_pre env):
  python ref_inventory_for_cell2location.py
  python ref_inventory_for_cell2location.py --ref /abs/adata_final.h5ad

Author: Jake Lehle, Texas Biomedical Research Institute
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REF_DEFAULT = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/adata_final.h5ad"
ALT_DEFAULT = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/adata_v_pp.h5ad"
SPATIAL_DEFAULT = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/03_qc/all_pucks_merged_QC.h5ad"
OUT_DEFAULT = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/TROUBLESHOOTING/reference_check"

MIN_CELLS_PER_TYPE = 25  # cell2location reference regression rule-of-thumb floor

ANNOTATION_CANDIDATES = [
    'final_annotation', 'cell_type', 'celltype', 'cell_types',
    'annotation', 'CellType', 'leiden_annotation',
]

IMMUNE_KEYWORDS = {
    'T cell':      ['t cell', 't-cell', 'cd8', 'cd4', 'treg', 'regulatory t',
                    'cytotoxic', 'helper', 'tfh', 'mait', 'gamma', 'nkt', 'exhausted'],
    'B / plasma':  ['b cell', 'b-cell', 'plasma', 'plasmablast', 'germinal'],
    'Myeloid/APC': ['macrophage', 'monocyte', 'dendritic', ' dc', 'cdc', 'pdc',
                    'langerhans', 'myeloid', 'kupffer', 'tam'],
    'NK':          ['nk cell', 'natural killer', ' nk'],
}


def banner(m):
    print("\n" + "=" * 72); print(m); print("=" * 72)


def is_int_counts(X, n=2000):
    vals = X.data[:n] if sparse.issparse(X) else np.asarray(X).ravel()
    vals = vals[vals != 0][:n]
    return (len(vals) == 0) or (float(np.mean(np.isclose(vals, np.round(vals)))) > 0.999)


def find_raw_counts(adata):
    """Report where raw integer counts live (X / layer / raw)."""
    found = []
    if is_int_counts(adata.X):
        found.append("X")
    for lyr in adata.layers.keys():
        if is_int_counts(adata.layers[lyr]):
            found.append(f"layers['{lyr}']")
    if adata.raw is not None:
        try:
            if is_int_counts(adata.raw.X):
                found.append(f".raw (n_vars={adata.raw.n_vars})")
        except Exception:
            pass
    return found


def find_annotation_col(adata, prefer=None):
    if prefer and prefer in adata.obs.columns:
        return prefer
    for c in ANNOTATION_CANDIDATES:
        if c in adata.obs.columns:
            return c
    # fallback: a categorical obs col with a reasonable number of categories
    for c in adata.obs.columns:
        s = adata.obs[c]
        if (s.dtype == object or isinstance(s.dtype, pd.CategoricalDtype)):
            nun = s.nunique()
            if 3 <= nun <= 100:
                return c
    return None


def bucket_immune(categories):
    cats = [str(c) for c in categories]
    print("\n  Immune subtype coverage:")
    for group, kws in IMMUNE_KEYWORDS.items():
        hits = [c for c in cats if any(k in c.lower() for k in kws)]
        flag = "" if hits else "   <-- none found"
        print(f"    {group:<12}: {len(hits)} types{flag}")
        for h in hits:
            print(f"        - {h}")


def inventory(path, label, spatial_genes=None, prefer_col=None, full=True):
    if not os.path.exists(path):
        print(f"\n[{label}] MISSING: {path}")
        return None
    banner(f"[{label}] {path}")
    adata = sc.read_h5ad(path)
    print(f"  shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    ann = find_annotation_col(adata, prefer=prefer_col)
    if ann is None:
        print("  no annotation-like obs column found; obs columns:")
        print(f"    {list(adata.obs.columns)}")
        return adata
    print(f"  annotation column: '{ann}'  ({adata.obs[ann].nunique()} types)")

    vc = adata.obs[ann].value_counts()
    print(f"\n  Cell-type counts (flagging < {MIN_CELLS_PER_TYPE} cells):")
    for ct, n in vc.items():
        flag = "   <-- too few" if n < MIN_CELLS_PER_TYPE else ""
        print(f"    {str(ct):<45} {n:>7,}{flag}")
    n_thin = int((vc < MIN_CELLS_PER_TYPE).sum())

    bucket_immune(vc.index)

    if full:
        # raw counts
        print("\n  Raw integer counts located in:", find_raw_counts(adata) or "NONE FOUND")
        # gene space
        n_ens = int(adata.var_names.to_series().astype(str).str.startswith('ENSG').sum())
        print(f"  var_names: {n_ens}/{adata.n_vars} Ensembl "
              f"({'SYMBOLS ok' if n_ens == 0 else 'CONTAINS ENSEMBL'})  "
              f"dups={int(adata.var_names.duplicated().sum())}")
        if spatial_genes is not None:
            ov = len(set(adata.var_names.astype(str)) & spatial_genes)
            print(f"  gene overlap with spatial: {ov:,} / spatial {len(spatial_genes):,} "
                  f"({100*ov/len(spatial_genes):.1f}% of spatial genes covered)")

    return adata, ann, vc, n_thin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", default=REF_DEFAULT)
    ap.add_argument("--alt", default=ALT_DEFAULT)
    ap.add_argument("--spatial", default=SPATIAL_DEFAULT)
    ap.add_argument("--col", default="final_annotation")
    ap.add_argument("--out", default=OUT_DEFAULT)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    banner("REFERENCE INVENTORY FOR cell2location")

    # spatial gene set (var only; cheap via backed read)
    spatial_genes = None
    if os.path.exists(args.spatial):
        try:
            sp = sc.read_h5ad(args.spatial, backed='r')
            spatial_genes = set(sp.var_names.astype(str))
            print(f"  spatial gene set: {len(spatial_genes):,} genes "
                  f"from {os.path.basename(args.spatial)}")
        except Exception as e:
            print(f"  could not read spatial genes ({e}); overlap check skipped")
    else:
        print(f"  spatial object not found; overlap check skipped")

    # primary reference
    res = inventory(args.ref, "adata_final (candidate reference)",
                    spatial_genes=spatial_genes, prefer_col=args.col, full=True)

    # peek at the alternate object's annotation only
    if os.path.exists(args.alt):
        inventory(args.alt, "adata_v_pp (alternate, annotation peek)",
                  spatial_genes=spatial_genes, prefer_col=args.col, full=False)

    # save table + barplot for the primary reference
    if isinstance(res, tuple):
        _, ann, vc, n_thin = res
        vc.to_csv(os.path.join(args.out, "reference_celltype_counts.tsv"), sep='\t',
                  header=['n_cells'])
        fig, axx = plt.subplots(figsize=(10, max(4, 0.3 * len(vc))))
        vc.sort_values().plot(kind='barh', ax=axx, color='#3b6fb5')
        axx.axvline(MIN_CELLS_PER_TYPE, color='#c0392b', ls='--', lw=1,
                    label=f'{MIN_CELLS_PER_TYPE}-cell floor')
        axx.set_xlabel('cells'); axx.set_title(f"reference cell types ({ann})")
        axx.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(args.out, "reference_celltype_counts.png"),
                    dpi=150, bbox_inches='tight')
        plt.close()

        banner("VERDICT (reference suitability for cell2location)")
        print(f"  cell types: {len(vc)}   below {MIN_CELLS_PER_TYPE}-cell floor: {n_thin}")
        print("  Usable as a cell2location reference if: it has the fine immune")
        print("  subtypes you need (see immune coverage above), most types clear the")
        print("  cell floor, and raw counts were located. If the immune buckets are")
        print("  empty or thin, this object was likely basal-subset; use adata_v_pp")
        print("  (or the fuller upstream object) instead.")
    print(f"\n  outputs: {args.out}")


if __name__ == "__main__":
    main()
