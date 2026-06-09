#!/usr/bin/env python3
"""
build_cell2location_reference.py

Combine the raw-count object (adata_tmp.h5ad: QC'd, bad cells dropped, NOT
normalized) with the labeled object (adata_final.h5ad: has final_annotation
but no raw counts) into a single cell2location reference that carries raw
counts AND cell-type labels.

The join is strictly by barcode (obs_names) and gene name (var_names), never
positional, and the script refuses to build unless every labeled cell is
found in the count source.

Two build modes:
  tmp_reference (default, recommended): reference = adata_tmp subset to the
      labeled cells, with final_annotation copied over by barcode. Keeps
      tmp's native (cleaner) gene space; raw counts stay in X.
  layer_on_final: your original idea -- align tmp counts to adata_final's
      obs/var by name and attach as adata_final.layers['counts']. Inherits
      final's mixed Ensembl/symbol var space; genes in final but absent from
      tmp become zero columns (reported).

Default is VERIFY ONLY. Pass --build to write a new file (originals untouched).

Run (sc_pre env):
  python build_cell2location_reference.py                       # verify only
  python build_cell2location_reference.py --build               # tmp_reference
  python build_cell2location_reference.py --build --mode layer_on_final

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

DIR = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input"
COUNTS_DEFAULT = os.path.join(DIR, "adata_tmp.h5ad")        # raw counts, QC'd
RAWEST_DEFAULT = os.path.join(DIR, "adata.h5ad")            # fallback (pre-QC)
LABELED_DEFAULT = os.path.join(DIR, "adata_final.h5ad")     # labels, normalized
SPATIAL_DEFAULT = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/03_qc/all_pucks_merged_QC.h5ad"
OUT_DIR = DIR
ANN_COL = "final_annotation"


def banner(m):
    print("\n" + "=" * 72); print(m); print("=" * 72)


def is_int_counts(X, n=2000):
    vals = X.data[:n] if sparse.issparse(X) else np.asarray(X).ravel()
    vals = vals[vals != 0][:n]
    return (len(vals) == 0) or (float(np.mean(np.isclose(vals, np.round(vals)))) > 0.999)


def namespace_report(var_names, label, spatial_genes=None):
    vn = var_names.astype(str)
    n_ens = int(vn.str.startswith('ENSG').sum())
    msg = f"  {label}: {len(vn):,} genes, {n_ens} Ensembl, {int(vn.duplicated().sum())} dups"
    ov = None
    if spatial_genes is not None:
        ov = len(set(vn) & spatial_genes)
        msg += f", spatial symbol overlap={ov:,}"
    print(msg)
    return ov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=COUNTS_DEFAULT)
    ap.add_argument("--labeled", default=LABELED_DEFAULT)
    ap.add_argument("--rawest", default=RAWEST_DEFAULT)
    ap.add_argument("--spatial", default=SPATIAL_DEFAULT)
    ap.add_argument("--col", default=ANN_COL)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--mode", choices=["tmp_reference", "layer_on_final"], default="tmp_reference")
    args = ap.parse_args()

    banner("BUILD cell2location REFERENCE (verify%s)" % ("" if not args.build else " + build"))

    # spatial gene set
    spatial_genes = None
    if os.path.exists(args.spatial):
        sp = sc.read_h5ad(args.spatial, backed='r')
        spatial_genes = set(sp.var_names.astype(str))
        print(f"  spatial genes: {len(spatial_genes):,}")

    # load
    for p, name in [(args.counts, "counts source"), (args.labeled, "labeled source")]:
        if not os.path.exists(p):
            print(f"  FATAL: {name} not found: {p}"); sys.exit(1)
    counts = sc.read_h5ad(args.counts)
    labeled = sc.read_h5ad(args.labeled)
    print(f"  counts  (adata_tmp):   {counts.n_obs:,} cells x {counts.n_vars:,} genes")
    print(f"  labeled (adata_final): {labeled.n_obs:,} cells x {labeled.n_vars:,} genes")

    # ---- VERIFY ----------------------------------------------------------
    banner("VERIFY")
    # 1. raw counts present in counts source
    int_ok = is_int_counts(counts.X)
    print(f"  counts source X is raw integer counts: {int_ok}")
    if not int_ok:
        print("  WARNING: counts source X is not integer; checking layers/raw...")
        for lyr in counts.layers:
            if is_int_counts(counts.layers[lyr]):
                print(f"    raw counts found in layers['{lyr}'] (use --counts layer manually)")

    # 2. annotation present
    if args.col not in labeled.obs.columns:
        print(f"  FATAL: '{args.col}' not in labeled.obs"); sys.exit(1)

    # 3. cell alignment (every labeled cell must exist in counts source)
    lab_cells = set(labeled.obs_names.astype(str))
    cnt_cells = set(counts.obs_names.astype(str))
    missing = lab_cells - cnt_cells
    print(f"  labeled cells: {len(lab_cells):,}   present in counts source: "
          f"{len(lab_cells & cnt_cells):,}   MISSING: {len(missing):,}")
    if missing:
        print(f"    e.g. {list(missing)[:3]}")
        if os.path.exists(args.rawest):
            rawest = sc.read_h5ad(args.rawest, backed='r')
            still = missing - set(rawest.obs_names.astype(str))
            print(f"    of those, still missing from adata.h5ad (rawest): {len(still):,}")
        print("    -> not all labeled cells have raw counts; resolve before building")

    # 4. gene namespaces + overlaps
    print()
    ov_tmp = namespace_report(counts.var_names, "counts source (tmp)", spatial_genes)
    ov_fin = namespace_report(labeled.var_names, "labeled (final)", spatial_genes)
    fin_in_tmp = len(set(labeled.var_names.astype(str)) & set(counts.var_names.astype(str)))
    print(f"  final genes also in tmp (for layer mode): {fin_in_tmp:,} / {labeled.n_vars:,}")

    # recommendation
    banner("RECOMMENDATION")
    can_build = (len(missing) == 0) and int_ok
    print(f"  buildable: {can_build}")
    if spatial_genes is not None and ov_tmp is not None and ov_fin is not None:
        better = "tmp_reference" if ov_tmp >= ov_fin else "layer_on_final"
        print(f"  spatial overlap: tmp={ov_tmp:,}  final={ov_fin:,}  -> prefer {better}")
    if not args.build:
        print("\n  verify only. Re-run with --build (and optional --mode) to write the reference.")
        return
    if not can_build:
        print("\n  refusing to build: fix the issues above first (raw counts / missing cells).")
        sys.exit(1)

    # ---- BUILD -----------------------------------------------------------
    banner(f"BUILD ({args.mode})")
    lab_map = labeled.obs[args.col].copy()
    lab_map.index = labeled.obs_names.astype(str)

    if args.mode == "tmp_reference":
        keep = counts.obs_names.astype(str).isin(lab_cells)
        ref = counts[keep].copy()
        ref.obs[args.col] = lab_map.reindex(ref.obs_names.astype(str)).values
        n_nan = pd.isna(ref.obs[args.col]).sum()
        assert n_nan == 0, f"{n_nan} cells failed label assignment"
        # store counts explicitly in a layer too, for cell2location layer= arg
        ref.layers['counts'] = ref.X.copy()
        print(f"  reference: {ref.n_obs:,} cells x {ref.n_vars:,} genes")
        print(f"  label counts:\n{ref.obs[args.col].value_counts().to_string()}")
        assert is_int_counts(ref.X), "reference X is not integer counts!"
        out = os.path.join(args.out, "reference_for_cell2location.h5ad")
        ref.write_h5ad(out)

    else:  # layer_on_final
        # align tmp counts to final's (obs, var) strictly by name
        tmp_sub = counts[counts.obs_names.astype(str).isin(lab_cells)].copy()
        tmp_sub = tmp_sub[labeled.obs_names.astype(str)]          # final cell order
        gene_pos = {g: i for i, g in enumerate(tmp_sub.var_names.astype(str))}
        cols = np.array([gene_pos.get(g, -1) for g in labeled.var_names.astype(str)])
        present = cols >= 0
        Xc = tmp_sub.X.tocsc()
        new = sparse.lil_matrix((labeled.n_obs, labeled.n_vars), dtype=Xc.dtype)
        new[:, np.where(present)[0]] = Xc[:, cols[present]]
        out_adata = labeled.copy()
        out_adata.layers['counts'] = new.tocsr()
        print(f"  attached counts layer to final: {labeled.n_obs:,} x {labeled.n_vars:,}")
        print(f"  final genes with counts: {int(present.sum()):,}  "
              f"zero-filled (absent in tmp): {int((~present).sum()):,}")
        assert is_int_counts(out_adata.layers['counts']), "counts layer not integer!"
        out = os.path.join(args.out, "adata_final_with_counts.h5ad")
        out_adata.write_h5ad(out)

    print(f"\n  wrote: {out}")
    print("  (cell2location: point setup_anndata at layer='counts')")


if __name__ == "__main__":
    main()
