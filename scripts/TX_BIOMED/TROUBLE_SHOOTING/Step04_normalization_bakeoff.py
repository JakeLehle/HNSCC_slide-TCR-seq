#!/usr/bin/env python3
"""
Step04_normalization_bakeoff.py

Compare normalization strategies for the sparse Slide-seq embedding, on the
data we already have. Reuses the popV checkpoint (raw counts + popv_prediction
labels) so NO popV rerun is needed.

Methods compared (all from raw counts; no normalization is double-applied):
  cpm_log1p   standard scRNA-seq: normalize_total + log1p + HVG + scale + PCA
  pearson     analytic Pearson residuals (scanpy's SCTransform-equivalent):
              pearson_residuals HVG -> normalize_pearson_residuals -> PCA
  cpm_regress (optional, --regress) cpm_log1p with total_counts regressed out

For each method it computes UMAP and two objective scores so the choice is
evidence-based rather than visual:
  depth_umap_corr   |Spearman(depth, UMAP axis)|, max over the 2 axes.
                    LOWER is better (embedding is less library-size driven).
  label_silhouette  silhouette of popv_prediction on PCA (subsampled).
                    HIGHER is better (cell types separate more).

It also writes a grid figure: rows = methods, cols = colored by
popv_prediction / depth / puck_id.

IMPORTANT context (from the spatial-normalization literature):
  No bead-level normalization is expected to break the depth ceiling for
  per-bead cell typing on ~350-gene Slide-seq. This bake-off is to confirm
  whether Pearson residuals buy meaningful separation before deciding
  between (a) using per-bead popv_prediction as-is, (b) re-mapping against a
  matched HNSCC single-cell reference, or (c) spatial binning.

Run (sc_pre env):
  python Step04_normalization_bakeoff.py
  python Step04_normalization_bakeoff.py --regress          # add depth-regression arm
  python Step04_normalization_bakeoff.py --subsample 40000  # faster on full data

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
N_TOP_GENES = 2000
N_PCS = 40


def banner(msg):
    print("\n" + "=" * 72)
    print(msg)
    print("=" * 72)


DEFAULT_OUTPUTS = "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs"


def _has_real_data(outputs_dir):
    """True only if this dir holds actual data files (not a stray empty tree)."""
    markers = [
        os.path.join(outputs_dir, "04_annotation", "all_pucks_popv_annotated.h5ad"),
        os.path.join(outputs_dir, "03_qc", "all_pucks_merged_QC.h5ad"),
    ]
    return any(os.path.isfile(m) for m in markers)


def resolve_outputs(arg):
    """Deterministic: explicit arg > known-good default > data-bearing walk-up."""
    if arg:
        return os.path.abspath(os.path.expanduser(arg))
    if _has_real_data(DEFAULT_OUTPUTS):
        return DEFAULT_OUTPUTS
    d = THIS_DIR
    for _ in range(8):
        cand = os.path.join(d, "data", "outputs")
        if _has_real_data(cand):
            return cand
        if os.path.basename(d) == "outputs" and _has_real_data(d):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    return DEFAULT_OUTPUTS


def is_integer_counts(X, n=2000):
    vals = X.data[:n] if sparse.issparse(X) else np.asarray(X).ravel()
    vals = vals[vals != 0][:n]
    if len(vals) == 0:
        return True
    return float(np.mean(np.isclose(vals, np.round(vals)))) > 0.999


def ensure_depth(adata):
    if 'total_counts' not in adata.obs:
        cnt = adata.layers['counts'] if 'counts' in adata.layers else adata.X
        adata.obs['total_counts'] = np.asarray(cnt.sum(axis=1)).ravel()
    adata.obs['log10_counts'] = np.log10(adata.obs['total_counts'].values + 1)
    return adata


# -----------------------------------------------------------------------------
# methods: each takes a raw-count copy, returns adata with X_pca + X_umap
# -----------------------------------------------------------------------------
def method_cpm_log1p(a):
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_TOP_GENES, flavor='seurat_v3', layer='counts')
    a = a[:, a.var['highly_variable']].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50, svd_solver='arpack')
    return a


def method_pearson(a):
    a.layers['counts'] = a.X.copy()
    sc.experimental.pp.highly_variable_genes(a, flavor='pearson_residuals', n_top_genes=N_TOP_GENES)
    a = a[:, a.var['highly_variable']].copy()
    sc.experimental.pp.normalize_pearson_residuals(a)
    sc.tl.pca(a, n_comps=50, svd_solver='arpack')
    return a


def method_cpm_regress(a):
    a.layers['counts'] = a.X.copy()
    if 'total_counts' not in a.obs:
        a.obs['total_counts'] = np.asarray(a.layers['counts'].sum(axis=1)).ravel()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=N_TOP_GENES, flavor='seurat_v3', layer='counts')
    a = a[:, a.var['highly_variable']].copy()
    sc.pp.regress_out(a, ['total_counts'])
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=50, svd_solver='arpack')
    return a


METHODS = {
    'cpm_log1p': method_cpm_log1p,
    'pearson': method_pearson,
}


def finalize_embedding(a):
    sc.pp.neighbors(a, n_pcs=N_PCS, use_rep='X_pca')
    sc.tl.umap(a, random_state=0)
    return a


def score_embedding(a, label_key='popv_prediction', subsample=10000):
    depth = a.obs['total_counts'].values
    u = a.obsm['X_umap']
    r_u = max(abs(spearmanr(depth, u[:, 0]).correlation),
              abs(spearmanr(depth, u[:, 1]).correlation))
    sil = np.nan
    try:
        from sklearn.metrics import silhouette_score
        rng = np.random.RandomState(0)
        idx = rng.choice(a.n_obs, min(subsample, a.n_obs), replace=False)
        lab = pd.Series(a.obs[label_key].astype(str).values[idx])
        keep = lab.isin(lab.value_counts()[lambda s: s >= 10].index)
        if keep.sum() > 50 and lab[keep].nunique() > 1:
            sil = float(silhouette_score(a.obsm['X_pca'][idx][keep.values], lab[keep].values))
    except Exception as e:
        print(f"    (silhouette skipped: {e})")
    return {'depth_umap_corr': float(r_u), 'label_silhouette': sil}


# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", default=None)
    ap.add_argument("--checkpoint", default=None, help="abs path to popv checkpoint h5ad")
    ap.add_argument("--out", default=None)
    ap.add_argument("--regress", action="store_true", help="add cpm+depth-regression arm")
    ap.add_argument("--subsample", type=int, default=0, help="subsample N cells for speed (0=all)")
    args = ap.parse_args()

    outputs = resolve_outputs(args.outputs)
    ckpt = (os.path.abspath(os.path.expanduser(args.checkpoint)) if args.checkpoint
            else os.path.join(outputs, "04_annotation", "all_pucks_popv_annotated.h5ad"))
    fig_dir = (os.path.abspath(os.path.expanduser(args.out)) if args.out
               else os.path.join(outputs, "TROUBLESHOOTING", "normalization_bakeoff"))

    methods = dict(METHODS)
    if args.regress:
        methods['cpm_regress'] = method_cpm_regress

    banner("NORMALIZATION BAKE-OFF")
    print("  [resolved paths] (all absolute)")
    for k, v in [("outputs", outputs), ("checkpoint", ckpt), ("figures", fig_dir)]:
        tag = "ok " if os.path.exists(v) else "MISSING"
        print(f"    {k:<11} [{tag}] {v}")
    print(f"  methods    : {list(methods)}")

    if not os.path.exists(ckpt):
        print("\nFATAL: checkpoint not found at the resolved path above. Pass it explicitly:")
        print("  python Step04_normalization_bakeoff.py \\")
        print(f"    --checkpoint {os.path.join(DEFAULT_OUTPUTS, '04_annotation', 'all_pucks_popv_annotated.h5ad')}")
        sys.exit(1)

    os.makedirs(fig_dir, exist_ok=True)

    base = sc.read_h5ad(ckpt)
    print(f"  loaded {base.n_obs:,} cells x {base.n_vars:,} genes")
    if not is_integer_counts(base.X):
        print("  FATAL: checkpoint X is not raw integer counts; this bake-off needs raw "
              "counts. Re-point --checkpoint at the pre-normalization object.")
        sys.exit(1)
    print("  X confirmed raw integer counts")

    if args.subsample and args.subsample < base.n_obs:
        rng = np.random.RandomState(0)
        idx = rng.choice(base.n_obs, args.subsample, replace=False)
        base = base[idx].copy()
        print(f"  subsampled to {base.n_obs:,} cells")

    base = ensure_depth(base)
    label_key = 'popv_prediction' if 'popv_prediction' in base.obs else None
    if label_key is None:
        print("  WARNING: popv_prediction not in checkpoint; silhouette/label panels disabled")

    # run each method
    results = {}
    metrics = {}
    for name, fn in methods.items():
        banner(f"METHOD: {name}")
        try:
            a = fn(base.copy())
            a.obs['total_counts'] = base.obs['total_counts'].values
            a.obs['log10_counts'] = base.obs['log10_counts'].values
            if label_key:
                a.obs[label_key] = base.obs[label_key].values
            a.obs['puck_id'] = base.obs['puck_id'].values if 'puck_id' in base.obs else 'NA'
            a = finalize_embedding(a)
            m = score_embedding(a, label_key=label_key) if label_key else {
                'depth_umap_corr': float(max(
                    abs(spearmanr(a.obs['total_counts'], a.obsm['X_umap'][:, 0]).correlation),
                    abs(spearmanr(a.obs['total_counts'], a.obsm['X_umap'][:, 1]).correlation))),
                'label_silhouette': np.nan}
            results[name] = a
            metrics[name] = m
            print(f"  depth_umap_corr = {m['depth_umap_corr']:.3f}  "
                  f"label_silhouette = {m['label_silhouette']:.3f}")
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f"  METHOD {name} FAILED: {e}")

    # grid figure
    banner("WRITING GRID FIGURE")
    cols = ([label_key] if label_key else []) + ['log10_counts', 'puck_id']
    nrow, ncol = len(results), len(cols)
    fig, axes = plt.subplots(nrow, ncol, figsize=(8 * ncol, 7 * nrow), squeeze=False)
    for r, (name, a) in enumerate(results.items()):
        for c, col in enumerate(cols):
            ax = axes[r][c]
            cmap = 'viridis' if col == 'log10_counts' else None
            legend = None if (col == label_key) else 'right margin'
            sc.pl.umap(a, color=col, ax=ax, show=False, frameon=False, size=3,
                       cmap=cmap, legend_loc=legend)
            tag = f" (depth r={metrics[name]['depth_umap_corr']:.2f})" if col == 'log10_counts' else ""
            ax.set_title(f"{name}: {col}{tag}", fontsize=12)
    plt.tight_layout()
    grid_path = os.path.join(fig_dir, "normalization_bakeoff_grid.png")
    plt.savefig(grid_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  grid -> {grid_path}")

    # metrics table
    banner("SUMMARY (lower depth_umap_corr = better; higher label_silhouette = better)")
    mdf = pd.DataFrame(metrics).T
    print(mdf.to_string())
    mdf.to_csv(os.path.join(fig_dir, "normalization_bakeoff_metrics.tsv"), sep='\t')
    print()
    if len(mdf):
        best_depth = mdf['depth_umap_corr'].idxmin()
        print(f"  least depth-driven embedding : {best_depth}")
        if mdf['label_silhouette'].notna().any():
            best_sil = mdf['label_silhouette'].idxmax()
            print(f"  best cell-type separation    : {best_sil}")
        print()
        print("  Reminder: if even the best method leaves a high depth_umap_corr and low")
        print("  label_silhouette, normalization is not the lever. Next move is per-bead")
        print("  popv_prediction as-is, or reference-mapping against the HNSCC scRNA-seq.")
    print(f"\n  figures: {fig_dir}")


if __name__ == "__main__":
    main()
