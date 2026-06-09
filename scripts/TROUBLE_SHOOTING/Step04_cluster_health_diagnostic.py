#!/usr/bin/env python3
"""
Step04_cluster_health_diagnostic.py

Investigate the "spiral horn / silly-string" UMAP and the coarse cell-type
calls, using the adata objects already on disk (no popV rerun). Tests three
independent hypotheses:

  H1  MERGE CORRUPTION  -- did merging the 3 pucks drop or scramble the
                           gene-expression matrix? (per-puck QC vs merged)
  H2  TOO-FEW-GENES /    -- did the batch-aware HVG step collapse to a thin
      DEPTH-DRIVEN          gene set, or is the manifold just library-size
                            driven? (HVG count + depth-colored UMAP)
  H3  GENE-NAME / MODEL  -- were var_names the right symbols with enough
      OVERLAP               overlap for popV? (markers + Ensembl check)

Plus a direct test of the "keep pucks separate" escape hatch: recompute an
independent UMAP for each puck from its own QC object and eyeball the shape.

Reads (relative to data/outputs):
  02_anndata/<puck>.h5ad                  raw, pre-QC
  03_qc/<puck>_QC.h5ad                    per-puck QC
  03_qc/all_pucks_merged_QC.h5ad          merged (popV input)
  04_annotation/all_pucks_annotated.h5ad  final (HVG, clusters, X_umap)

Writes diagnostic PNGs + a text report to:
  <outputs>/TROUBLESHOOTING/cluster_diagnostic/   (override with --out)

Run (sc_pre env):
  python Step04_cluster_health_diagnostic.py
  python Step04_cluster_health_diagnostic.py --outputs /abs/data/outputs
  python Step04_cluster_health_diagnostic.py --skip-puck-umaps   # faster

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
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

PUCKS = ['Puck_211214_29', 'Puck_211214_37', 'Puck_211214_40']

# Canonical markers for HNSCC + immune compartments. Used only to confirm the
# right gene symbols are present in var_names (H3); not for annotation.
MARKER_PANEL = {
    'basal_epithelial':   ['KRT5', 'KRT14', 'KRT15', 'TP63', 'KRT6A', 'KRT17'],
    'diff_epithelial':    ['KRT1', 'KRT10', 'IVL', 'SBSN', 'KRT13'],
    'epithelial_general': ['EPCAM', 'SFN', 'CDH1'],
    'T_cell':             ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD4', 'TRAC'],
    'B_cell':             ['CD19', 'MS4A1', 'CD79A', 'CD79B'],
    'myeloid_macrophage': ['CD68', 'LYZ', 'CD14', 'ITGAM', 'CSF1R', 'C1QA'],
    'fibroblast':         ['COL1A1', 'COL1A2', 'DCN', 'PDGFRA', 'LUM'],
    'endothelial':        ['PECAM1', 'VWF', 'CDH5', 'CLDN5'],
}


# -----------------------------------------------------------------------------
def banner(msg):
    print("\n" + "=" * 72)
    print(msg)
    print("=" * 72)


def _autodetect_outputs():
    d = THIS_DIR
    for _ in range(7):
        cand = os.path.join(d, "data", "outputs")
        if os.path.isdir(os.path.join(cand, "03_qc")):
            return cand
        if os.path.isdir(os.path.join(d, "03_qc")) and os.path.basename(d) == "outputs":
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    # last-resort known path
    return "/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs"


def is_integer_counts(X, n=2000):
    """Heuristic: are the first n nonzero values integers (i.e., raw counts)?"""
    if sparse.issparse(X):
        vals = X.data[:n]
    else:
        vals = np.asarray(X).ravel()
        vals = vals[vals != 0][:n]
    if len(vals) == 0:
        return True, 0.0
    frac_int = float(np.mean(np.isclose(vals, np.round(vals))))
    return frac_int > 0.999, frac_int


def summarize(adata, label):
    """Print a compact inventory and return a dict of key facts."""
    X = adata.X
    nnz = X.nnz if sparse.issparse(X) else int(np.count_nonzero(X))
    total = float(X.sum())
    intish, fi = is_integer_counts(X)
    # depth
    if sparse.issparse(X):
        per_cell_counts = np.asarray(X.sum(axis=1)).ravel()
        per_cell_genes = np.asarray((X > 0).sum(axis=1)).ravel()
    else:
        per_cell_counts = X.sum(axis=1)
        per_cell_genes = (X > 0).sum(axis=1)
    n_ensembl = int(adata.var_names.to_series().astype(str).str.startswith('ENSG').sum())

    print(f"\n[{label}]")
    print(f"    shape            : {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"    X integer counts : {intish} (int-fraction={fi:.3f})")
    print(f"    sparsity         : {100*(1 - nnz/(adata.n_obs*adata.n_vars)):.2f}% zeros")
    print(f"    median genes/cell: {np.median(per_cell_genes):.0f}   "
          f"median counts/cell: {np.median(per_cell_counts):.0f}")
    print(f"    var_names Ensembl: {n_ensembl}/{adata.n_vars}  "
          f"({'SYMBOLS ok' if n_ensembl == 0 else 'CONTAINS ENSEMBL IDs'})")
    has_feat = 'feature_name' in adata.var.columns
    if has_feat:
        match = (adata.var['feature_name'].astype(str).values == adata.var_names.astype(str).values).mean()
        print(f"    var['feature_name'] present, matches var_names: {100*match:.1f}%")
    print(f"    obs columns      : {list(adata.obs.columns)[:12]}"
          f"{' ...' if adata.obs.shape[1] > 12 else ''}")

    return {
        'label': label, 'n_obs': adata.n_obs, 'n_vars': adata.n_vars,
        'int_counts': intish, 'median_genes': float(np.median(per_cell_genes)),
        'median_counts': float(np.median(per_cell_counts)),
        'n_ensembl': n_ensembl, 'total': total,
        'var_names': set(adata.var_names.astype(str)),
        'gene_totals': pd.Series(
            np.asarray(X.sum(axis=0)).ravel(), index=adata.var_names.astype(str)
        ),
    }


# -----------------------------------------------------------------------------
def check_marker_overlap(var_names, label):
    print(f"\n[H3] Marker-panel presence in {label}")
    vn = set(var_names)
    verdict_ok = True
    for cat, genes in MARKER_PANEL.items():
        present = [g for g in genes if g in vn]
        missing = [g for g in genes if g not in vn]
        frac = len(present) / len(genes)
        flag = "" if frac >= 0.5 else "  <-- LOW"
        if frac < 0.5:
            verdict_ok = False
        print(f"    {cat:<20} {len(present)}/{len(genes)} present{flag}"
              f"{'   missing: ' + ','.join(missing) if missing else ''}")
    return verdict_ok


# -----------------------------------------------------------------------------
def recompute_puck_umap(qc_path, puck, fig_dir, n_top=2000):
    """Independent UMAP for a single puck from its own QC object."""
    a = sc.read_h5ad(qc_path)
    n0 = a.n_obs
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e6)
    sc.pp.log1p(a)
    try:
        sc.pp.highly_variable_genes(a, flavor='seurat_v3', n_top_genes=n_top, layer='counts')
    except Exception:
        sc.pp.highly_variable_genes(a, min_mean=0.0125, max_mean=3, min_disp=0.5)
    n_hvg = int(a.var['highly_variable'].sum())
    ah = a[:, a.var['highly_variable']].copy()
    n_comps = min(50, ah.n_vars - 1, ah.n_obs - 1)
    sc.tl.pca(ah, svd_solver='arpack', n_comps=n_comps)
    a.obsm['X_pca'] = ah.obsm['X_pca']
    sc.pp.neighbors(a, n_pcs=min(40, n_comps), use_rep='X_pca')
    sc.tl.umap(a, random_state=0)
    sc.tl.leiden(a, key_added='leiden_indep', resolution=0.5, random_state=0)

    counts = np.asarray(a.layers['counts'].sum(axis=1)).ravel()
    a.obs['log10_counts'] = np.log10(counts + 1)

    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    sc.pl.umap(a, color='leiden_indep', ax=axes[0], show=False, frameon=False, size=4)
    axes[0].set_title(f'{puck}: independent Leiden ({n_hvg} HVGs, {a.n_obs:,} cells)')
    sc.pl.umap(a, color='log10_counts', ax=axes[1], show=False, frameon=False,
               size=4, cmap='viridis')
    axes[1].set_title(f'{puck}: depth (log10 counts)')
    plt.tight_layout()
    out = os.path.join(fig_dir, f'INDEP_UMAP_{puck}.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    {puck}: {n0:,} cells, {n_hvg} HVGs -> {out}")
    return {'puck': puck, 'n_hvg': n_hvg, 'n_obs': n0}


def plot_combined_depth_umap(annotated, fig_dir):
    """Color the EXISTING combined UMAP by depth to test the depth hypothesis."""
    if 'X_umap' not in annotated.obsm:
        print("    (annotated object has no X_umap; skipping)")
        return
    a = annotated
    # depth from counts layer if present, else from raw
    if 'counts' in a.layers:
        counts = np.asarray(a.layers['counts'].sum(axis=1)).ravel()
    else:
        counts = np.asarray(a.X.sum(axis=1)).ravel()
    a.obs['log10_counts'] = np.log10(counts + 1)

    color_cols = ['log10_counts']
    if 'puck_id' in a.obs:
        color_cols.append('puck_id')
    if 'clusters' in a.obs:
        color_cols.append('clusters')

    fig, axes = plt.subplots(1, len(color_cols), figsize=(9 * len(color_cols), 8))
    if len(color_cols) == 1:
        axes = [axes]
    for ax, c in zip(axes, color_cols):
        cmap = 'viridis' if c == 'log10_counts' else None
        sc.pl.umap(a, color=c, ax=ax, show=False, frameon=False, size=3, cmap=cmap)
        ax.set_title(f'combined UMAP: {c}')
    plt.tight_layout()
    out = os.path.join(fig_dir, 'COMBINED_UMAP_depth_check.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    combined depth-colored UMAP -> {out}")


# -----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", default=None, help="absolute path to data/outputs")
    ap.add_argument("--out", default=None, help="diagnostic output dir")
    ap.add_argument("--skip-puck-umaps", action="store_true")
    args = ap.parse_args()

    outputs = os.path.abspath(args.outputs) if args.outputs else _autodetect_outputs()
    fig_dir = os.path.abspath(args.out) if args.out else os.path.join(
        outputs, "TROUBLESHOOTING", "cluster_diagnostic")
    os.makedirs(fig_dir, exist_ok=True)

    banner("CLUSTER HEALTH DIAGNOSTIC")
    print(f"  outputs : {outputs}")
    print(f"  figures : {fig_dir}")

    p_raw = {p: os.path.join(outputs, "02_anndata", f"{p}.h5ad") for p in PUCKS}
    p_qc = {p: os.path.join(outputs, "03_qc", f"{p}_QC.h5ad") for p in PUCKS}
    merged_path = os.path.join(outputs, "03_qc", "all_pucks_merged_QC.h5ad")
    annot_path = os.path.join(outputs, "04_annotation", "all_pucks_annotated.h5ad")

    # =====================================================================
    banner("[STEP 1] INVENTORY")
    # =====================================================================
    facts = {}
    for p in PUCKS:
        if os.path.exists(p_qc[p]):
            facts[p] = summarize(sc.read_h5ad(p_qc[p]), f"QC {p}")
        else:
            print(f"  MISSING: {p_qc[p]}")
    merged = sc.read_h5ad(merged_path) if os.path.exists(merged_path) else None
    if merged is not None:
        facts['merged'] = summarize(merged, "MERGED (popV input)")

    # =====================================================================
    banner("[STEP 2 / H1] MERGE INTEGRITY")
    # =====================================================================
    h1_ok = True
    if merged is not None and all(p in facts for p in PUCKS):
        sum_cells = sum(facts[p]['n_obs'] for p in PUCKS)
        print(f"  cell count: sum(per-puck QC)={sum_cells:,}  merged={facts['merged']['n_obs']:,}"
              f"   {'MATCH' if sum_cells == facts['merged']['n_obs'] else 'MISMATCH <--'}")
        if sum_cells != facts['merged']['n_obs']:
            h1_ok = False

        # gene set relationships
        sets = [facts[p]['var_names'] for p in PUCKS]
        inter = set.intersection(*sets)
        union = set.union(*sets)
        mset = facts['merged']['var_names']
        print(f"  genes: per-puck intersection={len(inter):,}  union={len(union):,}  "
              f"merged={len(mset):,}")
        print(f"         merged == union?        {mset == union}")
        print(f"         merged == intersection? {mset == inter}")
        missing_from_merged = union - mset
        if missing_from_merged:
            print(f"         {len(missing_from_merged)} genes present per-puck but DROPPED in merge")
            if len(missing_from_merged) > len(union) * 0.05:
                h1_ok = False

        # per-gene total-count reconciliation (naming-robust integrity check)
        if 'puck_id' in merged.obs.columns:
            print(f"  per-gene count reconciliation (per-puck QC vs merged subset):")
            for p in PUCKS:
                sub = merged[merged.obs['puck_id'].astype(str) == p]
                if sub.n_obs == 0:
                    print(f"    {p}: 0 cells in merged subset  <-- MISMATCH")
                    h1_ok = False
                    continue
                msub = pd.Series(np.asarray(sub.X.sum(axis=0)).ravel(),
                                 index=sub.var_names.astype(str))
                shared = facts[p]['gene_totals'].index.intersection(msub.index)
                a_tot = facts[p]['gene_totals'].reindex(shared).fillna(0).values
                b_tot = msub.reindex(shared).fillna(0).values
                denom = np.maximum(a_tot, 1)
                max_rel = float(np.max(np.abs(a_tot - b_tot) / denom)) if len(shared) else 0.0
                agree = max_rel < 1e-6
                print(f"    {p}: {sub.n_obs:,} cells, {len(shared):,} shared genes, "
                      f"max rel-diff={max_rel:.2e}  {'OK' if agree else 'MISMATCH <--'}")
                if not agree:
                    h1_ok = False

        # all-zero genes / cells + duplicate names in merged
        Xg = np.asarray(merged.X.sum(axis=0)).ravel()
        Xc = np.asarray(merged.X.sum(axis=1)).ravel()
        print(f"  merged all-zero genes: {(Xg == 0).sum():,}   all-zero cells: {(Xc == 0).sum():,}")
        dup_v = merged.var_names.duplicated().sum()
        dup_o = merged.obs_names.duplicated().sum()
        print(f"  merged duplicate var_names: {dup_v}   duplicate obs_names: {dup_o}")
        if dup_v or dup_o:
            h1_ok = False
    else:
        print("  cannot run (missing per-puck QC or merged object)")
        h1_ok = None

    # =====================================================================
    banner("[STEP 3 / H3] GENE NAMES & MODEL OVERLAP")
    # =====================================================================
    h3_ok = None
    if merged is not None:
        intish, fi = is_integer_counts(merged.X)
        print(f"  merged X is raw counts (popV requirement): {intish} (int-frac={fi:.3f})"
              f"{'' if intish else '   <-- popV needs raw counts!'}")
        h3_ok = check_marker_overlap(facts['merged']['var_names'], "MERGED (popV input)")
        if facts['merged']['n_ensembl'] > 0:
            print(f"  WARNING: merged var_names contain {facts['merged']['n_ensembl']} Ensembl IDs; "
                  f"if popV matched on these, symbol overlap with the model is degraded")
            h3_ok = False

    # =====================================================================
    banner("[STEP 4 / H2] HVG COUNT + DEPTH (production combined embedding)")
    # =====================================================================
    h2_note = ""
    if os.path.exists(annot_path):
        annotated = sc.read_h5ad(annot_path)
        n_hvg_prod = int(annotated.var['highly_variable'].sum()) if 'highly_variable' in annotated.var else None
        print(f"  production HVGs that fed PCA: {n_hvg_prod}")
        if n_hvg_prod is not None and n_hvg_prod < 500:
            h2_note = f"LOW HVG COUNT ({n_hvg_prod}) -- thin feature set is a likely cause of the stringy UMAP"
            print(f"    <-- {h2_note}")
        if 'highly_variable_nbatches' in annotated.var:
            nb = annotated.var['highly_variable_nbatches']
            print(f"    HVG-in-all-3-pucks: {(nb >= 3).sum():,}   in>=2: {(nb >= 2).sum():,}   in>=1: {(nb >= 1).sum():,}")
        print(f"  median genes/cell (merged): {facts['merged']['median_genes']:.0f}"
              f"{'   <-- low depth; embedding may be library-size driven' if facts['merged']['median_genes'] < 400 else ''}")
        plot_combined_depth_umap(annotated, fig_dir)
    else:
        print(f"  MISSING: {annot_path}")

    # =====================================================================
    banner("[STEP 5] INDEPENDENT PER-PUCK UMAPS (the 'keep separate' test)")
    # =====================================================================
    puck_umap_info = []
    if args.skip_puck_umaps:
        print("  skipped (--skip-puck-umaps)")
    else:
        for p in PUCKS:
            if os.path.exists(p_qc[p]):
                puck_umap_info.append(recompute_puck_umap(p_qc[p], p, fig_dir))
            else:
                print(f"  MISSING: {p_qc[p]}")

    # =====================================================================
    banner("VERDICT")
    # =====================================================================
    def verdict(flag):
        return "PASS" if flag is True else ("FAIL" if flag is False else "N/A")
    print(f"  H1 merge integrity     : {verdict(h1_ok)}")
    print(f"  H3 gene names / overlap: {verdict(h3_ok)}")
    print(f"  H2 thin-genes / depth  : {h2_note or 'see HVG count + depth UMAP above'}")
    print()
    print("  Interpretation:")
    print("   - If H1 FAIL: the merge dropped/scrambled data; fix Step03 before anything else.")
    print("   - If H1 PASS but production HVGs are few (H2): the batch-aware HVG intersection")
    print("     collapsed; widen it (min_disp lower, or seurat_v3 top-N without batch_key) and")
    print("     re-embed. Check whether the combined depth UMAP shows a gradient along the snake.")
    print("   - If the INDEP_UMAP_*.png per-puck plots look normal (blobs, not strings), the")
    print("     problem is specific to combination -> per-puck analysis is a valid fallback.")
    print("   - If H3 FAIL: wrong var column / Ensembl IDs reached popV -> fix gene_symbols")
    print("     before trusting cell types; re-run popV.")
    print(f"\n  figures: {fig_dir}")


if __name__ == "__main__":
    main()
