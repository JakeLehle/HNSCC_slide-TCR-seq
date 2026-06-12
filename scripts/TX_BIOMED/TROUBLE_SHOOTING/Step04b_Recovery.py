#!/usr/bin/env python3
"""
Step04b_Recovery.py

Recovery script for Step04 after the h5ad write crash.

popV completed successfully (~5.5 hours) and saved predictions to
  scripts/tmp/popv_output/predictions.csv
but the script crashed when writing the annotated h5ad due to a
nullable string dtype incompatibility between popV's in-memory
categoricals and anndata/h5py.

This script:
  1. Loads the QC h5ad (all 106,561 cells)
  2. Loads predictions.csv (99,341 annotated cells)
  3. Merges popV annotations into adata
  4. Applies revert_nullable_strings as a safety net
  5. Saves the popV-annotated checkpoint
  6. Runs post-annotation processing (normalize, HVG, PCA, UMAP, Leiden)
  7. Runs consensus annotation (ClusterCatcher weighted scoring)
  8. Generates all plots (UMAP, spatial, stacked bar)
  9. Saves the final annotated h5ad

Run in sc_pre env. Can run on login node if memory is available,
or via SLURM with the Step04 wrapper.

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spatial_config import (
    DIR_03_QC, DIR_04_ANNOTATION,
    PUCK_NAMES, N_HVG, N_PCS, LEIDEN_RESOLUTION, RANDOM_SEED,
    FONT_SIZE, DPI, COLOR_PUCK_29, COLOR_PUCK_37, COLOR_PUCK_40,
    banner, log, ensure_dir, save_fig,
)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe


# =============================================================================
# PATHS
# =============================================================================

QC_H5AD = os.path.join(DIR_03_QC, "all_pucks_merged_QC.h5ad")
PREDICTIONS_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "tmp", "popv_output", "predictions.csv"
)
OUTPUT_DIR = ensure_dir(DIR_04_ANNOTATION)
FIG_DIR = ensure_dir(os.path.join(DIR_04_ANNOTATION, "figures"))


# =============================================================================
# NULLABLE STRING FIX
# =============================================================================

def revert_nullable_strings(adata):
    """
    Revert nullable/arrow string dtypes to object dtype for h5ad compatibility.

    popV and other scverse tools can create categorical columns whose
    underlying categories use pandas' nullable StringDtype. anndata 0.12.x
    cannot serialize these to h5ad. This function forces all string-like
    columns back to object dtype.

    Reference: https://github.com/scverse/anndata/issues/2221
    """
    pd.options.future.infer_string = False

    adata.obs.index = adata.obs.index.astype(object)
    adata.var.index = adata.var.index.astype(object)
    if adata.raw is not None:
        adata.raw.var.index = adata.raw.var.index.astype(object)

    for frame_name, frame in [('obs', adata.obs), ('var', adata.var)]:
        str_cols = frame.select_dtypes(include=['category', 'string[python]']).columns
        for col in str_cols:
            frame[col] = pd.Series(frame[col], dtype="object")
            if frame[col].nunique() < frame[col].notna().sum():
                frame[col] = frame[col].astype("category")

    return adata


def safe_write_h5ad(adata, path, compression='gzip'):
    """Write h5ad with nullable string safety net."""
    adata = revert_nullable_strings(adata)
    adata.write_h5ad(path, compression=compression)
    size_mb = os.path.getsize(path) / 1e6
    log(f"  Saved: {path} ({size_mb:.1f} MB)")


# =============================================================================
# STEP 1: LOAD AND MERGE
# =============================================================================

def load_and_merge():
    """Load QC h5ad and merge popV predictions from CSV."""
    banner("STEP 1: LOAD QC DATA + POPV PREDICTIONS")

    # Load QC h5ad
    log(f"  Loading: {QC_H5AD}")
    adata = sc.read_h5ad(QC_H5AD)
    log(f"  QC adata: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

    # Load predictions
    log(f"  Loading: {PREDICTIONS_CSV}")
    preds = pd.read_csv(PREDICTIONS_CSV, index_col=0)
    log(f"  Predictions: {preds.shape[0]:,} cells x {preds.shape[1]} columns")
    log(f"  Columns: {list(preds.columns)}")

    # Align
    common = adata.obs_names.intersection(preds.index)
    excluded = adata.n_obs - len(common)
    log(f"  Common barcodes: {len(common):,}")
    log(f"  Excluded (low expression): {excluded:,}")

    # Subset to annotated cells
    adata = adata[common].copy()

    # Add prediction columns as categoricals
    for col in preds.columns:
        vals = preds.loc[adata.obs_names, col]
        if vals.dtype == 'int64' or vals.dtype == 'float64':
            adata.obs[col] = vals.values
        else:
            adata.obs[col] = pd.Categorical(vals)

    # Convert score columns to numeric (not categorical)
    for score_col in ['popv_prediction_score', 'popv_majority_vote_score']:
        if score_col in adata.obs.columns:
            adata.obs[score_col] = pd.to_numeric(
                adata.obs[score_col], errors='coerce'
            )

    log(f"  Merged adata: {adata.shape[0]:,} cells, "
        f"{len(adata.obs.columns)} obs columns")

    # Report cell type distribution
    if 'popv_prediction' in adata.obs.columns:
        n_types = adata.obs['popv_prediction'].nunique()
        log(f"\n  Cell types detected: {n_types}")
        log(f"  Top 10:")
        for ct, count in adata.obs['popv_prediction'].value_counts().head(10).items():
            pct = 100 * count / adata.n_obs
            log(f"    {ct}: {count:,} ({pct:.1f}%)")

    return adata


# =============================================================================
# STEP 2: POST-ANNOTATION PROCESSING
# =============================================================================

def post_annotation_processing(adata):
    """
    Normalize, HVG, PCA, UMAP, Leiden.

    Follows ClusterCatcher flow with spatial adaptations:
      - Explicit HVG selection (batch-aware by puck_id)
      - PCA on HVGs before neighbors (better for 100K+ cells)

    Stores raw counts in adata.layers['counts'] before normalization.
    """
    banner("STEP 2: POST-ANNOTATION PROCESSING")

    # Store raw counts as layer
    log(f"  Storing raw counts in layers['counts']...")
    adata.layers['counts'] = adata.X.copy()

    # Normalize
    log(f"  Normalizing (CPM + log1p)...")
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    # HVG selection (batch-aware)
    log(f"  Selecting highly variable genes...")
    try:
        sc.pp.highly_variable_genes(
            adata, min_mean=0.0125, max_mean=3, min_disp=0.5,
            batch_key='puck_id',
        )
        n_hvg = adata.var['highly_variable'].sum()
        log(f"    Batch-aware HVGs: {n_hvg:,}")
    except Exception as e:
        log(f"    Batch-aware HVG failed ({e}), using seurat_v3 top {N_HVG}...")
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=N_HVG)
        n_hvg = adata.var['highly_variable'].sum()
        log(f"    Seurat v3 HVGs: {n_hvg:,}")

    if n_hvg < 100:
        log(f"    Too few HVGs ({n_hvg}), forcing top {N_HVG}...")
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=N_HVG)

    # PCA on HVGs
    log(f"  Computing PCA on HVGs...")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    n_comps = min(50, adata_hvg.n_vars - 1, adata_hvg.n_obs - 1)
    sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=n_comps)
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns['pca'] = adata_hvg.uns['pca']
    log(f"    PCA components: {n_comps}")

    # Neighbors
    n_pcs = min(N_PCS, n_comps)
    log(f"  Computing neighbors (n_pcs={n_pcs})...")
    sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca')

    # UMAP
    log(f"  Computing UMAP...")
    sc.tl.umap(adata, random_state=RANDOM_SEED)

    # Leiden clustering
    log(f"  Leiden clustering (resolution={LEIDEN_RESOLUTION})...")
    sc.tl.leiden(adata, key_added='clusters',
                 resolution=LEIDEN_RESOLUTION, random_state=RANDOM_SEED)
    n_clusters = adata.obs['clusters'].nunique()
    log(f"    Found {n_clusters} clusters")

    return adata


# =============================================================================
# STEP 3: CONSENSUS ANNOTATION (ClusterCatcher weighted scoring)
# =============================================================================

def assign_cluster_based_annotation(adata):
    """
    Assign final cell type labels using cluster-level weighted scoring.

    Follows ClusterCatcher's method:
      1. Min-max normalize popV scores within each Leiden cluster
      2. Linear weight by normalized score
      3. Weighted score = prediction_score x weight
      4. For each cluster, the cell type with the highest summed
         weighted score becomes the consensus label
    """
    banner("STEP 3: CONSENSUS ANNOTATION")

    prediction_col = 'popv_prediction'
    score_col = 'popv_prediction_score'
    cluster_col = 'clusters'

    if prediction_col not in adata.obs.columns:
        log(f"  WARNING: {prediction_col} not found, skipping consensus")
        return adata

    if score_col not in adata.obs.columns:
        log(f"  WARNING: {score_col} not found, skipping consensus")
        return adata

    df = adata.obs[[cluster_col, prediction_col, score_col]].copy()
    df[score_col] = pd.to_numeric(df[score_col], errors='coerce').fillna(0)

    # Min-max normalize scores within each cluster
    def min_max_normalize(x):
        r = x.max() - x.min()
        if r == 0:
            return pd.Series(np.ones(len(x)) / len(x), index=x.index)
        return (x - x.min()) / r

    df['normalized_score'] = df.groupby(cluster_col)[score_col].transform(
        min_max_normalize
    )

    # Linear weights
    def linear_weights(x):
        total = x.sum()
        if total == 0:
            return pd.Series(np.ones(len(x)) / len(x), index=x.index)
        return x / total

    df['weight'] = df.groupby(cluster_col)['normalized_score'].transform(
        linear_weights
    )

    # Weighted score
    df['weighted_score'] = df[score_col] * df['weight']

    # For each cluster, get dominant cell type
    cluster_type_scores = df.groupby(
        [cluster_col, prediction_col]
    )['weighted_score'].sum().reset_index()

    dominant_types = cluster_type_scores.loc[
        cluster_type_scores.groupby(cluster_col)['weighted_score'].idxmax()
    ].set_index(cluster_col)[prediction_col]

    # Map back
    adata.obs['final_annotation'] = adata.obs[cluster_col].map(dominant_types)
    adata.obs['final_annotation'] = adata.obs['final_annotation'].astype('category')

    # Report
    n_types = adata.obs['final_annotation'].nunique()
    log(f"  Final annotation: {n_types} cell types across "
        f"{adata.obs[cluster_col].nunique()} clusters")
    log(f"\n  Cell type distribution:")
    for ct, count in adata.obs['final_annotation'].value_counts().items():
        pct = 100 * count / adata.n_obs
        log(f"    {ct}: {count:,} ({pct:.1f}%)")

    return adata


# =============================================================================
# STEP 4: VISUALIZATION
# =============================================================================

def generate_plots(adata, fig_dir):
    """Generate UMAP, spatial, and composition plots."""
    banner("STEP 4: GENERATING PLOTS")

    # --- UMAP: popV prediction ---
    if 'popv_prediction' in adata.obs.columns:
        log(f"  UMAP: popv_prediction")
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color='popv_prediction', ax=ax, show=False,
                    legend_loc='right margin', frameon=False, size=3)
        plt.tight_layout()
        save_fig(fig, 'UMAP_popv_prediction', fig_dir)

    # --- UMAP: popV score ---
    if 'popv_prediction_score' in adata.obs.columns:
        log(f"  UMAP: popv_prediction_score")
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color='popv_prediction_score', ax=ax, show=False,
                    cmap='magma', frameon=False, size=3)
        plt.tight_layout()
        save_fig(fig, 'UMAP_popv_score', fig_dir)

    # --- UMAP: final annotation ---
    if 'final_annotation' in adata.obs.columns:
        log(f"  UMAP: final_annotation")
        fig, ax = plt.subplots(figsize=(14, 10))
        sc.pl.umap(adata, color='final_annotation', ax=ax, show=False,
                    legend_loc='right margin', frameon=False, size=3)
        plt.tight_layout()
        save_fig(fig, 'UMAP_final_annotation', fig_dir)

    # --- UMAP: clusters ---
    log(f"  UMAP: clusters")
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='clusters', ax=ax, show=False,
                legend_loc='right margin', frameon=False, size=3)
    plt.tight_layout()
    save_fig(fig, 'UMAP_clusters', fig_dir)

    # --- UMAP: puck_id ---
    log(f"  UMAP: puck_id")
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='puck_id', ax=ax, show=False,
                legend_loc='right margin', frameon=False, size=3)
    plt.tight_layout()
    save_fig(fig, 'UMAP_puck_id', fig_dir)

    # --- Stacked bar: cell type per cluster ---
    if 'popv_prediction' in adata.obs.columns:
        log(f"  Stacked bar: cluster composition")
        plot_stacked_bar(adata, fig_dir)

    # --- Spatial scatter per puck ---
    annotation_col = (
        'final_annotation' if 'final_annotation' in adata.obs.columns
        else 'clusters'
    )
    for puck_name in PUCK_NAMES:
        mask = adata.obs['puck_id'] == puck_name
        if mask.sum() == 0:
            continue
        log(f"  Spatial: {puck_name}")
        plot_spatial_scatter(
            adata[mask].copy(), puck_name, annotation_col, fig_dir
        )
        plot_spatial_scatter(
            adata[mask].copy(), puck_name, 'clusters', fig_dir
        )


def plot_stacked_bar(adata, fig_dir):
    """Stacked bar of cell type composition per Leiden cluster."""
    ct = pd.crosstab(
        adata.obs['clusters'],
        adata.obs['popv_prediction'],
        normalize='index'
    )
    fig, ax = plt.subplots(figsize=(20, 8))
    ct.plot(kind='bar', stacked=True, ax=ax, width=0.85,
            legend=False, edgecolor='none')
    ax.set_xlabel('Leiden cluster', fontsize=FONT_SIZE)
    ax.set_ylabel('Fraction', fontsize=FONT_SIZE)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
              fontsize=FONT_SIZE - 10, ncol=1)
    plt.tight_layout()
    save_fig(fig, 'stacked_bar_cluster_celltypes', fig_dir)


def plot_spatial_scatter(adata_sub, puck_name, color_col, fig_dir):
    """Spatial scatter plot for one puck colored by annotation or cluster."""
    if 'x_coord' not in adata_sub.obs.columns:
        return

    fig, ax = plt.subplots(figsize=(10, 10))
    x = adata_sub.obs['x_coord'].values
    y = adata_sub.obs['y_coord'].values

    if color_col in adata_sub.obs.columns:
        cats = adata_sub.obs[color_col].astype('category')
        codes = cats.cat.codes.values
        n_cats = cats.cat.categories.size
        cmap = plt.cm.get_cmap('tab20', max(n_cats, 20))
        colors = [cmap(c % 20) for c in codes]
        ax.scatter(x, y, c=colors, s=1, alpha=0.6, rasterized=True)
    else:
        ax.scatter(x, y, s=1, alpha=0.6, c='#999999', rasterized=True)

    ax.set_aspect('equal')
    ax.set_title(f'{puck_name}: {color_col}', fontsize=FONT_SIZE)
    ax.axis('off')
    plt.tight_layout()
    save_fig(fig, f'spatial_{puck_name}_{color_col}', fig_dir)


# =============================================================================
# MAIN
# =============================================================================

def main():
    banner("STEP 04b RECOVERY: CELL TYPE ANNOTATION (from predictions.csv)")
    log(f"  Skipping popV (already completed, ~5.5 hours saved)")
    log(f"  Using saved predictions: {PREDICTIONS_CSV}")

    # Verify predictions exist
    if not os.path.exists(PREDICTIONS_CSV):
        log(f"  ERROR: Predictions file not found: {PREDICTIONS_CSV}")
        log(f"  Cannot recover. Re-run Step04 with the nullable string fix.")
        sys.exit(1)

    # Step 1: Load and merge
    adata = load_and_merge()

    # Save popV checkpoint (the step that originally crashed)
    popv_path = os.path.join(OUTPUT_DIR, "adata_popv_annotated.h5ad")
    log(f"\n  Saving popV checkpoint...")
    safe_write_h5ad(adata, popv_path)

    # Step 2: Post-annotation processing
    adata = post_annotation_processing(adata)

    # Step 3: Consensus annotation
    adata = assign_cluster_based_annotation(adata)

    # Step 4: Plots
    generate_plots(adata, FIG_DIR)

    # Final save
    banner("SAVING FINAL ANNOTATED ADATA")
    final_path = os.path.join(OUTPUT_DIR, "all_pucks_annotated.h5ad")
    safe_write_h5ad(adata, final_path)

    # Summary report
    summary_path = os.path.join(OUTPUT_DIR, "Step04b_recovery_summary.tsv")
    summary = {
        'total_cells_QC': 106561,
        'annotated_cells': adata.n_obs,
        'excluded_low_expression': 106561 - adata.n_obs,
        'genes': adata.n_vars,
        'n_clusters': adata.obs['clusters'].nunique(),
        'n_cell_types_popv': adata.obs['popv_prediction'].nunique()
            if 'popv_prediction' in adata.obs.columns else 0,
        'n_cell_types_final': adata.obs['final_annotation'].nunique()
            if 'final_annotation' in adata.obs.columns else 0,
    }
    pd.DataFrame([summary]).to_csv(summary_path, sep='\t', index=False)
    log(f"  Summary: {summary_path}")

    banner("STEP 04b RECOVERY COMPLETE")
    for k, v in summary.items():
        log(f"  {k}: {v}")


if __name__ == "__main__":
    main()
