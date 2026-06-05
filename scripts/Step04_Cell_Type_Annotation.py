#!/usr/bin/env python3
"""
Step04_Cell_Type_Annotation.py

Automated cell type annotation and clustering for QC'd Slide-seq data.
Follows the ClusterCatcher scanpy_qc_annotation.py pipeline flow adapted
for spatial bead-based data.

Pipeline:
  1. Load merged QC h5ad (raw counts)
  2. popV annotation (Tabula Sapiens, on raw counts before normalization)
  3. Merge popV annotations back to original AnnData
  4. Normalize (CPM + log1p), store raw counts in layer
  5. HVG selection (batch-aware by puck_id)
  6. PCA, neighbors, UMAP
  7. Leiden clustering
  8. Cluster-based consensus annotation (weighted scoring, ClusterCatcher method)
  9. Visualization (UMAP + spatial)
 10. Save annotated h5ad

Differences from ClusterCatcher (Slide-seq adaptations):
  - No Scrublet doublet detection (beads do not capture doublets)
  - No BBKNN (batch handled via HVG selection with puck_id key)
  - Explicit HVG + PCA before neighbors (better for 100K+ cells x 19K genes)
  - Spatial scatter plots in addition to UMAP

Reads from:  data/outputs/03_qc/
Writes to:   data/outputs/04_annotation/

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
import popv
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


# =========================================================================
# POPV ANNOTATION (follows ClusterCatcher flow)
# =========================================================================

def run_popv_annotation(adata, cache_dir):
    """
    Run popV cell type annotation on raw counts.

    Follows the ClusterCatcher scanpy_qc_annotation.py workflow:
      1. Save original obs columns
      2. Ensure gene_symbol / feature_name column exists
      3. Pull Tabula Sapiens HubModel
      4. Annotate with popV (on raw, unnormalized counts)
      5. Filter to query cells (remove reference cells)
      6. Merge popV columns back to original AnnData

    Parameters
    ----------
    adata : AnnData
        QC'd merged AnnData with raw counts in X.
    cache_dir : str
        Directory for caching the popV model.

    Returns
    -------
    adata : AnnData
        Original AnnData with popV annotation columns added to obs.
    """
    banner("POPV CELL TYPE ANNOTATION")

    huggingface_repo = "popV/tabula_sapiens_All_Cells"
    query_batch_key = "puck_id"
    gene_symbol_key = "feature_name"

    log(f"  Model:     {huggingface_repo}")
    log(f"  Batch key: {query_batch_key}")
    log(f"  Cells:     {adata.n_obs:,}")
    log(f"  Genes:     {adata.n_vars:,}")

    # Ensure gene symbol column exists (ClusterCatcher pattern)
    if gene_symbol_key not in adata.var.columns:
        adata.var[gene_symbol_key] = adata.var_names
        log(f"  Created var['{gene_symbol_key}'] from var_names")

    # Save original obs columns for clean merge
    original_obs_columns = list(adata.obs.columns)

    # Pull model
    log(f"  Pulling model from HuggingFace...")
    ensure_dir(cache_dir)
    hmo = popv.hub.HubModel.pull_from_huggingface_hub(
        huggingface_repo, cache_dir=cache_dir
    )

    # Annotate (on raw counts, no normalization)
    log(f"  Running popV annotation (this may take 30-60 minutes)...")
    adata_annotated = hmo.annotate_data(
        adata,
        query_batch_key=query_batch_key,
        prediction_mode="inference",
        gene_symbols=gene_symbol_key,
    )

    # Filter to query cells only (remove reference cells added by popV)
    if '_dataset' in adata_annotated.obs.columns:
        log(f"  Filtering to query cells only...")
        adata_annotated = adata_annotated[
            adata_annotated.obs['_dataset'] == 'query'
        ].copy()

    # Merge popV annotation columns back to original AnnData
    # (ClusterCatcher pattern: only add new columns, preserve original structure)
    log(f"  Merging popV annotations back to original AnnData...")
    popv_obs = adata_annotated.obs.drop(columns=original_obs_columns, errors='ignore')
    merged_df = pd.merge(
        adata.obs, popv_obs,
        left_index=True, right_index=True, how='inner'
    )

    adata = adata[adata.obs_names.isin(merged_df.index)].copy()
    merged_df = merged_df.loc[adata.obs_names]
    adata.obs = merged_df

    assert all(adata.obs_names == merged_df.index), "Index mismatch after popV merge"

    log(f"  Merged {len(popv_obs.columns)} popV columns, {adata.n_obs:,} cells retained")

    # Report cell type distribution
    if 'popv_prediction' in adata.obs.columns:
        n_types = adata.obs['popv_prediction'].nunique()
        log(f"  Cell types detected: {n_types}")
        log(f"  Top 10:")
        for ct, count in adata.obs['popv_prediction'].value_counts().head(10).items():
            pct = 100 * count / adata.n_obs
            log(f"    {ct}: {count:,} ({pct:.1f}%)")

    return adata


# =========================================================================
# POST-ANNOTATION PROCESSING
# =========================================================================

def post_annotation_processing(adata):
    """
    Normalize, embed, and cluster after popV annotation.

    Follows ClusterCatcher flow with two additions for spatial data:
      - Explicit HVG selection (batch-aware by puck_id)
      - PCA on HVGs before neighbors (better for 100K+ cells)

    Stores raw counts in adata.layers['counts'] before normalization.
    """
    banner("POST-ANNOTATION PROCESSING")

    # Store raw counts as layer (ClusterCatcher pattern)
    log(f"  Storing raw counts in layers['counts']...")
    adata.layers['counts'] = adata.X.copy()

    # Normalize
    log(f"  Normalizing (CPM + log1p)...")
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    # HVG selection (batch-aware, spatial adaptation)
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


# =========================================================================
# CONSENSUS ANNOTATION (ClusterCatcher weighted scoring)
# =========================================================================

def assign_cluster_based_annotation(adata):
    """
    Assign final cell type labels using cluster-level weighted scoring.

    Exactly follows ClusterCatcher's assign_cluster_based_annotation():
      1. Min-max normalize popV scores within each Leiden cluster
      2. Linear weight by normalized score
      3. Weighted score = prediction_score x weight
      4. Dominant type per cluster = highest aggregated weighted score
    """
    banner("CLUSTER-BASED CONSENSUS ANNOTATION")

    has_scores = 'popv_prediction_score' in adata.obs.columns

    if has_scores:
        df = adata.obs[['clusters', 'popv_prediction', 'popv_prediction_score']].copy()
        df['popv_prediction_score'] = pd.to_numeric(
            df['popv_prediction_score'], errors='coerce'
        )

        def min_max_normalize(x):
            x_min, x_range = x.min(), x.max() - x.min()
            return np.zeros(len(x)) if x_range == 0 else (x - x_min) / x_range

        def linear_weights(x):
            total = x.sum()
            return np.ones(len(x)) / len(x) if total == 0 else x / total

        df['normalized_score'] = df.groupby('clusters')[
            'popv_prediction_score'
        ].transform(min_max_normalize)
        df['weight'] = df.groupby('clusters')[
            'normalized_score'
        ].transform(linear_weights)
        df['weighted_score'] = df['popv_prediction_score'] * df['weight']

        agg = df.groupby(
            ['clusters', 'popv_prediction']
        )['weighted_score'].sum().reset_index()
        dominant = agg.loc[
            agg.groupby('clusters')['weighted_score'].idxmax()
        ].set_index('clusters')['popv_prediction']
    else:
        log(f"  No popv_prediction_score found, using majority voting")
        counts = adata.obs.groupby(
            ['clusters', 'popv_prediction']
        ).size().reset_index(name='count')
        dominant = counts.loc[
            counts.groupby('clusters')['count'].idxmax()
        ].set_index('clusters')['popv_prediction']

    adata.obs['final_annotation'] = adata.obs['clusters'].map(dominant)

    log(f"  Final cell type assignments:")
    for ct, count in adata.obs['final_annotation'].value_counts().items():
        pct = 100 * count / adata.n_obs
        log(f"    {ct}: {count:,} ({pct:.1f}%)")

    return adata


# =========================================================================
# VISUALIZATION
# =========================================================================

def generate_annotation_plots(adata, fig_dir):
    """Generate all annotation plots following ClusterCatcher output set."""
    banner("GENERATING ANNOTATION PLOTS")

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
        adata.obs['popv_prediction_score'] = pd.to_numeric(
            adata.obs['popv_prediction_score'], errors='coerce'
        )
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color='popv_prediction_score', ax=ax, show=False,
                    cmap='magma', frameon=False, size=3)
        plt.tight_layout()
        save_fig(fig, 'UMAP_popv_score', fig_dir)

    # --- UMAP: final annotation with labels ---
    if 'final_annotation' in adata.obs.columns:
        log(f"  UMAP: final_annotation")
        adata.obs['final_annotation'] = adata.obs['final_annotation'].astype('category')

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

    # --- Stacked bar: cell type composition per cluster ---
    if 'popv_prediction' in adata.obs.columns:
        log(f"  Stacked bar: cluster composition")
        plot_stacked_bar(adata, fig_dir)

    # --- Spatial plots per puck ---
    annotation_col = 'final_annotation' if 'final_annotation' in adata.obs.columns else 'clusters'
    for puck_name in PUCK_NAMES:
        mask = adata.obs['puck_id'] == puck_name
        if mask.sum() == 0:
            continue
        log(f"  Spatial: {puck_name}")
        plot_spatial_scatter(adata[mask].copy(), puck_name, annotation_col, fig_dir)
        plot_spatial_scatter(adata[mask].copy(), puck_name, 'clusters', fig_dir)


def plot_stacked_bar(adata, fig_dir):
    """Stacked bar of cell type composition per Leiden cluster."""
    cluster_counts = adata.obs['clusters'].value_counts()
    clusters = sorted(
        [c for c in cluster_counts.index if cluster_counts[c] >= 50],
        key=lambda x: int(x) if x.isdigit() else x
    )

    fig, ax = plt.subplots(figsize=(max(16, len(clusters) * 1.2), 7))

    all_types = adata.obs['popv_prediction'].astype('category').cat.categories
    cmap = plt.colormaps['turbo']
    colors = {ct: cmap(i / max(len(all_types) - 1, 1))
              for i, ct in enumerate(all_types)}

    for ci, cluster in enumerate(clusters):
        sub = adata.obs[adata.obs['clusters'] == cluster]
        cts = sub['popv_prediction'].value_counts(normalize=True) * 100
        bottom = 0
        for ct in all_types:
            pct = cts.get(ct, 0)
            if pct > 0:
                ax.bar(ci, pct, bottom=bottom, color=colors[ct],
                       width=0.7, alpha=1.0 if pct >= 5 else 0.4,
                       edgecolor='white', linewidth=0.3)
                bottom += pct

    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels([f'C{c}' for c in clusters], rotation=45, ha='right')
    ax.set_ylabel('Percentage')
    ax.set_title('Cell Type Composition per Cluster')
    ax.set_ylim(0, 100)
    plt.tight_layout()
    save_fig(fig, 'Stacked_Bar_Cluster_Composition', fig_dir)


def plot_spatial_scatter(adata_puck, puck_name, color_by, fig_dir):
    """Spatial scatter plot for a single puck colored by annotation."""
    categories = sorted(adata_puck.obs[color_by].unique())
    cmap = plt.colormaps['turbo']
    colors = {cat: cmap(i / max(len(categories) - 1, 1))
              for i, cat in enumerate(categories)}

    fig, ax = plt.subplots(figsize=(12, 10))
    for cat in categories:
        m = adata_puck.obs[color_by] == cat
        ax.scatter(
            adata_puck.obs.loc[m, 'x_coord'],
            adata_puck.obs.loc[m, 'y_coord'],
            c=[colors[cat]], s=1, label=cat, alpha=0.7, rasterized=True,
        )
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Y (um)')
    ax.set_title(f'{puck_name}: {color_by}')
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
              fontsize=FONT_SIZE - 14, markerscale=5)
    plt.tight_layout()
    save_fig(fig, f'{puck_name}_spatial_{color_by}', fig_dir)


# =========================================================================
# EXPORT SUMMARIES
# =========================================================================

def export_summaries(adata, ann_dir):
    """Export annotation and QC summary tables."""
    banner("EXPORTING SUMMARIES")

    # Annotation summary: cell types x pucks
    if 'final_annotation' in adata.obs.columns:
        ct_summary = pd.crosstab(
            adata.obs['puck_id'], adata.obs['final_annotation']
        )
        ct_summary['total_cells'] = ct_summary.sum(axis=1)
        ct_path = os.path.join(ann_dir, "annotation_summary.tsv")
        ct_summary.to_csv(ct_path, sep='\t')
        log(f"  Saved: {ct_path}")

    # Cluster mapping: which cluster maps to which cell type
    if 'final_annotation' in adata.obs.columns:
        cluster_map = adata.obs.groupby('clusters')['final_annotation'].first()
        cluster_sizes = adata.obs['clusters'].value_counts().sort_index()
        map_df = pd.DataFrame({
            'cluster': cluster_map.index,
            'cell_type': cluster_map.values,
            'n_cells': [cluster_sizes.get(c, 0) for c in cluster_map.index],
        })
        map_path = os.path.join(ann_dir, "cluster_to_celltype_map.tsv")
        map_df.to_csv(map_path, sep='\t', index=False)
        log(f"  Saved: {map_path}")


# =========================================================================
# MAIN
# =========================================================================

def main():
    banner("STEP 04: CELL TYPE ANNOTATION AND CLUSTERING")
    log(f"  Input:  {DIR_03_QC}")
    log(f"  Output: {DIR_04_ANNOTATION}")

    ann_dir = ensure_dir(DIR_04_ANNOTATION)
    ann_fig_dir = ensure_dir(os.path.join(DIR_04_ANNOTATION, "figures"))
    cache_dir = os.path.join(ann_dir, "popv_cache")

    # --- Load merged QC data ---
    banner("LOADING MERGED QC DATA")
    merged_path = os.path.join(DIR_03_QC, "all_pucks_merged_QC.h5ad")
    if not os.path.exists(merged_path):
        log(f"  FATAL: {merged_path} not found. Run Step03 first.")
        sys.exit(1)

    adata = sc.read_h5ad(merged_path)
    log(f"  Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    log(f"  Pucks: {adata.obs['puck_id'].value_counts().to_dict()}")

    # --- popV annotation (on raw counts, ClusterCatcher flow) ---
    adata = run_popv_annotation(adata, cache_dir)

    # Save popV-annotated (pre-clustering, checkpoint)
    popv_path = os.path.join(ann_dir, "all_pucks_popv_annotated.h5ad")
    adata.write_h5ad(popv_path)
    log(f"  Checkpoint saved: {popv_path}")

    # --- Post-annotation processing ---
    adata = post_annotation_processing(adata)

    # --- Cluster-based consensus annotation ---
    adata = assign_cluster_based_annotation(adata)

    # --- Visualization ---
    generate_annotation_plots(adata, ann_fig_dir)

    # --- Export summaries ---
    export_summaries(adata, ann_dir)

    # --- Save final annotated data ---
    banner("SAVING FINAL ANNOTATED DATA")

    final_path = os.path.join(ann_dir, "all_pucks_annotated.h5ad")
    adata.write_h5ad(final_path)
    size_mb = os.path.getsize(final_path) / 1e6
    log(f"  Saved: {final_path} ({size_mb:.1f} MB)")

    # Per-puck annotated h5ad
    for puck_name in PUCK_NAMES:
        mask = adata.obs['puck_id'] == puck_name
        if mask.sum() == 0:
            continue
        puck_path = os.path.join(ann_dir, f"{puck_name}_annotated.h5ad")
        adata[mask].copy().write_h5ad(puck_path)
        log(f"  Saved: {puck_path}")

    # --- Final summary ---
    banner("STEP 04 SUMMARY")
    log(f"  Total cells:   {adata.n_obs:,}")
    log(f"  Total genes:   {adata.n_vars:,}")
    log(f"  Clusters:      {adata.obs['clusters'].nunique()}")
    if 'final_annotation' in adata.obs.columns:
        log(f"  Cell types:    {adata.obs['final_annotation'].nunique()}")
    log(f"  Layers:        {list(adata.layers.keys())}")

    log(f"\n  Cells per puck:")
    for puck_name in PUCK_NAMES:
        n = (adata.obs['puck_id'] == puck_name).sum()
        log(f"    {puck_name}: {n:,}")

    banner("STEP 04 COMPLETE")
    log(f"  Output: {DIR_04_ANNOTATION}/")
    log(f"  Next: Step06 (SComatic mutation calling) or")
    log(f"        Step05 (HPV detection with Kraken2)")


if __name__ == "__main__":
    main()
