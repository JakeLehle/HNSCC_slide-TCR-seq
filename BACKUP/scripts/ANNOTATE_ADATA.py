#!/usr/bin/env python3
"""
QC and Cell Type Annotation Pipeline for Slide-TCR-seq Data

This script performs:
1. Load all three puck AnnData objects
2. QC filtering (outliers, mitochondrial, etc.)
3. Automated cell type annotation with popV
4. Leiden clustering and UMAP visualization
5. Spatial visualization with Squidpy

Follows the initial data loading from load_slideseq_data.py

Author: Jake Lehle
Project: Slide-TCR-seq HPV+ HNSCC analysis
Date: December 2025
"""

#%% Imports
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sb
from scipy.stats import median_abs_deviation
from pathlib import Path
from typing import Dict, List, Optional
import warnings
import anndata as ad
import popv

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

#%% Configuration
# Update these paths for your system
ANNDATA_DIR = "/work/sdz852/WORKING/slide-TCR-seq/analysis/anndata"
OUTPUT_DIR = "/work/sdz852/WORKING/slide-TCR-seq/analysis/processed"
FIGURES_DIR = "/work/sdz852/WORKING/slide-TCR-seq/analysis/figures"

# Create output directories
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)

# Set scanpy figure directory
sc.settings.figdir = FIGURES_DIR
sc.settings.verbosity = 2

# Puck information
PUCK_INFO = {
    'Puck_211214_29': {'sample_barcode': 'AGATTTAA', 'description': 'HPV+ HNSCC sample 29'},
    'Puck_211214_37': {'sample_barcode': 'GGCGTCGA', 'description': 'HPV+ HNSCC sample 37'},
    'Puck_211214_40': {'sample_barcode': 'ATCACTCG', 'description': 'HPV+ HNSCC sample 40 (best RNA)'}
}

#%% Helper Functions

def is_outlier(adata: ad.AnnData, metric: str, nmads: int) -> pd.Series:
    """
    Identify outliers based on median absolute deviation (MAD).
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with QC metrics in .obs
    metric : str
        Column name in adata.obs to check for outliers
    nmads : int
        Number of MADs to use as threshold
        
    Returns
    -------
    outlier : pd.Series
        Boolean series indicating outlier cells
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def QC_spatial_adata(adata: ad.AnnData, puck_id: str, min_genes: int = 200, 
                     min_cells: int = 20, mt_threshold: float = 20.0) -> ad.AnnData:
    """
    Perform QC filtering on a spatial AnnData object.
    
    Adapted for Slide-seq spatial data with considerations for:
    - Lower UMI counts typical of spatial technologies
    - Preservation of spatial structure during filtering
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    puck_id : str
        Puck identifier for labeling plots
    min_genes : int
        Minimum genes per cell (default: 200, lower than standard scRNA-seq)
    min_cells : int
        Minimum cells per gene
    mt_threshold : float
        Maximum mitochondrial percentage
        
    Returns
    -------
    adata : AnnData
        Filtered AnnData object
    """
    print(f"\n{'='*60}")
    print(f"QC FILTERING: {puck_id}")
    print(f"{'='*60}")
    print(f"Starting cells: {adata.n_obs:,}")
    print(f"Starting genes: {adata.n_vars:,}")
    
    # Add gene annotations
    adata.var["mt"] = adata.var_names.str.startswith(('mt-', 'MT-', 'Mt-'))
    adata.var["ribo"] = adata.var_names.str.startswith(('Rps', 'Rpl', 'RPS', 'RPL'))
    adata.var["hb"] = adata.var_names.str.contains('^HB[^(P)]', regex=True)
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], 
        inplace=True, percent_top=[20], log1p=True
    )
    
    # Generate QC plots before filtering
    print(f"\nGenerating QC plots (before filtering)...")
    
    # Scatter plot: counts vs genes colored by MT%
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", 
                  color="pct_counts_mt", ax=axes[0], show=False, title=f"{puck_id}: Counts vs Genes")
    
    # Histogram of total counts
    axes[1].hist(adata.obs["total_counts"], bins=100, edgecolor='black', alpha=0.7)
    axes[1].set_xlabel("Total counts")
    axes[1].set_ylabel("Number of cells")
    axes[1].set_title(f"{puck_id}: UMI Distribution")
    axes[1].axvline(x=np.median(adata.obs["total_counts"]), color='red', linestyle='--', 
                    label=f'Median: {np.median(adata.obs["total_counts"]):.0f}')
    axes[1].legend()
    
    # Histogram of MT%
    axes[2].hist(adata.obs["pct_counts_mt"], bins=50, edgecolor='black', alpha=0.7)
    axes[2].set_xlabel("Mitochondrial %")
    axes[2].set_ylabel("Number of cells")
    axes[2].set_title(f"{puck_id}: MT% Distribution")
    axes[2].axvline(x=mt_threshold, color='red', linestyle='--', label=f'Threshold: {mt_threshold}%')
    axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f'{puck_id}_QC_before_filtering.pdf'), bbox_inches='tight')
    plt.close()
    
    # Spatial QC plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for ax, color, title in zip(axes, 
                                 ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'],
                                 ['Total UMIs', 'Genes Detected', 'MT %']):
        scatter = ax.scatter(adata.obs['x_coord'], adata.obs['y_coord'], 
                            c=adata.obs[color], s=1, cmap='viridis', alpha=0.7)
        ax.set_title(f"{puck_id}: {title}")
        ax.set_xlabel("X coordinate (μm)")
        ax.set_ylabel("Y coordinate (μm)")
        ax.set_aspect('equal')
        plt.colorbar(scatter, ax=ax)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f'{puck_id}_Spatial_QC_before_filtering.pdf'), bbox_inches='tight')
    plt.close()
    
    # Identify outliers using MAD method
    # More lenient for spatial data (using 5 MADs instead of 3)
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    
    # MT outliers (stricter, using 3 MADs + hard threshold)
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > mt_threshold
    )
    
    # Report outlier counts
    print(f"\nOutlier detection:")
    print(f"  General outliers: {adata.obs['outlier'].sum():,} ({100*adata.obs['outlier'].mean():.1f}%)")
    print(f"  MT outliers: {adata.obs['mt_outlier'].sum():,} ({100*adata.obs['mt_outlier'].mean():.1f}%)")
    
    # Filter cells
    n_before = adata.n_obs
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"\nAfter outlier filtering: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")
    
    # Filter by minimum genes
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print(f"After min_genes={min_genes} filter: {adata.n_obs:,} cells ({n_before - adata.n_obs:,} removed)")
    
    # Filter genes
    n_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print(f"After min_cells={min_cells} filter: {adata.n_vars:,} genes ({n_before - adata.n_vars:,} removed)")
    
    # Generate QC plots after filtering
    print(f"\nGenerating QC plots (after filtering)...")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", 
                  color="pct_counts_mt", ax=axes[0], show=False, title=f"{puck_id}: After QC")
    
    axes[1].hist(adata.obs["total_counts"], bins=100, edgecolor='black', alpha=0.7)
    axes[1].set_xlabel("Total counts")
    axes[1].set_ylabel("Number of cells")
    axes[1].set_title(f"{puck_id}: UMI Distribution (filtered)")
    
    axes[2].hist(adata.obs["pct_counts_mt"], bins=50, edgecolor='black', alpha=0.7)
    axes[2].set_xlabel("Mitochondrial %")
    axes[2].set_ylabel("Number of cells")
    axes[2].set_title(f"{puck_id}: MT% Distribution (filtered)")
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f'{puck_id}_QC_after_filtering.pdf'), bbox_inches='tight')
    plt.close()
    
    # Final summary
    print(f"\nFinal: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    print(f"Median UMIs/cell: {adata.obs['total_counts'].median():.0f}")
    print(f"Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"Median MT%: {adata.obs['pct_counts_mt'].median():.1f}%")
    
    return adata


def load_all_pucks(anndata_dir: str) -> Dict[str, ad.AnnData]:
    """
    Load all puck AnnData objects from h5ad files.
    
    Parameters
    ----------
    anndata_dir : str
        Directory containing *_anndata.h5ad files
        
    Returns
    -------
    adatas : dict
        Dictionary of AnnData objects keyed by puck ID
    """
    print("\n" + "=" * 60)
    print("LOADING PUCK DATA")
    print("=" * 60)
    
    adatas = {}
    anndata_dir = Path(anndata_dir)
    
    for puck_id in PUCK_INFO.keys():
        h5ad_file = anndata_dir / f"{puck_id}_anndata.h5ad"
        
        if h5ad_file.exists():
            print(f"\nLoading {puck_id}...")
            adata = sc.read_h5ad(h5ad_file)
            adata.obs['puck_id'] = puck_id
            adatas[puck_id] = adata
            print(f"  ✓ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        else:
            print(f"  ✗ File not found: {h5ad_file}")
    
    return adatas


def merge_pucks(adatas: Dict[str, ad.AnnData]) -> ad.AnnData:
    """
    Merge multiple puck AnnData objects into a single object.
    
    Parameters
    ----------
    adatas : dict
        Dictionary of AnnData objects
        
    Returns
    -------
    adata_merged : AnnData
        Merged AnnData object with puck_id in .obs
    """
    print("\n" + "=" * 60)
    print("MERGING PUCKS")
    print("=" * 60)
    
    # Concatenate all pucks
    adata_list = list(adatas.values())
    puck_ids = list(adatas.keys())
    
    # Use 'inner' join to only keep genes present in ALL pucks
    # This avoids NaN values that cause issues downstream
    adata_merged = ad.concat(
        adata_list,
        join='inner',  # Changed from 'outer' to avoid NaN
        label='puck_id',
        keys=puck_ids,
        index_unique='_'
    )
    
    # Make var_names unique
    adata_merged.var_names_make_unique()
    
    # Verify no NaN values
    if hasattr(adata_merged.X, 'toarray'):
        n_nan = np.isnan(adata_merged.X.toarray()).sum()
    else:
        n_nan = np.isnan(adata_merged.X).sum()
    
    if n_nan > 0:
        print(f"  WARNING: {n_nan:,} NaN values detected after merge")
    
    print(f"\nMerged AnnData:")
    print(f"  Total cells: {adata_merged.n_obs:,}")
    print(f"  Total genes: {adata_merged.n_vars:,}")
    print(f"  Cells per puck:")
    for puck_id in puck_ids:
        n_cells = (adata_merged.obs['puck_id'] == puck_id).sum()
        print(f"    {puck_id}: {n_cells:,}")
    
    return adata_merged


def annotate_with_popv(adata: ad.AnnData, 
                       huggingface_repo: str = "popV/tabula_sapiens_All_Cells",
                       batch_key: str = "puck_id") -> ad.AnnData:
    """
    Perform automated cell type annotation using popV.
    
    Parameters
    ----------
    adata : AnnData
        QC-filtered AnnData object (raw counts, not normalized)
    huggingface_repo : str
        HuggingFace model repository for popV
    batch_key : str
        Batch key for integration
        
    Returns
    -------
    adata : AnnData
        AnnData with popV annotations in .obs
    """
    print("\n" + "=" * 60)
    print("CELL TYPE ANNOTATION WITH popV")
    print("=" * 60)
    
    import popv
    
    # Pull the pre-trained model
    print(f"\nLoading model from: {huggingface_repo}")
    cache_dir = os.path.join(OUTPUT_DIR, "popv_cache")
    hmo = popv.hub.HubModel.pull_from_huggingface_hub(huggingface_repo, cache_dir=cache_dir)
    
    # Annotate
    print(f"\nAnnotating {adata.n_obs:,} cells...")
    adata_annotated = hmo.annotate_data(
        adata,
        query_batch_key=batch_key,
        prediction_mode="inference",
        gene_symbols='feature_name' if 'feature_name' in adata.var.columns else None
    )
    
    # Filter to query cells only (remove reference cells)
    adata_annotated = adata_annotated[adata_annotated.obs["_dataset"] == "query"].copy()
    
    print(f"\nAnnotation complete!")
    print(f"Unique cell types: {adata_annotated.obs['popv_prediction'].nunique()}")
    print(f"\nTop 10 cell types:")
    print(adata_annotated.obs['popv_prediction'].value_counts().head(10))
    
    return adata_annotated


def cluster_and_embed(adata: ad.AnnData, resolution: float = 1.0) -> ad.AnnData:
    """
    Normalize, cluster, and compute UMAP embedding.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object
    resolution : float
        Leiden clustering resolution
        
    Returns
    -------
    adata : AnnData
        Processed AnnData with clusters and UMAP
    """
    print("\n" + "=" * 60)
    print("CLUSTERING AND EMBEDDING")
    print("=" * 60)
    
    # Make a copy for processing
    adata_pp = adata.copy()
    
    # Handle NaN values from outer join during merge
    # Replace NaN with 0 (gene not detected in that puck)
    if hasattr(adata_pp.X, 'toarray'):
        # Sparse matrix
        X_dense = adata_pp.X.toarray()
    else:
        X_dense = adata_pp.X.copy()
    
    n_nan = np.isnan(X_dense).sum()
    if n_nan > 0:
        print(f"\nReplacing {n_nan:,} NaN values with 0...")
        X_dense = np.nan_to_num(X_dense, nan=0.0)
        adata_pp.X = X_dense
    
    # Also handle any infinite values
    n_inf = np.isinf(X_dense).sum()
    if n_inf > 0:
        print(f"Replacing {n_inf:,} infinite values with 0...")
        adata_pp.X = np.nan_to_num(adata_pp.X, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Filter genes with zero variance (can cause issues)
    print("\nFiltering genes with zero counts across all cells...")
    gene_counts = np.asarray(adata_pp.X.sum(axis=0)).flatten()
    nonzero_genes = gene_counts > 0
    n_removed = (~nonzero_genes).sum()
    if n_removed > 0:
        print(f"  Removing {n_removed:,} genes with zero counts")
        adata_pp = adata_pp[:, nonzero_genes].copy()
    
    print(f"  Remaining: {adata_pp.n_vars:,} genes")
    
    # Normalize
    print("\nNormalizing...")
    sc.pp.normalize_total(adata_pp, target_sum=1e6)
    sc.pp.log1p(adata_pp)
    
    # Store raw counts
    adata_pp.raw = adata_pp
    
    # Feature selection - use flavor='seurat_v3' which is more robust
    print("Selecting highly variable genes...")
    try:
        sc.pp.highly_variable_genes(
            adata_pp, 
            min_mean=0.0125, 
            max_mean=3, 
            min_disp=0.5,
            batch_key='puck_id'  # Account for batch effects
        )
    except Exception as e:
        print(f"  Batch-aware HVG selection failed: {e}")
        print("  Trying without batch correction...")
        sc.pp.highly_variable_genes(
            adata_pp, 
            flavor='seurat_v3',
            n_top_genes=2000
        )
    
    n_hvg = adata_pp.var['highly_variable'].sum()
    print(f"  HVGs: {n_hvg:,}")
    
    if n_hvg < 100:
        print("  WARNING: Very few HVGs detected. Using top 2000 variable genes instead.")
        sc.pp.highly_variable_genes(adata_pp, flavor='seurat_v3', n_top_genes=2000)
        print(f"  HVGs (seurat_v3): {adata_pp.var['highly_variable'].sum():,}")
    
    # Subset to HVGs for PCA
    adata_hvg = adata_pp[:, adata_pp.var['highly_variable']].copy()
    
    # PCA
    print("Computing PCA...")
    n_comps = min(50, adata_hvg.n_vars - 1, adata_hvg.n_obs - 1)
    sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=n_comps)
    
    # Copy PCA results back to full object
    adata_pp.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata_pp.uns['pca'] = adata_hvg.uns['pca']
    
    # Neighbors
    print("Computing neighbors...")
    n_pcs = min(30, n_comps)
    sc.pp.neighbors(adata_pp, n_pcs=n_pcs, use_rep='X_pca')
    
    # UMAP
    print("Computing UMAP...")
    sc.tl.umap(adata_pp)
    
    # Leiden clustering
    print(f"Clustering (resolution={resolution})...")
    sc.tl.leiden(adata_pp, key_added='clusters', resolution=resolution, random_state=42)
    print(f"  Found {adata_pp.obs['clusters'].nunique()} clusters")
    
    return adata_pp


def assign_cluster_based_annotation(adata: ad.AnnData) -> ad.AnnData:
    """
    Assign final cell type annotation based on cluster-level consensus.
    
    Uses weighted voting within each cluster to assign the dominant
    cell type based on popV prediction scores.
    
    Parameters
    ----------
    adata : AnnData
        AnnData with 'clusters' and 'popv_prediction' columns
        
    Returns
    -------
    adata : AnnData
        AnnData with 'final_annotation' column added
    """
    print("\n" + "=" * 60)
    print("ASSIGNING CLUSTER-BASED ANNOTATIONS")
    print("=" * 60)
    
    # Make a copy of the relevant columns
    df = adata.obs[['clusters', 'popv_prediction', 'popv_prediction_score']].copy()
    
    # Ensure scores are numeric
    df['popv_prediction_score'] = pd.to_numeric(df['popv_prediction_score'], errors='coerce')
    
    # Min-Max normalization within each cluster
    def min_max_normalize(x):
        x_min = x.min()
        x_range = x.max() - x_min
        if x_range == 0:
            return np.zeros(len(x))
        return (x - x_min) / x_range
    
    df['normalized_score'] = df.groupby('clusters')['popv_prediction_score'].transform(min_max_normalize)
    
    # Linear weighting
    def linear_weights(x):
        total = x.sum()
        if total == 0:
            return np.ones(len(x)) / len(x)
        return x / total
    
    df['weight'] = df.groupby('clusters')['normalized_score'].transform(linear_weights)
    
    # Calculate weighted score
    df['weighted_score'] = df['popv_prediction_score'] * df['weight']
    
    # Aggregate scores by cluster and predicted cell type
    cluster_type_scores = df.groupby(['clusters', 'popv_prediction'])['weighted_score'].sum().reset_index()
    
    # For each cluster, get the cell type with highest aggregated weighted score
    dominant_types = cluster_type_scores.loc[
        cluster_type_scores.groupby('clusters')['weighted_score'].idxmax()
    ].set_index('clusters')['popv_prediction']
    
    # Map the dominant type back to all cells
    adata.obs['final_annotation'] = adata.obs['clusters'].map(dominant_types)
    
    print(f"\nFinal cell type assignments:")
    print(adata.obs['final_annotation'].value_counts())
    
    return adata


def plot_spatial_annotations(adata: ad.AnnData, puck_id: str, color_by: str = 'final_annotation'):
    """
    Create spatial scatter plot colored by cell type annotation.
    
    Parameters
    ----------
    adata : AnnData
        AnnData with spatial coordinates and annotations
    puck_id : str
        Puck identifier for title and filename
    color_by : str
        Column in .obs to color by
    """
    # Get unique categories and create color palette
    categories = adata.obs[color_by].astype('category').cat.categories
    cmap = plt.get_cmap('turbo')
    values = np.linspace(0, 1, len(categories))
    colors = {cat: cmap(val) for cat, val in zip(categories, values)}
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    for cat in categories:
        mask = adata.obs[color_by] == cat
        ax.scatter(
            adata.obs.loc[mask, 'x_coord'],
            adata.obs.loc[mask, 'y_coord'],
            c=[colors[cat]],
            s=2,
            label=cat,
            alpha=0.7
        )
    
    ax.set_xlabel('X coordinate (μm)', fontsize=12)
    ax.set_ylabel('Y coordinate (μm)', fontsize=12)
    ax.set_title(f'{puck_id}: {color_by}', fontsize=14)
    ax.set_aspect('equal')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, markerscale=3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f'{puck_id}_Spatial_{color_by}.pdf'), 
                bbox_inches='tight', dpi=150)
    plt.close()


def plot_umap_annotations(adata: ad.AnnData, color_by: str = 'final_annotation', 
                          filename_suffix: str = ''):
    """
    Create UMAP plot colored by annotation.
    """
    import matplotlib.patheffects as pe
    from textwrap import wrap
    
    categories = adata.obs[color_by].astype('category').cat.categories
    cmap = plt.get_cmap('turbo')
    values = np.linspace(0, 1, len(categories))
    palette = [cmap(val) for val in values]
    
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color=color_by, palette=palette, ax=ax, show=False, 
               legend_loc='right margin', frameon=False, size=10)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f'UMAP_{color_by}{filename_suffix}.pdf'), 
                bbox_inches='tight', dpi=150)
    plt.close()


#%% Main Pipeline

def main():
    """
    Main pipeline execution.
    """
    print("\n" + "=" * 70)
    print("SLIDE-TCR-SEQ QC AND ANNOTATION PIPELINE")
    print("HPV+ HNSCC - Sophia Liu (Ragon Institute)")
    print("=" * 70)
    
    # Step 1: Load all pucks
    adatas = load_all_pucks(ANNDATA_DIR)
    
    if not adatas:
        print("\n✗ ERROR: No pucks loaded!")
        sys.exit(1)
    
    # Step 2: QC each puck individually
    print("\n" + "=" * 70)
    print("STEP 2: QUALITY CONTROL")
    print("=" * 70)
    
    adatas_qc = {}
    for puck_id, adata in adatas.items():
        # Use more lenient parameters for spatial data
        adata_qc = QC_spatial_adata(
            adata, 
            puck_id, 
            min_genes=100,      # Lower threshold for spatial
            min_cells=10,       # Lower threshold for spatial
            mt_threshold=25.0   # Slightly higher MT threshold
        )
        adatas_qc[puck_id] = adata_qc
        
        # Save QC'd individual puck
        output_file = os.path.join(OUTPUT_DIR, f"{puck_id}_QC.h5ad")
        adata_qc.write_h5ad(output_file)
        print(f"  Saved: {output_file}")
    
    # Step 3: Merge pucks
    adata_merged = merge_pucks(adatas_qc)
    
    # Save merged (pre-annotation)
    merged_file = os.path.join(OUTPUT_DIR, "all_pucks_merged_QC.h5ad")
    adata_merged.write_h5ad(merged_file)
    print(f"\nSaved merged QC data: {merged_file}")
    
    # Step 4: Cell type annotation with popV
    print("\n" + "=" * 70)
    print("STEP 4: CELL TYPE ANNOTATION")
    print("=" * 70)
    
    try:
        adata_annotated = annotate_with_popv(adata_merged, batch_key='puck_id')
        
        # Save annotated
        annotated_file = os.path.join(OUTPUT_DIR, "all_pucks_popv_annotated.h5ad")
        adata_annotated.write_h5ad(annotated_file)
        print(f"\nSaved popV annotated data: {annotated_file}")
        
    except ImportError:
        print("\n⚠ WARNING: popV not installed. Skipping automated annotation.")
        print("  Install with: pip install popv")
        print("  Continuing with clustering only...")
        adata_annotated = adata_merged
    except Exception as e:
        print(f"\n⚠ WARNING: popV annotation failed: {e}")
        print("  Continuing with clustering only...")
        adata_annotated = adata_merged
    
    # Step 5: Clustering and embedding
    adata_pp = cluster_and_embed(adata_annotated, resolution=1.0)
    
    # Step 6: Assign cluster-based annotations (if popV was successful)
    if 'popv_prediction' in adata_pp.obs.columns:
        adata_pp = assign_cluster_based_annotation(adata_pp)
        annotation_col = 'final_annotation'
    else:
        annotation_col = 'clusters'
    
    # Step 7: Visualization
    print("\n" + "=" * 70)
    print("STEP 7: VISUALIZATION")
    print("=" * 70)
    
    # UMAP plots
    plot_umap_annotations(adata_pp, 'clusters', '_clusters')
    plot_umap_annotations(adata_pp, 'puck_id', '_by_puck')
    
    if 'popv_prediction' in adata_pp.obs.columns:
        plot_umap_annotations(adata_pp, 'popv_prediction', '_popv')
        plot_umap_annotations(adata_pp, 'final_annotation', '_final')
    
    # Spatial plots for each puck
    for puck_id in PUCK_INFO.keys():
        puck_mask = adata_pp.obs['puck_id'] == puck_id
        if puck_mask.sum() > 0:
            adata_puck = adata_pp[puck_mask].copy()
            plot_spatial_annotations(adata_puck, puck_id, annotation_col)
            plot_spatial_annotations(adata_puck, puck_id, 'clusters')
    
    # Step 8: Save final processed data
    print("\n" + "=" * 70)
    print("STEP 8: SAVING FINAL DATA")
    print("=" * 70)
    
    # Save combined processed
    final_file = os.path.join(OUTPUT_DIR, "all_pucks_processed.h5ad")
    adata_pp.write_h5ad(final_file)
    print(f"Saved final processed data: {final_file}")
    
    # Save individual processed pucks
    for puck_id in PUCK_INFO.keys():
        puck_mask = adata_pp.obs['puck_id'] == puck_id
        if puck_mask.sum() > 0:
            adata_puck = adata_pp[puck_mask].copy()
            puck_file = os.path.join(OUTPUT_DIR, f"{puck_id}_processed.h5ad")
            adata_puck.write_h5ad(puck_file)
            print(f"Saved: {puck_file}")
    
    # Summary
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE!")
    print("=" * 70)
    
    print(f"\nFinal dataset summary:")
    print(f"  Total cells: {adata_pp.n_obs:,}")
    print(f"  Total genes: {adata_pp.n_vars:,}")
    print(f"  Clusters: {adata_pp.obs['clusters'].nunique()}")
    if 'final_annotation' in adata_pp.obs.columns:
        print(f"  Cell types: {adata_pp.obs['final_annotation'].nunique()}")
    
    print(f"\nCells per puck:")
    for puck_id in PUCK_INFO.keys():
        n_cells = (adata_pp.obs['puck_id'] == puck_id).sum()
        print(f"  {puck_id}: {n_cells:,}")
    
    print(f"\nOutput files saved to: {OUTPUT_DIR}")
    print(f"Figures saved to: {FIGURES_DIR}")
    
    print("\nNext steps:")
    print("  1. Review QC plots and adjust thresholds if needed")
    print("  2. Examine spatial distribution of cell types")
    print("  3. Align matched.bam to HPV16 to identify HPV+ regions")
    print("  4. Compare cell type composition between HPV+ and HPV- regions")
    print("  5. Run spatial neighborhood analysis with Squidpy")
    
    return adata_pp


#%% Run pipeline
if __name__ == "__main__":
    adata_final = main()
