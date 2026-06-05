#!/usr/bin/env python3
"""
Step03_QC_Filter.py

Quality control filtering for Slide-seq spatial data.

Per-puck QC:
  - MAD-based outlier detection (5 MADs general, 3 MADs for MT%)
  - Hard mitochondrial threshold (25%)
  - Minimum genes per cell (100) and minimum cells per gene (10)
  - Spatial QC plots (before and after filtering)

Then merge all pucks (inner join on genes) and convert to sparse format.

Reads from:  data/outputs/02_anndata/
Writes to:   data/outputs/03_qc/

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
from scipy.stats import median_abs_deviation
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spatial_config import (
    DIR_02_ANNDATA, DIR_03_QC, DIR_FIGURES, DIR_TROUBLESHOOTING,
    PUCK_NAMES, MIN_GENES, MIN_CELLS, MT_THRESHOLD, N_HVG,
    FONT_SIZE, DPI, COLOR_PUCK_29, COLOR_PUCK_37, COLOR_PUCK_40,
    banner, log, ensure_dir, save_fig,
)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# =========================================================================
# QC FUNCTIONS
# =========================================================================

def is_outlier(adata, metric, nmads):
    """Identify outliers using median absolute deviation."""
    M = adata.obs[metric]
    med = np.median(M)
    mad = median_abs_deviation(M)
    return (M < med - nmads * mad) | (M > med + nmads * mad)


def qc_single_puck(adata, puck_name, qc_fig_dir):
    """
    Run QC filtering on a single puck. Returns filtered AnnData.

    Steps:
      1. Add gene annotations (MT, ribo, hemoglobin)
      2. Recompute QC metrics with full annotation set
      3. Generate pre-filtering spatial QC plots
      4. MAD-based outlier detection + hard MT threshold
      5. min_genes / min_cells filtering
      6. Generate post-filtering plots
      7. Convert to sparse CSR matrix
    """
    banner(f"QC: {puck_name}", char="-")
    log(f"  Starting: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # --- Gene annotations ---
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]', regex=True)

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt', 'ribo', 'hb'],
        inplace=True, percent_top=[20], log1p=True,
    )

    # --- Pre-filtering spatial QC plots ---
    plot_spatial_qc(adata, puck_name, qc_fig_dir, suffix='_before_QC')

    # --- Pre-filtering distribution plots ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    axes[0].hist(adata.obs['total_counts'], bins=100, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Total UMIs')
    axes[0].set_ylabel('Cells')
    axes[0].set_title(f'{puck_name}: UMI distribution')
    med_umi = adata.obs['total_counts'].median()
    axes[0].axvline(med_umi, color='red', ls='--', label=f'Median: {med_umi:.0f}')
    axes[0].legend(fontsize=FONT_SIZE - 10)

    axes[1].hist(adata.obs['n_genes_by_counts'], bins=100, edgecolor='black', alpha=0.7)
    axes[1].set_xlabel('Genes detected')
    axes[1].set_ylabel('Cells')
    axes[1].set_title(f'{puck_name}: Gene distribution')

    axes[2].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7)
    axes[2].set_xlabel('Mitochondrial %')
    axes[2].set_ylabel('Cells')
    axes[2].set_title(f'{puck_name}: MT% distribution')
    axes[2].axvline(MT_THRESHOLD, color='red', ls='--',
                    label=f'Threshold: {MT_THRESHOLD}%')
    axes[2].legend(fontsize=FONT_SIZE - 10)

    plt.tight_layout()
    save_fig(fig, f'{puck_name}_distributions_before_QC', qc_fig_dir)

    # --- Outlier detection ---
    adata.obs['outlier'] = (
        is_outlier(adata, 'log1p_total_counts', 5)
        | is_outlier(adata, 'log1p_n_genes_by_counts', 5)
        | is_outlier(adata, 'pct_counts_in_top_20_genes', 5)
    )
    adata.obs['mt_outlier'] = (
        is_outlier(adata, 'pct_counts_mt', 3)
        | (adata.obs['pct_counts_mt'] > MT_THRESHOLD)
    )

    n_general = adata.obs['outlier'].sum()
    n_mt = adata.obs['mt_outlier'].sum()
    log(f"  Outliers detected:")
    log(f"    General (5 MAD): {n_general:,} ({100 * n_general / adata.n_obs:.1f}%)")
    log(f"    MT (3 MAD + {MT_THRESHOLD}%): {n_mt:,} ({100 * n_mt / adata.n_obs:.1f}%)")

    # --- Apply filters ---
    n0 = adata.n_obs
    adata = adata[(~adata.obs['outlier']) & (~adata.obs['mt_outlier'])].copy()
    log(f"  After outlier removal: {adata.n_obs:,} cells ({n0 - adata.n_obs:,} removed)")

    n1 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    log(f"  After min_genes={MIN_GENES}: {adata.n_obs:,} cells ({n1 - adata.n_obs:,} removed)")

    n2 = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    log(f"  After min_cells={MIN_CELLS}: {adata.n_vars:,} genes ({n2 - adata.n_vars:,} removed)")

    # --- Post-filtering plots ---
    plot_spatial_qc(adata, puck_name, qc_fig_dir, suffix='_after_QC')

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    axes[0].hist(adata.obs['total_counts'], bins=100, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Total UMIs')
    axes[0].set_ylabel('Cells')
    axes[0].set_title(f'{puck_name}: UMI (after QC)')

    axes[1].hist(adata.obs['n_genes_by_counts'], bins=100, edgecolor='black', alpha=0.7)
    axes[1].set_xlabel('Genes detected')
    axes[1].set_ylabel('Cells')
    axes[1].set_title(f'{puck_name}: Genes (after QC)')

    axes[2].hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7)
    axes[2].set_xlabel('Mitochondrial %')
    axes[2].set_ylabel('Cells')
    axes[2].set_title(f'{puck_name}: MT% (after QC)')

    plt.tight_layout()
    save_fig(fig, f'{puck_name}_distributions_after_QC', qc_fig_dir)

    # --- Convert to sparse ---
    if not sparse.issparse(adata.X):
        log(f"  Converting to sparse CSR matrix...")
        dense_bytes = adata.X.nbytes
        adata.X = sparse.csr_matrix(adata.X)
        sparse_bytes = adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes
        log(f"    Dense: {dense_bytes / 1e9:.1f} GB -> Sparse: {sparse_bytes / 1e6:.0f} MB "
            f"({100 * sparse_bytes / dense_bytes:.1f}%)")

    # --- Final summary ---
    log(f"  Final: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    log(f"    Median UMIs/cell: {adata.obs['total_counts'].median():.0f}")
    log(f"    Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    log(f"    Median MT%: {adata.obs['pct_counts_mt'].median():.1f}%")

    return adata


def plot_spatial_qc(adata, puck_name, fig_dir, suffix=''):
    """Generate spatial scatter plots for key QC metrics."""
    metrics = [
        ('total_counts', 'Total UMIs'),
        ('n_genes_by_counts', 'Genes Detected'),
        ('pct_counts_mt', 'MT %'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))

    for ax, (col, title) in zip(axes, metrics):
        sc_obj = ax.scatter(
            adata.obs['x_coord'], adata.obs['y_coord'],
            c=adata.obs[col], s=0.5, cmap='viridis', alpha=0.7, rasterized=True,
        )
        ax.set_title(f'{puck_name}: {title}')
        ax.set_xlabel('X (um)')
        ax.set_ylabel('Y (um)')
        ax.set_aspect('equal')
        plt.colorbar(sc_obj, ax=ax, shrink=0.8)

    plt.tight_layout()
    save_fig(fig, f'{puck_name}_spatial{suffix}', fig_dir)


# =========================================================================
# MERGE FUNCTION
# =========================================================================

def merge_pucks(adatas_qc):
    """Merge all QC'd pucks using inner join on genes."""
    banner("MERGING PUCKS")

    adata_list = list(adatas_qc.values())
    puck_ids = list(adatas_qc.keys())

    # Ensure puck_id is in obs for each
    for pid, adata in adatas_qc.items():
        adata.obs['puck_id'] = pid

    # Inner join: only genes present in ALL pucks (avoids NaN)
    adata_merged = ad.concat(
        adata_list, join='inner',
        label='puck_id', keys=puck_ids, index_unique='_',
    )
    adata_merged.var_names_make_unique()

    log(f"  Merged AnnData: {adata_merged.n_obs:,} cells, {adata_merged.n_vars:,} genes")
    for pid in puck_ids:
        n = (adata_merged.obs['puck_id'] == pid).sum()
        log(f"    {pid}: {n:,} cells")

    # Verify sparse
    if not sparse.issparse(adata_merged.X):
        adata_merged.X = sparse.csr_matrix(adata_merged.X)

    # Verify no NaN
    if sparse.issparse(adata_merged.X):
        n_nan = np.isnan(adata_merged.X.data).sum()
    else:
        n_nan = np.isnan(adata_merged.X).sum()
    if n_nan > 0:
        log(f"  WARNING: {n_nan:,} NaN values in merged matrix")

    return adata_merged


# =========================================================================
# MAIN
# =========================================================================

def main():
    banner("STEP 03: QUALITY CONTROL AND FILTERING")
    log(f"  Input:  {DIR_02_ANNDATA}")
    log(f"  Output: {DIR_03_QC}")

    qc_dir = ensure_dir(DIR_03_QC)
    qc_fig_dir = ensure_dir(os.path.join(DIR_03_QC, "figures"))
    ensure_dir(DIR_TROUBLESHOOTING)

    log(f"  Parameters:")
    log(f"    MIN_GENES:    {MIN_GENES}")
    log(f"    MIN_CELLS:    {MIN_CELLS}")
    log(f"    MT_THRESHOLD: {MT_THRESHOLD}%")

    # --- Load per-puck h5ad ---
    banner("LOADING PER-PUCK H5AD FILES")
    adatas = {}
    for puck_name in PUCK_NAMES:
        h5ad_path = os.path.join(DIR_02_ANNDATA, f"{puck_name}.h5ad")
        if not os.path.exists(h5ad_path):
            log(f"  SKIP: {h5ad_path} not found")
            continue
        log(f"  Loading {puck_name}...")
        adata = sc.read_h5ad(h5ad_path)
        log(f"    {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        adatas[puck_name] = adata

    if not adatas:
        log("FATAL: No pucks loaded.")
        sys.exit(1)

    # --- Per-puck QC ---
    adatas_qc = {}
    qc_rows = []
    for puck_name, adata in adatas.items():
        adata_qc = qc_single_puck(adata, puck_name, qc_fig_dir)
        adatas_qc[puck_name] = adata_qc

        # Save per-puck QC'd h5ad
        out_path = os.path.join(qc_dir, f"{puck_name}_QC.h5ad")
        adata_qc.write_h5ad(out_path)
        size_mb = os.path.getsize(out_path) / 1e6
        log(f"  Saved: {out_path} ({size_mb:.1f} MB)")

        qc_rows.append({
            'puck': puck_name,
            'cells_before': adata.n_obs,
            'cells_after': adata_qc.n_obs,
            'pct_retained': 100 * adata_qc.n_obs / adata.n_obs,
            'genes_after': adata_qc.n_vars,
            'median_umis': adata_qc.obs['total_counts'].median(),
            'median_genes': adata_qc.obs['n_genes_by_counts'].median(),
            'median_mt_pct': adata_qc.obs['pct_counts_mt'].median(),
        })

    # --- QC summary table ---
    qc_df = pd.DataFrame(qc_rows)
    qc_path = os.path.join(qc_dir, "Step03_QC_summary.tsv")
    qc_df.to_csv(qc_path, sep='\t', index=False)
    log(f"  QC summary: {qc_path}")

    # --- Merge pucks ---
    adata_merged = merge_pucks(adatas_qc)

    merged_path = os.path.join(qc_dir, "all_pucks_merged_QC.h5ad")
    adata_merged.write_h5ad(merged_path)
    size_mb = os.path.getsize(merged_path) / 1e6
    log(f"  Saved merged: {merged_path} ({size_mb:.1f} MB)")

    # --- Merged QC overview figure ---
    banner("MERGED DATASET OVERVIEW")
    puck_colors = {
        'Puck_211214_29': COLOR_PUCK_29,
        'Puck_211214_37': COLOR_PUCK_37,
        'Puck_211214_40': COLOR_PUCK_40,
    }

    fig, axes = plt.subplots(1, 3, figsize=(21, 6))
    for puck_name, adata_qc in adatas_qc.items():
        c = puck_colors.get(puck_name, '#999999')
        axes[0].hist(adata_qc.obs['total_counts'], bins=80, alpha=0.5,
                     label=puck_name, color=c, edgecolor='none')
        axes[1].hist(adata_qc.obs['n_genes_by_counts'], bins=80, alpha=0.5,
                     label=puck_name, color=c, edgecolor='none')
        axes[2].hist(adata_qc.obs['pct_counts_mt'], bins=40, alpha=0.5,
                     label=puck_name, color=c, edgecolor='none')

    axes[0].set_xlabel('Total UMIs'); axes[0].set_title('UMI distribution (post-QC)')
    axes[1].set_xlabel('Genes detected'); axes[1].set_title('Genes distribution (post-QC)')
    axes[2].set_xlabel('MT %'); axes[2].set_title('MT% distribution (post-QC)')
    for ax in axes:
        ax.legend(fontsize=FONT_SIZE - 10)
        ax.set_ylabel('Cells')
    plt.tight_layout()
    save_fig(fig, 'merged_QC_distributions', qc_fig_dir)

    # --- Final summary ---
    banner("STEP 03 SUMMARY")
    total_before = sum(r['cells_before'] for r in qc_rows)
    total_after = sum(r['cells_after'] for r in qc_rows)
    log(f"  Total cells: {total_before:,} -> {total_after:,} "
        f"({100 * total_after / total_before:.1f}% retained)")
    for r in qc_rows:
        log(f"    {r['puck']}: {r['cells_before']:,} -> {r['cells_after']:,} "
            f"({r['pct_retained']:.1f}%)")
    log(f"  Merged genes (inner join): {adata_merged.n_vars:,}")

    banner("STEP 03 COMPLETE")
    log(f"  Output: {DIR_03_QC}/")
    log(f"  Next: Run Step04_Cell_Type_Annotation.py")


if __name__ == "__main__":
    main()
