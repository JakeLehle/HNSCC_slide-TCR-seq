#!/usr/bin/env python3
"""
Step02_AnnData_Prep.py

Load Sophia Liu's pre-processed Slide-seq data into AnnData format.

Since STAR re-alignment from BCL is not possible (missing H52J2DMXY flow cell
data for barcode correction), this script reads directly from Sophia's
Drop-seq tools outputs:
  - Expression: *.matched.digital_expression.txt.gz (gene x barcode matrix)
  - Coordinates: barcode_matching/*_barcode_matching.txt.gz (barcode to x,y)

Outputs per-puck h5ad files with spatial coordinates to 02_anndata/.

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pathlib import Path
import warnings

# Add scripts directory to path for config import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spatial_config import (
    FASTQ_DIR, DIR_02_ANNDATA, DIR_TROUBLESHOOTING,
    PUCKS, PUCK_NAMES, puck_input_dir, puck_barcode_matching, puck_expression_file,
    puck_bam_file, RANDOM_SEED,
    banner, log, ensure_dir, save_fig
)


# =========================================================================
# DATA LOADING FUNCTIONS
# =========================================================================

def load_digital_expression(filepath):
    """
    Load Drop-seq tools digital expression matrix.

    Format: tab-separated, first row = GENE + cell barcodes,
    subsequent rows = gene_name + counts. Matrix is genes (rows) x cells (cols).
    Returns transposed cells x genes for AnnData.
    """
    log(f"  Loading expression: {Path(filepath).name}")

    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)

    genes = df.index.tolist()
    barcodes = df.columns.tolist()

    log(f"    Shape: {len(genes):,} genes x {len(barcodes):,} cells")
    log(f"    First 3 genes: {genes[:3]}")
    log(f"    First 3 barcodes: {barcodes[:3]}")

    return df, genes, barcodes


def load_barcode_coordinates(filepath):
    """
    Load barcode to spatial coordinate mapping.

    barcode_matching.txt.gz format (4 tab-separated columns, no header):
      Col1: Observed barcode (may have errors)
      Col2: Corrected barcode (matched to puck position, has -1 suffix)
      Col3: X coordinate (microns)
      Col4: Y coordinate (microns)

    Corrected barcode (col2) matches the expression matrix barcode format.
    """
    log(f"  Loading coordinates: {Path(filepath).name}")

    coords_df = pd.read_csv(
        filepath, sep='\t', compression='gzip', header=None,
        names=['observed_barcode', 'corrected_barcode', 'x', 'y']
    )

    log(f"    Raw rows: {len(coords_df):,}")

    # Corrected barcode (with -1 suffix) matches expression matrix format
    coords_df['barcode'] = coords_df['corrected_barcode']

    # Deduplicate: multiple observed barcodes can map to same corrected barcode
    coords_df = coords_df.drop_duplicates(subset='barcode', keep='first')
    coords_df = coords_df.set_index('barcode')

    log(f"    Unique spatial positions: {len(coords_df):,}")

    return coords_df[['x', 'y', 'observed_barcode']]


def load_single_puck(puck_name):
    """
    Load one Slide-seq puck into an AnnData object with spatial coordinates.
    """
    puck_dir = Path(puck_input_dir(puck_name))
    info = PUCKS[puck_name]

    banner(f"Loading puck: {puck_name}", char="-")
    log(f"  Directory: {puck_dir}")
    log(f"  Description: {info['description']}")

    # Resolve file paths
    expr_file = puck_expression_file(puck_name)
    coords_file = puck_barcode_matching(puck_name)
    bam_file = puck_bam_file(puck_name)

    for label, path in [("Expression", expr_file), ("Coordinates", coords_file)]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} file not found: {path}")

    # --- Load expression matrix ---
    expr_df, genes, barcodes = load_digital_expression(expr_file)

    # --- Load spatial coordinates ---
    coords_df = load_barcode_coordinates(coords_file)

    # --- Build AnnData (cells x genes) ---
    log(f"  Creating AnnData (transposing genes x cells to cells x genes)...")
    adata = ad.AnnData(
        X=expr_df.T.values,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes),
    )
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # --- Match barcodes between expression and coordinates ---
    expr_bc = set(adata.obs_names)
    coord_bc = set(coords_df.index)
    common_bc = expr_bc & coord_bc

    log(f"  Barcode matching:")
    log(f"    Expression barcodes: {len(expr_bc):,}")
    log(f"    Coordinate barcodes: {len(coord_bc):,}")
    log(f"    Common barcodes:     {len(common_bc):,}")

    if len(common_bc) == 0:
        log(f"    Expression examples: {list(expr_bc)[:5]}")
        log(f"    Coordinate examples: {list(coord_bc)[:5]}")
        raise ValueError(f"No matching barcodes for {puck_name}")

    match_rate = len(common_bc) / len(expr_bc)
    if match_rate < 0.95:
        warnings.warn(f"Only {match_rate:.1%} of expression barcodes have coordinates")
    log(f"    Match rate: {match_rate:.1%}")

    # Filter to common barcodes (preserve expression matrix order)
    common_list = [bc for bc in adata.obs_names if bc in common_bc]
    adata = adata[common_list, :].copy()

    # --- Add spatial coordinates ---
    spatial = np.array([[coords_df.loc[bc, 'x'], coords_df.loc[bc, 'y']]
                        for bc in adata.obs_names])
    adata.obsm['spatial'] = spatial
    adata.obs['x_coord'] = spatial[:, 0]
    adata.obs['y_coord'] = spatial[:, 1]
    adata.obs['observed_barcode'] = [coords_df.loc[bc, 'observed_barcode']
                                      for bc in adata.obs_names]

    # --- Puck metadata ---
    adata.obs['puck_id'] = puck_name
    adata.uns['puck_id'] = puck_name
    adata.uns['sample_barcode'] = info['sample_barcode']
    adata.uns['description'] = info['description']
    adata.uns['data_source'] = 'Sophia Liu, Ragon Institute'
    adata.uns['technology'] = 'Slide-TCR-seq'
    adata.uns['spatial'] = {
        puck_name: {
            'images': {},
            'scalefactors': {
                'tissue_hires_scalef': 1.0,
                'spot_diameter_fullres': 10.0,   # Slide-seq beads ~10 um
            }
        }
    }
    adata.uns['files'] = {
        'expression': str(expr_file),
        'coordinates': str(coords_file),
        'bam': str(bam_file) if os.path.exists(bam_file) else None,
    }

    # --- QC metrics ---
    log(f"  Computing QC metrics...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None,
                                log1p=False, inplace=True)

    # --- Summary ---
    log(f"  Final AnnData: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log(f"    Spatial X: [{spatial[:, 0].min():.0f}, {spatial[:, 0].max():.0f}] um")
    log(f"    Spatial Y: [{spatial[:, 1].min():.0f}, {spatial[:, 1].max():.0f}] um")
    log(f"    Median UMIs/cell: {adata.obs['total_counts'].median():.0f}")
    log(f"    Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    log(f"    Median MT%: {adata.obs['pct_counts_mt'].median():.1f}%")

    return adata


# =========================================================================
# MAIN
# =========================================================================

def main():
    banner("STEP 02: ANNDATA PREPARATION FROM SOPHIA'S PROCESSED DATA")
    log(f"  Input:  {FASTQ_DIR}")
    log(f"  Output: {DIR_02_ANNDATA}")

    ensure_dir(DIR_02_ANNDATA)
    ensure_dir(DIR_TROUBLESHOOTING)

    # --- Load each puck ---
    adatas = {}
    for puck_name in PUCK_NAMES:
        try:
            adata = load_single_puck(puck_name)
            adatas[puck_name] = adata
        except Exception as e:
            log(f"  ERROR loading {puck_name}: {e}")
            import traceback
            traceback.print_exc()

    if not adatas:
        log("FATAL: No pucks loaded. Check input paths.")
        sys.exit(1)

    # --- Save per-puck h5ad files ---
    banner("SAVING PER-PUCK H5AD FILES")
    for puck_name, adata in adatas.items():
        out_path = os.path.join(DIR_02_ANNDATA, f"{puck_name}.h5ad")
        adata.write_h5ad(out_path)
        size_mb = os.path.getsize(out_path) / 1e6
        log(f"  Saved: {out_path} ({size_mb:.1f} MB)")

    # --- Summary report ---
    banner("SUMMARY: ALL PUCKS")

    total_cells = 0
    total_umis = 0.0
    rows = []

    for puck_name, adata in adatas.items():
        n = adata.n_obs
        total_cells += n
        umi_sum = float(adata.X.sum())
        total_umis += umi_sum

        row = {
            'puck': puck_name,
            'cells': n,
            'genes': adata.n_vars,
            'total_umis': umi_sum,
            'median_umis': adata.obs['total_counts'].median(),
            'median_genes': adata.obs['n_genes_by_counts'].median(),
            'median_mt_pct': adata.obs['pct_counts_mt'].median(),
            'x_min': adata.obs['x_coord'].min(),
            'x_max': adata.obs['x_coord'].max(),
            'y_min': adata.obs['y_coord'].min(),
            'y_max': adata.obs['y_coord'].max(),
        }
        rows.append(row)

        log(f"  {puck_name}:")
        log(f"    {n:,} cells, {adata.n_vars:,} genes, {umi_sum:,.0f} UMIs")
        log(f"    Median: {row['median_umis']:.0f} UMIs, "
            f"{row['median_genes']:.0f} genes, {row['median_mt_pct']:.1f}% MT")

    log(f"")
    log(f"  Combined: {total_cells:,} cells, {total_umis:,.0f} UMIs")

    # Save summary TSV
    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(DIR_02_ANNDATA, "Step02_puck_summary.tsv")
    summary_df.to_csv(summary_path, sep='\t', index=False)
    log(f"  Summary table: {summary_path}")

    banner("STEP 02 COMPLETE")
    log(f"  Output: {DIR_02_ANNDATA}/")
    log(f"  Next: Run Step03_QC_Filter.py")


if __name__ == "__main__":
    main()
