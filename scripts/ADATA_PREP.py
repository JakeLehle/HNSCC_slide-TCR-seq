#!/usr/bin/env python3
"""
Load Slide-seq data into AnnData format for Squidpy analysis.

Based on data from Sophia Liu (Ragon Institute) - HPV+ HNSCC spatial transcriptomics.

File specifications from Sophia:
- Expression: XYZ.matched.digital_expression.txt.gz (gene x barcode matrix)
- Coordinates: barcode_matching/XYZ_barcode_matching.txt.gz (barcode to x,y mapping)
- BAM: XYZ.matched.bam (for HPV alignment)

Pucks:
- Puck_211214_29 (sample barcode: AGATTTAA)
- Puck_211214_37 (sample barcode: GGCGTCGA)
- Puck_211214_40 (sample barcode: ATCACTCG) - Best RNA recovery per Sophia

Author: Jake Lehle
Project: Slide-TCR-seq HPV+ HNSCC analysis
Date: December 2025
"""

import os
import gzip
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
from typing import Tuple, Optional, Dict
import warnings


def load_digital_expression(filepath: str) -> Tuple[pd.DataFrame, list, list]:
    """
    Load the digital expression matrix from Drop-seq tools output.
    
    The .digital_expression.txt.gz format from Drop-seq tools is:
    - First row: GENE followed by cell barcodes (tab-separated)
    - Subsequent rows: gene_name followed by counts
    - Matrix is genes (rows) x cells (columns)
    
    Parameters
    ----------
    filepath : str
        Path to the .matched.digital_expression.txt.gz file
        
    Returns
    -------
    matrix : pd.DataFrame
        Expression matrix (genes x cells)
    genes : list
        Gene names
    barcodes : list
        Cell barcodes
    """
    print(f"  Loading expression matrix from: {Path(filepath).name}")
    
    # Read the file - it's a tab-separated matrix with genes as rows
    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)
    
    genes = df.index.tolist()
    barcodes = df.columns.tolist()
    
    print(f"    Shape: {len(genes)} genes × {len(barcodes)} cells")
    print(f"    First 3 genes: {genes[:3]}")
    print(f"    First 3 barcodes: {barcodes[:3]}")
    
    return df, genes, barcodes


def load_barcode_coordinates(filepath: str) -> pd.DataFrame:
    """
    Load barcode to spatial coordinate mapping.
    
    The barcode_matching.txt.gz format contains 4 tab-separated columns:
    Column 1: Observed barcode (may have errors)
    Column 2: Corrected barcode (matched to puck position) with -1 suffix
    Column 3: X coordinate (microns)
    Column 4: Y coordinate (microns)
    
    We use the corrected barcode (column 2) as the unique identifier.
    
    Parameters
    ----------
    filepath : str
        Path to the *_barcode_matching.txt.gz file
        
    Returns
    -------
    coords_df : pd.DataFrame
        DataFrame with corrected barcode (without -1 suffix) as index and x, y columns
    """
    print(f"  Loading spatial coordinates from: {Path(filepath).name}")
    
    # Read the coordinates file - no header
    coords_df = pd.read_csv(
        filepath, 
        sep='\t', 
        compression='gzip',
        header=None,
        names=['observed_barcode', 'corrected_barcode', 'x', 'y']
    )
    
    print(f"    Shape: {coords_df.shape}")
    print(f"    Columns: {coords_df.columns.tolist()}")
    print(f"    Preview:")
    print(coords_df.head(3).to_string())
    
    # The corrected_barcode already has -1 suffix, which matches the expression matrix format
    # So we use it directly as the barcode identifier
    coords_df['barcode'] = coords_df['corrected_barcode']
    
    # Keep only unique corrected barcodes (some observed barcodes map to same corrected barcode)
    # This is the error-correction step that Drop-seq tools performed
    coords_df = coords_df.drop_duplicates(subset='barcode', keep='first')
    
    # Set barcode as index
    coords_df = coords_df.set_index('barcode')
    
    print(f"    Unique spatial positions: {len(coords_df)}")
    print(f"    First 5 barcodes: {coords_df.index.tolist()[:5]}")
    
    return coords_df[['x', 'y', 'observed_barcode']]


def load_slideseq_puck(puck_dir: str, puck_name: str) -> ad.AnnData:
    """
    Load a single Slide-seq puck into an AnnData object.
    
    Uses file specifications from Sophia Liu:
    - Expression: *.matched.digital_expression.txt.gz
    - Coordinates: barcode_matching/*_barcode_matching.txt.gz
    
    Parameters
    ----------
    puck_dir : str
        Path to the puck directory (e.g., '2022-01-28_Puck_211214_29')
    puck_name : str
        Puck identifier (e.g., 'Puck_211214_29')
        
    Returns
    -------
    adata : AnnData
        AnnData object with expression matrix and spatial coordinates
    """
    puck_dir = Path(puck_dir)
    
    print("=" * 60)
    print(f"Loading puck: {puck_name}")
    print("=" * 60)
    
    # Define file paths per Sophia's specification
    expr_file = puck_dir / f"{puck_name}.matched.digital_expression.txt.gz"
    coords_file = puck_dir / "barcode_matching" / f"{puck_name}_barcode_matching.txt.gz"
    bam_file = puck_dir / f"{puck_name}.matched.bam"
    
    # Check files exist
    if not expr_file.exists():
        raise FileNotFoundError(f"Expression file not found: {expr_file}")
    if not coords_file.exists():
        raise FileNotFoundError(f"Coordinates file not found: {coords_file}")
    
    # Load expression matrix
    expr_df, genes, barcodes = load_digital_expression(expr_file)
    
    # Load spatial coordinates
    coords_df = load_barcode_coordinates(coords_file)
    
    # Create AnnData (cells x genes, so transpose the expression matrix)
    # expr_df is genes x cells, we need cells x genes for AnnData
    print(f"\n  Creating AnnData object...")
    X = expr_df.T.values  # Transpose to cells x genes
    
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes)
    )
    
    # Make names unique (shouldn't be necessary but good practice)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    
    # Match barcodes between expression and coordinates
    expr_barcodes = set(adata.obs_names)
    coord_barcodes = set(coords_df.index)
    common_barcodes = expr_barcodes.intersection(coord_barcodes)
    
    print(f"\n  Barcode matching:")
    print(f"    Expression barcodes: {len(expr_barcodes):,}")
    print(f"    Coordinate barcodes: {len(coord_barcodes):,}")
    print(f"    Common barcodes: {len(common_barcodes):,}")
    
    if len(common_barcodes) == 0:
        print("\n  ERROR: No direct barcode matches!")
        print(f"    Expression barcode examples: {list(expr_barcodes)[:5]}")
        print(f"    Coordinate barcode examples: {list(coord_barcodes)[:5]}")
        raise ValueError("No matching barcodes - check barcode format differences")
    
    # Check for significant mismatch
    match_rate = len(common_barcodes) / len(expr_barcodes)
    if match_rate < 0.95:
        warnings.warn(f"Only {match_rate:.1%} of expression barcodes have coordinates")
    
    # Filter to common barcodes and maintain order
    common_barcodes_list = [bc for bc in adata.obs_names if bc in common_barcodes]
    adata = adata[common_barcodes_list, :].copy()
    
    # Add spatial coordinates to obsm
    print(f"\n  Adding spatial coordinates...")
    spatial_coords = np.zeros((adata.n_obs, 2))
    for i, bc in enumerate(adata.obs_names):
        spatial_coords[i, 0] = coords_df.loc[bc, 'x']
        spatial_coords[i, 1] = coords_df.loc[bc, 'y']
    
    adata.obsm['spatial'] = spatial_coords
    
    # Add coordinate metadata to obs
    adata.obs['x_coord'] = spatial_coords[:, 0]
    adata.obs['y_coord'] = spatial_coords[:, 1]
    adata.obs['observed_barcode'] = [coords_df.loc[bc, 'observed_barcode'] for bc in adata.obs_names]
    
    # Add puck metadata to uns
    adata.uns['puck_id'] = puck_name
    adata.uns['data_source'] = 'Sophia Liu, Ragon Institute'
    adata.uns['technology'] = 'Slide-TCR-seq'
    adata.uns['spatial'] = {
        puck_name: {
            'images': {},  # No images provided
            'scalefactors': {
                'tissue_hires_scalef': 1.0,
                'spot_diameter_fullres': 10.0  # Slide-seq beads are ~10 microns
            }
        }
    }
    
    # Store file paths for downstream analysis
    adata.uns['files'] = {
        'expression': str(expr_file),
        'coordinates': str(coords_file),
        'bam': str(bam_file) if bam_file.exists() else None
    }
    
    # Basic QC metrics
    print(f"\n  Computing QC metrics...")
    
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=['mt'], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    # Summary statistics
    print(f"\n  Final AnnData summary:")
    print(f"    Cells: {adata.n_obs:,}")
    print(f"    Genes: {adata.n_vars:,}")
    print(f"    Spatial extent:")
    print(f"      X: [{spatial_coords[:, 0].min():.1f}, {spatial_coords[:, 0].max():.1f}] microns")
    print(f"      Y: [{spatial_coords[:, 1].min():.1f}, {spatial_coords[:, 1].max():.1f}] microns")
    print(f"    Total counts per cell:")
    print(f"      Median: {adata.obs['total_counts'].median():.0f}")
    print(f"      Mean: {adata.obs['total_counts'].mean():.0f}")
    print(f"    Genes detected per cell:")
    print(f"      Median: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"      Mean: {adata.obs['n_genes_by_counts'].mean():.0f}")
    print(f"    Mitochondrial %:")
    print(f"      Median: {adata.obs['pct_counts_mt'].median():.1f}%")
    
    return adata


def load_all_pucks(base_dir: str) -> Dict[str, ad.AnnData]:
    """
    Load all three pucks into separate AnnData objects.
    
    Parameters
    ----------
    base_dir : str
        Base directory containing puck subdirectories
        
    Returns
    -------
    adatas : dict
        Dictionary with puck IDs as keys and AnnData objects as values
        Keys: 'Puck_211214_29', 'Puck_211214_37', 'Puck_211214_40'
    """
    base_dir = Path(base_dir)
    adatas = {}
    
    # Expected puck directories with metadata
    puck_info = {
        'Puck_211214_29': {
            'directory': '2022-01-28_Puck_211214_29',
            'sample_barcode': 'AGATTTAA',
            'description': 'HPV+ HNSCC sample 29'
        },
        'Puck_211214_37': {
            'directory': '2022-01-28_Puck_211214_37',
            'sample_barcode': 'GGCGTCGA',
            'description': 'HPV+ HNSCC sample 37'
        },
        'Puck_211214_40': {
            'directory': '2022-01-28_Puck_211214_40',
            'sample_barcode': 'ATCACTCG',
            'description': 'HPV+ HNSCC sample 40 (best RNA recovery)'
        }
    }
    
    print("\n" + "=" * 70)
    print("LOADING ALL PUCKS")
    print("=" * 70 + "\n")
    
    for puck_name, info in puck_info.items():
        puck_dir = base_dir / info['directory']
        
        if not puck_dir.exists():
            print(f"\n✗ Directory not found: {puck_dir}")
            continue
        
        try:
            adata = load_slideseq_puck(puck_dir, puck_name)
            
            # Add sample metadata
            adata.uns['sample_barcode'] = info['sample_barcode']
            adata.uns['description'] = info['description']
            
            adatas[puck_name] = adata
            print(f"\n✓ Successfully loaded {puck_name}")
            
        except Exception as e:
            print(f"\n✗ ERROR loading {puck_name}: {e}")
            import traceback
            traceback.print_exc()
    
    return adatas


def save_pucks(adatas: Dict[str, ad.AnnData], output_dir: str = ".") -> None:
    """
    Save each puck as a separate h5ad file.
    
    Parameters
    ----------
    adatas : dict
        Dictionary of AnnData objects
    output_dir : str
        Directory to save files to
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n" + "=" * 70)
    print("SAVING FILES")
    print("=" * 70)
    
    for puck_id, adata in adatas.items():
        output_file = output_dir / f"{puck_id}_anndata.h5ad"
        adata.write_h5ad(output_file)
        print(f"  ✓ Saved: {output_file}")
        print(f"    Size: {output_file.stat().st_size / 1e6:.1f} MB")


def print_summary(adatas: Dict[str, ad.AnnData]) -> None:
    """
    Print summary statistics for all loaded pucks.
    
    Parameters
    ----------
    adatas : dict
        Dictionary of AnnData objects
    """
    print("\n" + "=" * 70)
    print("SUMMARY - ALL PUCKS")
    print("=" * 70)
    
    for puck_id, adata in adatas.items():
        print(f"\n{puck_id}:")
        print(f"  Description: {adata.uns['description']}")
        print(f"  Sample barcode: {adata.uns['sample_barcode']}")
        print(f"  Cells: {adata.n_obs:,}")
        print(f"  Genes: {adata.n_vars:,}")
        print(f"  Total UMIs: {adata.X.sum():.0f}")
        print(f"  Median UMIs/cell: {adata.obs['total_counts'].median():.0f}")
        print(f"  Median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
        print(f"  Spatial extent: X=[{adata.obs['x_coord'].min():.0f}, {adata.obs['x_coord'].max():.0f}], "
              f"Y=[{adata.obs['y_coord'].min():.0f}, {adata.obs['y_coord'].max():.0f}]")
    
    # Combined statistics
    total_cells = sum(adata.n_obs for adata in adatas.values())
    total_umis = sum(adata.X.sum() for adata in adatas.values())
    
    print(f"\nCombined:")
    print(f"  Total cells: {total_cells:,}")
    print(f"  Total UMIs: {total_umis:.0f}")


# ============================================================================
# MAIN ANALYSIS WORKFLOW
# ============================================================================

if __name__ == "__main__":
    
    # UPDATE THIS PATH to your data location
    DATA_DIR = "/work/sdz852/WORKING/slide-TCR-seq/fastq"
    OUTPUT_DIR = "/work/sdz852/WORKING/slide-TCR-seq/analysis/anndata"
    
    print("\n" + "=" * 70)
    print("SLIDE-TCR-SEQ DATA LOADING - HPV+ HNSCC")
    print("Sophia Liu (Ragon Institute)")
    print("=" * 70)
    
    # Load all three pucks as separate objects
    adatas = load_all_pucks(DATA_DIR)
    
    if not adatas:
        print("\n✗ ERROR: No pucks were successfully loaded!")
        print("Check that DATA_DIR is correct and contains the expected directories.")
        exit(1)
    
    # Print summary
    print_summary(adatas)
    
    # Save each puck separately
    save_pucks(adatas, OUTPUT_DIR)
    
    print("\n" + "=" * 70)
    print("COMPLETE!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Load saved h5ad files in Jupyter/Python:")
    print("     adata = sc.read_h5ad('Puck_211214_40_anndata.h5ad')")
    print("  2. Visualize spatial distribution:")
    print("     import squidpy as sq")
    print("     sq.pl.spatial_scatter(adata, color='total_counts')")
    print("  3. Cluster and annotate cell types")
    print("  4. Align matched.bam files to HPV16 genome to identify HPV+ regions")
    print("  5. Spatial analysis of HPV+ vs HPV- tumor regions")
    print("\nData sources:")
    for puck_id, adata in adatas.items():
        print(f"  {puck_id}:")
        print(f"    Expression: {adata.uns['files']['expression']}")
        print(f"    BAM: {adata.uns['files']['bam']}")