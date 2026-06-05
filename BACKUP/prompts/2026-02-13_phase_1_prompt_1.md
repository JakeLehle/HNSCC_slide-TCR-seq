# Phase 1: Cell Type Annotation with popV

## Context

This phase performs automated cell type annotation on Slide-TCR-seq spatial transcriptomics 
data from HPV+ head and neck squamous cell carcinoma (HNSCC). Cell type annotation is the 
**BLOCKING** foundation for all downstream analyses - no other phase can proceed until this 
completes successfully.

We use popV with the pre-trained Tabula Sapiens model from HuggingFace, which provides 
comprehensive human cell type annotations without requiring tissue-specific reference datasets.
After popV annotation, we perform Leiden clustering and assign the dominant cell type to each 
cluster using a weighted scoring approach.

**CRITICAL DESIGN DECISION:** Process each puck INDIVIDUALLY, not the combined dataset.
- There are 3 pucks: Puck_211214_29, Puck_211214_37, Puck_211214_40
- Each puck has ~30K cells after QC
- Individual processing is essential for:
  * Spatial visualization (spatial coordinates are puck-specific)
  * Avoiding batch effects during annotation
  * Easier debugging and quality control
  * Downstream spatial analyses that require per-puck organization
- Run the complete pipeline on each puck sequentially
- Merge annotations at the end if needed for cross-puck comparisons

**Reference Script:** `scripts/example_reference_scripts/SC_Cluster_Annotation.py`
This script demonstrates the exact workflow for scRNA-seq data. Adapt it for spatial data 
by preserving spatial coordinates and using `puck_id` context appropriately.

## Environment

```yaml
name: phase1_celltype
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - scanpy>=1.9.0
  - squidpy>=1.2.0
  - anndata>=0.9.0
  - popv>=0.9.0
  - bbknn>=1.5.0
  - pandas>=2.0.0
  - numpy>=1.24.0
  - matplotlib>=3.7.0
  - seaborn>=0.12.0
  - scipy>=1.10.0
  - numba>=0.57.0
  - adjustText>=0.8
```

**Data Location:** `data/outputs/analysis/processed/`
**Output Location:** `data/outputs/analysis/annotated/`
**Figures Location:** `data/outputs/analysis/figures/`
**Scripts Location:** `scripts/`

## Input Files

| File | Path | Description |
|------|------|-------------|
| Puck 29 | `data/outputs/analysis/processed/Puck_211214_29_processed.h5ad` | QC'd spatial data (~30K cells) |
| Puck 37 | `data/outputs/analysis/processed/Puck_211214_37_processed.h5ad` | QC'd spatial data (~30K cells) |
| Puck 40 | `data/outputs/analysis/processed/Puck_211214_40_processed.h5ad` | QC'd spatial data (~30K cells) |
| Reference Script | `scripts/example_reference_scripts/SC_Cluster_Annotation.py` | Working annotation pipeline |

**Input Data Structure Expected:**
- `adata.X` contains raw counts (required for popV)
- `adata.obsm['spatial']` contains x,y coordinates
- `adata.var` has gene symbols (check column name: 'gene_symbol', 'feature_name', or index)
- ~30K cells per puck after prior QC filtering

## Steps

### Step 1: Validate Input Data and Setup Environment

**Goal:** Confirm input files exist, have correct structure, and set up output directories.

**Input:** 
- `data/outputs/analysis/processed/Puck_211214_29_processed.h5ad`
- `data/outputs/analysis/processed/Puck_211214_37_processed.h5ad`
- `data/outputs/analysis/processed/Puck_211214_40_processed.h5ad`

**Output:**
- `data/outputs/analysis/annotated/` directory created
- `data/outputs/analysis/figures/phase1/` directory created
- Validation log confirming data structure

**Approach:**
```python
import scanpy as sc
import squidpy as sq
import os
import sys

# Define puck IDs to process
PUCK_IDS = ['Puck_211214_29', 'Puck_211214_37', 'Puck_211214_40']
INPUT_DIR = 'data/outputs/analysis/processed'
OUTPUT_DIR = 'data/outputs/analysis/annotated'
FIGURES_DIR = 'data/outputs/analysis/figures/phase1'

# Create output directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Validate each puck file
for puck_id in PUCK_IDS:
    filepath = os.path.join(INPUT_DIR, f'{puck_id}_processed.h5ad')
    
    if not os.path.exists(filepath):
        print(f"ERROR: Missing file {filepath}")
        sys.exit(1)
    
    adata = sc.read_h5ad(filepath)
    print(f"\n{puck_id}:")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")
    print(f"  Has spatial coords: {'spatial' in adata.obsm}")
    print(f"  X matrix type: {type(adata.X)}")
    print(f"  X contains raw counts: {adata.X.max() > 100}")  # Raw counts typically have high max
    
    # Check for gene symbol column
    if 'gene_symbol' in adata.var.columns:
        print(f"  Gene symbols in: var['gene_symbol']")
    elif 'feature_name' in adata.var.columns:
        print(f"  Gene symbols in: var['feature_name']")
    else:
        print(f"  Gene symbols in: var.index")
```

**Success Criteria:**
- All 3 puck files exist and are readable
- Each puck has >20,000 cells
- Spatial coordinates present in `adata.obsm['spatial']`
- Raw counts present (max value > 100)

---

### Step 2: Run popV Annotation (Per Puck)

**Goal:** Annotate cell types using the pre-trained Tabula Sapiens model.

**Input:** `data/outputs/analysis/processed/{puck_id}_processed.h5ad`

**Output:** `data/outputs/analysis/annotated/{puck_id}_popv_raw.h5ad`

**CRITICAL:** popV requires RAW counts in `adata.X`. If your data is normalized, you must 
either reload the raw data or use `adata.raw` if it was preserved during QC.

**Approach:**
```python
import popv
import scanpy as sc

def run_popv_annotation(puck_id):
    """
    Run popV cell type annotation on a single puck.
    """
    input_path = f'data/outputs/analysis/processed/{puck_id}_processed.h5ad'
    output_path = f'data/outputs/analysis/annotated/{puck_id}_popv_raw.h5ad'
    
    print(f"Processing {puck_id}...")
    
    # Load data
    adata = sc.read_h5ad(input_path)
    
    # IMPORTANT: Ensure we have raw counts
    # If adata.X is normalized, check for adata.raw
    if hasattr(adata, 'raw') and adata.raw is not None:
        print("  Using adata.raw for popV (contains raw counts)")
        adata_for_popv = adata.raw.to_adata()
        # Preserve spatial coordinates
        adata_for_popv.obsm['spatial'] = adata.obsm['spatial']
    else:
        print("  Using adata.X directly (assuming raw counts)")
        adata_for_popv = adata.copy()
    
    # Determine gene symbol column
    if 'gene_symbol' in adata_for_popv.var.columns:
        gene_col = 'gene_symbol'
    elif 'feature_name' in adata_for_popv.var.columns:
        gene_col = 'feature_name'
    else:
        gene_col = None  # Use index
        adata_for_popv.var['gene_symbol'] = adata_for_popv.var_names
        gene_col = 'gene_symbol'
    
    print(f"  Using gene column: {gene_col}")
    print(f"  Cells: {adata_for_popv.n_obs}, Genes: {adata_for_popv.n_vars}")
    
    # Load pre-trained Tabula Sapiens model
    huggingface_repo = "popV/tabula_sapiens_All_Cells"
    
    print("  Downloading/loading Tabula Sapiens model...")
    hmo = popv.hub.HubModel.pull_from_huggingface_hub(
        huggingface_repo,
        cache_dir="tmp/tabula_sapiens"
    )
    
    # Run annotation
    # Note: For spatial data, we don't have a batch key in the traditional sense
    # We'll use a dummy batch key or skip batch correction in popV
    print("  Running popV annotation (this may take 10-30 minutes)...")
    
    adata_annotated = hmo.annotate_data(
        adata_for_popv,
        query_batch_key=None,  # No batch correction within a single puck
        prediction_mode="inference",  # Full integration mode for best accuracy
        gene_symbols=gene_col
    )
    
    # Filter to keep only query data (exclude reference cells if any)
    if '_dataset' in adata_annotated.obs.columns:
        adata_annotated = adata_annotated[adata_annotated.obs["_dataset"] == "query"].copy()
    
    # Verify popV outputs
    assert 'popv_prediction' in adata_annotated.obs.columns, "popV prediction missing!"
    assert 'popv_prediction_score' in adata_annotated.obs.columns, "popV score missing!"
    
    # Add puck_id to obs for tracking
    adata_annotated.obs['puck_id'] = puck_id
    
    # Save intermediate result
    adata_annotated.write_h5ad(output_path)
    print(f"  Saved: {output_path}")
    
    # Print summary
    print(f"\n  popV Results for {puck_id}:")
    print(f"  Unique cell types: {adata_annotated.obs['popv_prediction'].nunique()}")
    print(f"  Mean confidence: {adata_annotated.obs['popv_prediction_score'].mean():.3f}")
    print(f"  Top 5 cell types:")
    print(adata_annotated.obs['popv_prediction'].value_counts().head())
    
    return adata_annotated

# Process each puck
for puck_id in PUCK_IDS:
    run_popv_annotation(puck_id)
```

**Success Criteria:**
- `popv_prediction` column present in each output
- `popv_prediction_score` column present (values 0-1)
- Mean prediction score > 0.5
- Multiple cell types detected (>5 unique types)

---

### Step 3: Clustering with Leiden Algorithm (Per Puck)

**Goal:** Perform unsupervised clustering to group similar cells.

**Input:** `data/outputs/analysis/annotated/{puck_id}_popv_raw.h5ad`

**Output:** `data/outputs/analysis/annotated/{puck_id}_clustered.h5ad`

**Approach:**
```python
import scanpy as sc
from os import cpu_count

def cluster_puck(puck_id):
    """
    Normalize, reduce dimensions, and cluster a single puck.
    """
    input_path = f'data/outputs/analysis/annotated/{puck_id}_popv_raw.h5ad'
    output_path = f'data/outputs/analysis/annotated/{puck_id}_clustered.h5ad'
    
    print(f"Clustering {puck_id}...")
    
    # Load popV-annotated data
    adata = sc.read_h5ad(input_path)
    
    # Store raw counts before normalization (if not already stored)
    if adata.raw is None:
        adata.raw = adata.copy()
    
    # Normalize and log-transform
    print("  Normalizing...")
    sc.pp.normalize_total(adata, target_sum=1e6)  # CPM normalization
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    print("  Finding HVGs...")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', 
                                 span=0.3, subset=False)
    
    # PCA on HVGs
    print("  Running PCA...")
    n_pcs = min(50, adata.n_obs - 1, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_pcs, use_highly_variable=True)
    
    # Compute neighborhood graph
    print("  Computing neighbors...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, n_pcs))
    
    # UMAP embedding
    print("  Computing UMAP...")
    sc.tl.umap(adata, random_state=42)
    
    # Leiden clustering at multiple resolutions
    print("  Leiden clustering...")
    for resolution in [0.5, 1.0, 1.5]:
        sc.tl.leiden(adata, resolution=resolution, 
                     key_added=f'leiden_{resolution}', random_state=42)
    
    # Set default clusters to resolution 1.0
    adata.obs['clusters'] = adata.obs['leiden_1.0']
    
    # Print cluster summary
    print(f"\n  Clustering Results for {puck_id}:")
    print(f"  Number of clusters (res=1.0): {adata.obs['clusters'].nunique()}")
    print(f"  Cluster sizes:")
    print(adata.obs['clusters'].value_counts().head(10))
    
    # Save
    adata.write_h5ad(output_path)
    print(f"  Saved: {output_path}")
    
    return adata

# Process each puck
for puck_id in PUCK_IDS:
    cluster_puck(puck_id)
```

**Success Criteria:**
- 10-30 clusters identified (resolution 1.0)
- UMAP coordinates in `adata.obsm['X_umap']`
- `clusters` column in `adata.obs`
- No cluster contains >50% of all cells (avoid single dominant cluster)

---

### Step 4: Assign Cluster-Level Cell Types (Per Puck)

**Goal:** Assign the dominant popV cell type to each cluster using weighted scoring.

**Input:** `data/outputs/analysis/annotated/{puck_id}_clustered.h5ad`

**Output:** `data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad`

**Logic:** 
For each cluster, aggregate popV predictions weighted by confidence scores.
The cell type with highest weighted score becomes the cluster's annotation.
This smooths over noisy individual predictions while respecting confidence.

**Approach:**
```python
import scanpy as sc
import pandas as pd
import numpy as np

def assign_cluster_based_annotation(adata):
    """
    Assigns cell type to each cluster based on weighted popV predictions.
    
    Algorithm:
    1. Min-Max normalize popV scores within each cluster
    2. Apply linear weighting (scores sum to 1 per cluster)
    3. Calculate weighted score per cell type per cluster
    4. Assign dominant cell type to cluster
    5. Map back to all cells
    """
    df = adata.obs[['clusters', 'popv_prediction', 'popv_prediction_score']].copy()
    df['popv_prediction_score'] = pd.to_numeric(df['popv_prediction_score'], errors='coerce')
    
    # Min-Max normalization within each cluster
    def min_max_normalize(x):
        x_min, x_max = x.min(), x.max()
        if x_max - x_min == 0:
            return np.zeros(len(x))
        return (x - x_min) / (x_max - x_min)
    
    df['normalized_score'] = df.groupby('clusters')['popv_prediction_score'].transform(min_max_normalize)
    
    # Linear weighting (normalize to sum to 1 within each cluster)
    def linear_weights(x):
        total = x.sum()
        if total == 0:
            return np.ones(len(x)) / len(x)
        return x / total
    
    df['weight'] = df.groupby('clusters')['normalized_score'].transform(linear_weights)
    df['weighted_score'] = df['popv_prediction_score'] * df['weight']
    
    # Aggregate by cluster and cell type, find dominant type
    cluster_type_scores = df.groupby(['clusters', 'popv_prediction'])['weighted_score'].sum().reset_index()
    
    dominant_types = (
        cluster_type_scores
        .sort_values('weighted_score', ascending=False)
        .drop_duplicates('clusters')
        .set_index('clusters')['popv_prediction']
    )
    
    # Map dominant type back to all cells
    adata.obs['final_annotation'] = adata.obs['clusters'].map(dominant_types)
    
    # Calculate cluster-level confidence (mean score for dominant type)
    cluster_confidence = df.groupby('clusters').apply(
        lambda x: x[x['popv_prediction'] == dominant_types[x.name]]['popv_prediction_score'].mean()
    )
    adata.obs['cluster_confidence'] = adata.obs['clusters'].map(cluster_confidence)
    
    return adata


def annotate_puck(puck_id):
    """
    Apply cluster-based annotation to a single puck.
    """
    input_path = f'data/outputs/analysis/annotated/{puck_id}_clustered.h5ad'
    output_path = f'data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad'
    
    print(f"Annotating clusters for {puck_id}...")
    
    # Load clustered data
    adata = sc.read_h5ad(input_path)
    
    # Apply cluster-based annotation
    adata = assign_cluster_based_annotation(adata)
    
    # Convert final_annotation to categorical for plotting
    adata.obs['final_annotation'] = adata.obs['final_annotation'].astype('category')
    
    # Print annotation summary
    print(f"\n  Final Annotations for {puck_id}:")
    print(f"  Unique cell types: {adata.obs['final_annotation'].nunique()}")
    print(f"  Cell type distribution:")
    print(adata.obs['final_annotation'].value_counts())
    
    # Save
    adata.write_h5ad(output_path)
    print(f"  Saved: {output_path}")
    
    return adata

# Process each puck
for puck_id in PUCK_IDS:
    annotate_puck(puck_id)
```

**Success Criteria:**
- `final_annotation` column present with categorical cell types
- Expected cell types detected: T cells, B cells, myeloid, epithelial, fibroblasts
- No cluster left without annotation
- Cluster confidence scores reasonable (mean > 0.4)

---

### Step 5: Generate Visualization Figures (Per Puck)

**Goal:** Create QC and validation figures for each puck.

**Input:** `data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad`

**Output:** 
- `data/outputs/analysis/figures/phase1/{puck_id}_umap_celltypes.pdf`
- `data/outputs/analysis/figures/phase1/{puck_id}_spatial_celltypes.pdf`
- `data/outputs/analysis/figures/phase1/{puck_id}_popv_scores.pdf`
- `data/outputs/analysis/figures/phase1/{puck_id}_cluster_composition.pdf`
- `data/outputs/analysis/figures/phase1/{puck_id}_marker_dotplot.pdf`

**Approach:**
```python
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def generate_figures(puck_id):
    """
    Generate all visualization figures for a single puck.
    """
    input_path = f'data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad'
    figures_dir = f'data/outputs/analysis/figures/phase1'
    
    print(f"Generating figures for {puck_id}...")
    
    # Load annotated data
    adata = sc.read_h5ad(input_path)
    
    # Set up color palette
    n_types = adata.obs['final_annotation'].nunique()
    palette = sns.color_palette('turbo', n_types)
    
    # Figure 1: UMAP colored by cell type
    print("  Creating UMAP plot...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    sc.pl.umap(adata, color='final_annotation', ax=axes[0], show=False,
               title=f'{puck_id}: Cell Types', frameon=False)
    sc.pl.umap(adata, color='clusters', ax=axes[1], show=False,
               title=f'{puck_id}: Leiden Clusters', frameon=False, legend_loc='on data')
    
    plt.tight_layout()
    plt.savefig(f'{figures_dir}/{puck_id}_umap_celltypes.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 2: Spatial plot colored by cell type
    print("  Creating spatial plot...")
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    sq.pl.spatial_scatter(
        adata,
        color='final_annotation',
        size=1.5,
        ax=ax,
        title=f'{puck_id}: Spatial Cell Types'
    )
    
    plt.savefig(f'{figures_dir}/{puck_id}_spatial_celltypes.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 3: popV confidence scores
    print("  Creating confidence score plot...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # UMAP colored by score
    sc.pl.umap(adata, color='popv_prediction_score', ax=axes[0], show=False,
               title='popV Confidence Score', cmap='magma', frameon=False)
    
    # Histogram of scores
    axes[1].hist(adata.obs['popv_prediction_score'], bins=50, edgecolor='black', alpha=0.7)
    axes[1].axvline(x=0.5, color='red', linestyle='--', label='Threshold (0.5)')
    axes[1].set_xlabel('popV Prediction Score')
    axes[1].set_ylabel('Number of Cells')
    axes[1].set_title('Distribution of Confidence Scores')
    axes[1].legend()
    
    plt.tight_layout()
    plt.savefig(f'{figures_dir}/{puck_id}_popv_scores.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 4: Cluster composition stacked bar
    print("  Creating cluster composition plot...")
    
    # Calculate proportions
    cluster_celltype = pd.crosstab(
        adata.obs['clusters'], 
        adata.obs['popv_prediction'],
        normalize='index'
    ) * 100
    
    # Sort by most common cell type in cluster
    cluster_order = adata.obs.groupby('clusters')['final_annotation'].first().sort_values().index
    
    fig, ax = plt.subplots(figsize=(14, 8))
    cluster_celltype.loc[cluster_order].plot(
        kind='bar', stacked=True, ax=ax, 
        colormap='turbo', edgecolor='white', linewidth=0.5
    )
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Percentage (%)')
    ax.set_title(f'{puck_id}: Cell Type Composition per Cluster')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f'{figures_dir}/{puck_id}_cluster_composition.pdf', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Figure 5: Marker gene validation
    print("  Creating marker dotplot...")
    
    marker_genes = {
        'T_cells': ['CD3D', 'CD3E', 'CD4', 'CD8A'],
        'B_cells': ['CD19', 'MS4A1', 'CD79A'],
        'Myeloid': ['CD14', 'CD68', 'LYZ'],
        'Dendritic': ['CLEC9A', 'CD1C'],
        'Fibroblasts': ['COL1A1', 'DCN'],
        'Endothelial': ['PECAM1', 'VWF'],
        'Epithelial': ['EPCAM', 'KRT5', 'KRT14']
    }
    
    # Filter to genes present in data
    available_markers = {}
    for cell_type, genes in marker_genes.items():
        available = [g for g in genes if g in adata.var_names]
        if available:
            available_markers[cell_type] = available
    
    if available_markers:
        sc.pl.dotplot(
            adata,
            var_names=available_markers,
            groupby='final_annotation',
            save=f'_{puck_id}_marker_validation.pdf',
            show=False
        )
        # Move scanpy's auto-saved file to our directory
        import shutil
        import os
        src = f'figures/dotplot_{puck_id}_marker_validation.pdf'
        if os.path.exists(src):
            shutil.move(src, f'{figures_dir}/{puck_id}_marker_dotplot.pdf')
    
    print(f"  Figures saved to {figures_dir}/")
    
    return True

# Process each puck
for puck_id in PUCK_IDS:
    generate_figures(puck_id)
```

**Success Criteria:**
- All 5 figure types generated for each puck
- Spatial plots show biologically plausible cell type distributions
- Marker genes validate expected cell types
- Confidence scores predominantly > 0.5

---

### Step 6: Generate Summary Statistics and Report

**Goal:** Create summary tables and a markdown report for Phase 1.

**Input:** All `{puck_id}_celltype_annotated.h5ad` files

**Output:**
- `data/outputs/analysis/annotated/celltype_summary_all_pucks.csv`
- `reports/phase1_celltype_annotation_report.md`

**Approach:**
```python
import scanpy as sc
import pandas as pd
from datetime import datetime

def generate_summary():
    """
    Generate summary statistics across all pucks.
    """
    PUCK_IDS = ['Puck_211214_29', 'Puck_211214_37', 'Puck_211214_40']
    
    all_summaries = []
    
    for puck_id in PUCK_IDS:
        filepath = f'data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad'
        adata = sc.read_h5ad(filepath)
        
        # Per-puck summary
        summary = {
            'puck_id': puck_id,
            'total_cells': adata.n_obs,
            'n_clusters': adata.obs['clusters'].nunique(),
            'n_celltypes': adata.obs['final_annotation'].nunique(),
            'mean_popv_score': adata.obs['popv_prediction_score'].mean(),
            'median_popv_score': adata.obs['popv_prediction_score'].median(),
            'pct_high_confidence': (adata.obs['popv_prediction_score'] > 0.5).mean() * 100
        }
        
        # Add cell type counts
        for ct in adata.obs['final_annotation'].unique():
            summary[f'n_{ct}'] = (adata.obs['final_annotation'] == ct).sum()
            summary[f'pct_{ct}'] = (adata.obs['final_annotation'] == ct).mean() * 100
        
        all_summaries.append(summary)
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(all_summaries)
    summary_df.to_csv('data/outputs/analysis/annotated/celltype_summary_all_pucks.csv', index=False)
    
    # Generate markdown report
    report = f"""# Phase 1: Cell Type Annotation Report

**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
**Pipeline:** Slide-TCR-seq Spatial Transcriptomics Analysis

## Executive Summary

Cell type annotation completed for {len(PUCK_IDS)} pucks using popV with Tabula Sapiens model.

## Per-Puck Statistics

| Puck ID | Total Cells | Clusters | Cell Types | Mean Score | High Conf (%) |
|---------|-------------|----------|------------|------------|---------------|
"""
    
    for _, row in summary_df.iterrows():
        report += f"| {row['puck_id']} | {row['total_cells']:,} | {row['n_clusters']} | {row['n_celltypes']} | {row['mean_popv_score']:.3f} | {row['pct_high_confidence']:.1f}% |\n"
    
    report += """
## Cell Type Distribution

"""
    
    # Get all unique cell types across pucks
    all_celltypes = set()
    for puck_id in PUCK_IDS:
        filepath = f'data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad'
        adata = sc.read_h5ad(filepath)
        all_celltypes.update(adata.obs['final_annotation'].unique())
    
    report += "| Cell Type | " + " | ".join(PUCK_IDS) + " |\n"
    report += "|-----------|" + "|".join(["------" for _ in PUCK_IDS]) + "|\n"
    
    for ct in sorted(all_celltypes):
        row = f"| {ct} |"
        for puck_id in PUCK_IDS:
            col = f'pct_{ct}'
            if col in summary_df.columns:
                pct = summary_df[summary_df['puck_id'] == puck_id][col].values
                if len(pct) > 0 and not pd.isna(pct[0]):
                    row += f" {pct[0]:.1f}% |"
                else:
                    row += " - |"
            else:
                row += " - |"
        report += row + "\n"
    
    report += """
## Output Files

| File | Description |
|------|-------------|
"""
    
    for puck_id in PUCK_IDS:
        report += f"| `{puck_id}_celltype_annotated.h5ad` | Annotated AnnData with cell types |\n"
    
    report += """| `celltype_summary_all_pucks.csv` | Summary statistics table |

## Quality Control

### Confidence Score Distribution
- Target: >70% cells with score > 0.5
- See individual puck figures for detailed distributions

### Expected Cell Types
- ✓ T cells (CD3D+, CD4+ or CD8A+)
- ✓ B cells (CD19+, MS4A1+)
- ✓ Myeloid cells (CD14+, CD68+)
- ✓ Epithelial/Tumor cells (EPCAM+, KRT5+)
- ✓ Fibroblasts (COL1A1+, DCN+)
- ✓ Endothelial cells (PECAM1+, VWF+)

## Next Steps

Phase 1 outputs are required for:
- **Phase 2:** Spatial neighborhood analysis, TCR calling, TLS identification
- **Phase 3:** Mutation calling with SComatic (requires cell type labels)

## Files Generated

```
data/outputs/analysis/annotated/
├── Puck_211214_29_celltype_annotated.h5ad
├── Puck_211214_37_celltype_annotated.h5ad
├── Puck_211214_40_celltype_annotated.h5ad
└── celltype_summary_all_pucks.csv

data/outputs/analysis/figures/phase1/
├── Puck_211214_29_umap_celltypes.pdf
├── Puck_211214_29_spatial_celltypes.pdf
├── Puck_211214_29_popv_scores.pdf
├── Puck_211214_29_cluster_composition.pdf
├── Puck_211214_29_marker_dotplot.pdf
├── [same for Puck_211214_37]
└── [same for Puck_211214_40]
```
"""
    
    # Save report
    with open('reports/phase1_celltype_annotation_report.md', 'w') as f:
        f.write(report)
    
    print("Summary and report generated!")
    print("  - data/outputs/analysis/annotated/celltype_summary_all_pucks.csv")
    print("  - reports/phase1_celltype_annotation_report.md")

generate_summary()
```

**Success Criteria:**
- Summary CSV contains all pucks
- Report includes per-puck statistics
- All expected cell types documented

---

### Step 7: Create Checkpoint File and Verify Completion

**Goal:** Create checkpoint file indicating Phase 1 completion and verify all outputs.

**Input:** All outputs from previous steps

**Output:**
- `data/outputs/working/checkpoints/phase1_celltype_annotation.done`
- Verification log

**Approach:**
```python
import os
import scanpy as sc
from datetime import datetime

def verify_phase1_completion():
    """
    Verify all Phase 1 outputs exist and meet quality criteria.
    """
    PUCK_IDS = ['Puck_211214_29', 'Puck_211214_37', 'Puck_211214_40']
    
    errors = []
    warnings = []
    
    print("=" * 60)
    print("PHASE 1 COMPLETION VERIFICATION")
    print("=" * 60)
    
    # Check annotated h5ad files
    print("\n1. Checking annotated AnnData files...")
    for puck_id in PUCK_IDS:
        filepath = f'data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad'
        
        if not os.path.exists(filepath):
            errors.append(f"MISSING: {filepath}")
            continue
        
        adata = sc.read_h5ad(filepath)
        
        # Check required columns
        required_cols = ['popv_prediction', 'popv_prediction_score', 
                         'clusters', 'final_annotation', 'puck_id']
        for col in required_cols:
            if col not in adata.obs.columns:
                errors.append(f"{puck_id}: Missing obs column '{col}'")
        
        # Check cell counts
        if adata.n_obs < 20000:
            warnings.append(f"{puck_id}: Low cell count ({adata.n_obs})")
        
        # Check confidence scores
        mean_score = adata.obs['popv_prediction_score'].mean()
        if mean_score < 0.5:
            warnings.append(f"{puck_id}: Low mean confidence ({mean_score:.3f})")
        
        # Check spatial coordinates
        if 'spatial' not in adata.obsm:
            errors.append(f"{puck_id}: Missing spatial coordinates")
        
        print(f"  ✓ {puck_id}: {adata.n_obs} cells, {adata.obs['final_annotation'].nunique()} cell types")
    
    # Check figures
    print("\n2. Checking figure outputs...")
    expected_figures = [
        'umap_celltypes.pdf',
        'spatial_celltypes.pdf', 
        'popv_scores.pdf',
        'cluster_composition.pdf',
        'marker_dotplot.pdf'
    ]
    
    for puck_id in PUCK_IDS:
        for fig in expected_figures:
            filepath = f'data/outputs/analysis/figures/phase1/{puck_id}_{fig}'
            if not os.path.exists(filepath):
                warnings.append(f"MISSING FIGURE: {filepath}")
            else:
                print(f"  ✓ {puck_id}_{fig}")
    
    # Check summary files
    print("\n3. Checking summary files...")
    summary_files = [
        'data/outputs/analysis/annotated/celltype_summary_all_pucks.csv',
        'reports/phase1_celltype_annotation_report.md'
    ]
    
    for filepath in summary_files:
        if not os.path.exists(filepath):
            errors.append(f"MISSING: {filepath}")
        else:
            print(f"  ✓ {filepath}")
    
    # Report results
    print("\n" + "=" * 60)
    
    if errors:
        print("ERRORS (must fix):")
        for e in errors:
            print(f"  ✗ {e}")
    
    if warnings:
        print("\nWARNINGS (review recommended):")
        for w in warnings:
            print(f"  ⚠ {w}")
    
    if not errors:
        print("\n✓ PHASE 1 COMPLETE - All critical checks passed")
        
        # Create checkpoint file
        checkpoint_dir = 'data/outputs/working/checkpoints'
        os.makedirs(checkpoint_dir, exist_ok=True)
        
        checkpoint_file = f'{checkpoint_dir}/phase1_celltype_annotation.done'
        with open(checkpoint_file, 'w') as f:
            f.write(f"Phase 1 completed: {datetime.now().isoformat()}\n")
            f.write(f"Pucks processed: {', '.join(PUCK_IDS)}\n")
            for puck_id in PUCK_IDS:
                filepath = f'data/outputs/analysis/annotated/{puck_id}_celltype_annotated.h5ad'
                adata = sc.read_h5ad(filepath)
                f.write(f"{puck_id}: {adata.n_obs} cells, {adata.obs['final_annotation'].nunique()} cell types\n")
        
        print(f"\nCheckpoint created: {checkpoint_file}")
        return True
    else:
        print("\n✗ PHASE 1 INCOMPLETE - Fix errors before proceeding to Phase 2")
        return False

# Run verification
success = verify_phase1_completion()
```

**Success Criteria:**
- All 3 puck h5ad files exist with required columns
- All figures generated
- Summary files exist
- Checkpoint file created
- Zero critical errors

---

## Expected Outputs

### Primary Data Files
| File | Path | Description |
|------|------|-------------|
| Annotated Puck 29 | `data/outputs/analysis/annotated/Puck_211214_29_celltype_annotated.h5ad` | ~30K cells with annotations |
| Annotated Puck 37 | `data/outputs/analysis/annotated/Puck_211214_37_celltype_annotated.h5ad` | ~30K cells with annotations |
| Annotated Puck 40 | `data/outputs/analysis/annotated/Puck_211214_40_celltype_annotated.h5ad` | ~30K cells with annotations |
| Summary CSV | `data/outputs/analysis/annotated/celltype_summary_all_pucks.csv` | Statistics table |

### AnnData Structure (per puck)
```
adata.obs columns:
├── popv_prediction          # Raw popV cell type (string)
├── popv_prediction_score    # Confidence 0-1 (float)
├── clusters                 # Leiden cluster ID (string)
├── leiden_0.5              # Low resolution clusters
├── leiden_1.0              # Medium resolution clusters  
├── leiden_1.5              # High resolution clusters
├── final_annotation        # Cluster-level cell type (category)
├── cluster_confidence      # Mean confidence for cluster annotation
└── puck_id                 # Puck identifier

adata.obsm:
├── spatial                 # x,y coordinates (preserved from input)
├── X_pca                  # PCA embedding
└── X_umap                 # UMAP embedding

adata.var:
├── highly_variable        # HVG flag
└── gene_symbol           # Gene names
```

### Figures (per puck)
| Figure | Description |
|--------|-------------|
| `{puck}_umap_celltypes.pdf` | UMAP colored by cell type and clusters |
| `{puck}_spatial_celltypes.pdf` | Spatial scatter plot of cell types |
| `{puck}_popv_scores.pdf` | Confidence score distribution |
| `{puck}_cluster_composition.pdf` | Stacked bar of cell types per cluster |
| `{puck}_marker_dotplot.pdf` | Marker gene validation |

### Reports
| File | Description |
|------|-------------|
| `reports/phase1_celltype_annotation_report.md` | Full summary report |

### Checkpoints
| File | Description |
|------|-------------|
| `data/outputs/working/checkpoints/phase1_celltype_annotation.done` | Completion marker |

---

## Dependencies Between Steps

```
Step 1 (Validate) ──────────────────────────────────┐
                                                     │
Step 2 (popV) ◄──────────────────────────────────────┤
     │                                               │
     │ Requires raw counts                           │
     ▼                                               │
Step 3 (Clustering) ◄─────────────── Requires Step 2 │
     │                                               │
     │ Requires normalized + clustered               │
     ▼                                               │
Step 4 (Annotation) ◄─────────────── Requires Step 3 │
     │                                               │
     │ Requires final_annotation                     │
     ▼                                               │
Step 5 (Figures) ◄────────────────── Requires Step 4 │
     │                                               │
     │ Requires all pucks annotated                  │
     ▼                                               │
Step 6 (Summary) ◄────────────────── Requires Step 5 │
     │                                               │
     │ Requires all outputs                          │
     ▼                                               │
Step 7 (Verify) ◄─────────────────── Requires Step 6 ┘
```

**Parallel Execution Note:** 
Steps 2-5 can be parallelized ACROSS pucks (run Puck 29, 37, 40 simultaneously).
Steps 6-7 must wait for all pucks to complete.

---

## Success Criteria (Phase 1 Complete)

| Criterion | Target | Verification |
|-----------|--------|--------------|
| All pucks processed | 3/3 | Checkpoint file lists all |
| Cells per puck | >20,000 | Summary CSV |
| Mean confidence score | >0.5 | Summary CSV |
| High confidence cells | >70% | Summary CSV |
| Cell types detected | >5 unique | Summary CSV |
| Expected types present | T, B, myeloid, epithelial | Marker dotplot |
| Spatial coords preserved | Yes | adata.obsm['spatial'] |
| All figures generated | 15 total (5 per puck) | File check |
| Report generated | Yes | File exists |
| Checkpoint created | Yes | File exists |

---

## Notes

1. **DO NOT use R or Seurat** - This entire pipeline uses Python (popV, scanpy, squidpy)

2. **Process pucks individually** - Do not combine into all_pucks until Phase 2 if needed

3. **Raw counts required** - popV needs unnormalized data; check adata.raw if X is normalized

4. **Reference the example script** - `scripts/example_reference_scripts/SC_Cluster_Annotation.py` 
   contains the exact workflow adapted for scRNA-seq

5. **GPU optional** - popV can use GPU for faster inference but CPU works fine

6. **Cache the model** - Tabula Sapiens model is cached in `tmp/tabula_sapiens/` after first download

7. **Intermediate files** - `_popv_raw.h5ad` and `_clustered.h5ad` are intermediate; final output is `_celltype_annotated.h5ad`

8. **Checkpoint before Phase 2** - Verify `phase1_celltype_annotation.done` exists before starting Phase 2
