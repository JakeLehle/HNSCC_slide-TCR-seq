#!/usr/bin/env python3
"""
Step04_Cell_Type_Annotation.py

Automated cell type annotation and clustering for QC'd Slide-seq data.
Combined single-environment pipeline: popV + cell2location + Pearson
embedding, all in the unified sc_pre env.

Pipeline:
  1. Load merged QC h5ad (raw counts)
  2. popV annotation (Tabula Sapiens, on raw counts)        [checkpointed]
  3. Store raw counts in layers['counts'] (canonical source)
  4. cell2location deconvolution                            [checkpointed]
       4a. RegressionModel on HNSCC single-cell reference (cached signatures)
       4b. Cell2location spatial mapping (cached per detection_alpha)
       4c. per-bead argmax label (c2l_argmax) + abundances
  5. Pearson-residual embedding (from raw counts) -> PCA/UMAP/Leiden
       log1p-CPM written to X for plotting/DE ONLY (never enters the embedding)
  6. BBKNN batch correction (puck_id)
  7. Leiden clustering (corrected graph)
  8. popV cluster-consensus annotation (final_annotation)
  9. popV + cell2location reconciliation (consensus_annotation + agreement)
 10. Visualization (popV / c2l / consensus / agreement UMAPs + proportions)
 11. Save annotated h5ad (sanitized writes)

-----------------------------------------------------------------------------
CHANGELOG (v5)
-----------------------------------------------------------------------------
vs v4:
  - cell2location deconvolution folded in (unified sc_pre env now hosts both
    popV and cell2location, so no cross-env hand-off). Two-stage:
    RegressionModel signatures from the HNSCC single-cell reference, then
    Cell2location spatial mapping on the raw merged Slide-seq counts.
  - Both heavy stages are CHECKPOINTED: popV reuses the existing
    all_pucks_popv_annotated.h5ad if present; reference signatures cache to
    c2l_reference_signatures.csv; spatial abundances cache per detection_alpha
    (c2l_abundances_alpha{N}.csv). The alpha sweep reruns only stage 4b.
  - NEW consensus_annotation: popV (popv_prediction) reconciled with
    cell2location (c2l_argmax) via an editable LABEL_HARMONIZE map, plus a
    consensus_agreement flag. A popv-vs-c2l crosstab is written to disk to
    drive completing the map. final_annotation (popV cluster-consensus) is
    kept unchanged for downstream compatibility (SComatic/neoantigen).
  - Device-aware: accelerator = 'gpu' if torch.cuda.is_available() else 'cpu'.

  AUDIT (log-CPM vs Pearson, per Jake's question): the embedding is built
  ONLY from layers['counts'] (raw) via Pearson residuals -> X_pca. neighbors
  uses use_rep='X_pca', so X is never read by PCA/neighbors/UMAP/Leiden.
  log1p-CPM is applied to X AFTER X_pca is set and serves gene plots/DE only.
  See post_annotation_processing() for the explicit assertion + log line.

-----------------------------------------------------------------------------
CHANGELOG (v4)
-----------------------------------------------------------------------------
vs v3:
  - Embedding normalization switched from CPM+log1p to analytic Pearson
    residuals (sc.experimental.pp). The bake-off confirmed this collapses the
    library-size "snake" (depth_umap_corr 0.83 -> 0.16). Residuals drive the
    PCA/UMAP/Leiden ONLY; adata.X stays log-normalized for gene plotting and
    .raw. BBKNN, the figure set, and the popV/whitelist logic are unchanged.

-----------------------------------------------------------------------------
CHANGELOG (v3)
-----------------------------------------------------------------------------
vs v2:
  - popV merge now WHITELISTS popv_* (+ _predict_cells, over_clustering)
    instead of carrying all ~85 reference columns. Diagnostic
    (Step04_popv_write_diagnostic.py) confirmed this alone clears the
    "Can't implicitly convert non-string objects to strings" write error:
    the single offending column was observation_joinid (all-NaN object),
    which is not popv_-prefixed and is therefore dropped.
  - safe_write() / sanitize_obs_for_h5ad() retained as a cheap guard so a
    future popV column-set change cannot re-break the write.

vs original Step04 (ClusterCatcher parity / batch correction):
  - BBKNN batch correction (puck_id) restored. NOTE: Leiden runs AFTER
    BBKNN here so clusters reflect the corrected graph. ClusterCatcher runs
    Leiden BEFORE BBKNN; set CLUSTER_AFTER_BBKNN=False to match it exactly.
  - var_names restored from feature_name after popV (popV reindexes var to
    Ensembl). Guarded: no-op if no Ensembl IDs are present.
  - Labeled final_annotation UMAP (adjustText + path effects).
  - Cell-type proportions per puck (stacked bar, hex colors).
  - Before/after batch-correction UMAP (colored by puck).

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

# Batch-correction flags (kept out of the main import so spatial_config does
# not need editing; override there if you want them centralized).
try:
    from spatial_config import RUN_BBKNN, BBKNN_BATCH_KEY, CLUSTER_AFTER_BBKNN
except ImportError:
    RUN_BBKNN = True            # restore ClusterCatcher batch correction
    BBKNN_BATCH_KEY = 'puck_id'
    CLUSTER_AFTER_BBKNN = True  # cluster on the corrected graph (see CHANGELOG)

# Pearson-residual embedding HVG count (validated in the normalization bake-off)
N_PEARSON_HVG = 2000

# ---- cell2location configuration ----------------------------------------
# Reference = NMF-paper HNSCC single-cell object built earlier
# (build_cell2location_reference.py): raw counts in X + layers['counts'],
# obs['final_annotation'] (12 types), batch column 'subject id'.
try:
    from spatial_config import C2L_REFERENCE_PATH
except ImportError:
    C2L_REFERENCE_PATH = (
        "/work/sdz852/WORKING/slide-TCR-seq-working/data/inputs/annotation/"
        "reference_for_cell2location.h5ad"
    )
C2L_REF_BATCH_KEY = "subject_id"        # 'subject id' is renamed to this
C2L_LABELS_KEY    = "final_annotation"  # reference cell-type column

C2L_CFG = {
    "n_cells_per_location": 1,   # near-single-cell Slide-seq beads
    "detection_alpha":      20,  # high within-slide variability; try 200 too
    "ref_epochs":           250,
    "map_epochs":           30000,
}

# popV (Tabula Sapiens vocab) <-> cell2location (HNSCC ref vocab) harmonization.
# Maps BOTH vocabularies onto a shared coarse label so the two annotators can
# be compared. The cell2location side is seeded from the reference's 12 types.
# The popV side is INTENTIONALLY INCOMPLETE: review the popv-vs-c2l crosstab
# written by build_consensus_annotation(), then fill in the popV strings.
LABEL_HARMONIZE = {
    # ---- cell2location reference labels (final_annotation) -> shared ----
    "basal cell":   "epithelial",
    "CD4 T":        "CD4_T",
    "CD8 T":        "CD8_T",
    "Treg":         "Treg",
    "B":            "B",
    "macrophage":   "macrophage",
    "fibroblast":   "fibroblast",
    "endothelial":  "endothelial",
    "smooth muscle":"smooth_muscle",
    "pDC":          "pDC",
    "mast":         "mast",
    "mDC":          "mDC",
    # ---- popV labels (Tabula Sapiens) -> shared  [FILL IN FROM CROSSTAB] ----
    # e.g. "cd4-positive, alpha-beta t cell": "CD4_T",
    #      "cd8-positive, alpha-beta t cell": "CD8_T",
    #      "b cell": "B", "macrophage": "macrophage", ...
}

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import to_hex


# =========================================================================
# H5AD WRITE SAFETY
# =========================================================================

def sanitize_obs_for_h5ad(adata):
    """
    Coerce obs/var into h5ad-writable dtypes before write_h5ad().

    With the popv_* whitelist this is mostly belt-and-suspenders, but it
    guards against future popV column drift:
      - var index (and raw.var index) -> str
      - object columns -> fillna('').astype(str)   [normalizes numpy.str_]
      - EMPTY categoricals (zero categories) -> filled str

    LEFT UNTOUCHED (so they stay usable downstream):
      - float / int numeric columns
      - int-valued categoricals (popv_prediction_score, popv_prediction_depth)
    """
    adata.var.index = adata.var.index.astype(str)
    if adata.raw is not None:
        try:
            adata.raw.var.index = adata.raw.var.index.astype(str)
        except Exception as e:
            log(f"    note: could not recast raw.var index ({e}); continuing")

    n_obj, n_empty_cat = 0, 0
    for col in adata.obs.columns:
        s = adata.obs[col]
        if s.dtype == object:
            adata.obs[col] = s.fillna('').astype(str)
            n_obj += 1
        elif isinstance(s.dtype, pd.CategoricalDtype) and len(s.cat.categories) == 0:
            filled = s.astype(object).where(s.notna(), '')
            adata.obs[col] = filled.astype(str)
            n_empty_cat += 1

    log(f"    sanitized obs: {n_obj} object cols -> str, "
        f"{n_empty_cat} empty categoricals -> str")
    return adata


def safe_write(adata, path):
    """Sanitize obs/var, then write h5ad. Use for EVERY write in this script."""
    sanitize_obs_for_h5ad(adata)
    adata.write_h5ad(path)
    return path


# =========================================================================
# POPV ANNOTATION
# =========================================================================

def run_popv_annotation(adata, cache_dir):
    """
    Run popV cell type annotation on raw counts.

    Steps:
      1. Save original obs columns
      2. Ensure feature_name (gene symbol) var column exists
      3. Pull Tabula Sapiens HubModel
      4. Annotate (raw counts)
      5. Restore var_names from feature_name (popV reindexes to Ensembl)
      6. Filter to query cells
      7. Merge popV columns back -- WHITELIST popv_* only
    """
    banner("POPV CELL TYPE ANNOTATION")

    huggingface_repo = "popV/tabula_sapiens_All_Cells"
    query_batch_key = "puck_id"
    gene_symbol_key = "feature_name"

    log(f"  Model:     {huggingface_repo}")
    log(f"  Batch key: {query_batch_key}")
    log(f"  Cells:     {adata.n_obs:,}")
    log(f"  Genes:     {adata.n_vars:,}")

    if gene_symbol_key not in adata.var.columns:
        adata.var[gene_symbol_key] = adata.var_names
        log(f"  Created var['{gene_symbol_key}'] from var_names")

    original_obs_columns = list(adata.obs.columns)

    log(f"  Pulling model from HuggingFace...")
    ensure_dir(cache_dir)
    hmo = popv.hub.HubModel.pull_from_huggingface_hub(
        huggingface_repo, cache_dir=cache_dir
    )

    log(f"  Running popV annotation (this may take 30-60 minutes)...")
    adata_annotated = hmo.annotate_data(
        adata,
        query_batch_key=query_batch_key,
        prediction_mode="inference",
        gene_symbols=gene_symbol_key,
    )

    # ---------------------------------------------------------------------
    # Restore var_names from feature_name (ClusterCatcher parity).
    # popV's annotate_data() re-indexes var_names to Ensembl IDs in place.
    # Downstream (SComatic / neoantigen) needs HUGO symbols. Guarded so this
    # is a no-op if no Ensembl IDs are present.
    # ---------------------------------------------------------------------
    if gene_symbol_key in adata.var.columns:
        n_ensembl = adata.var_names.to_series().astype(str).str.startswith('ENSG').sum()
        if n_ensembl > 0:
            log(f"  Restoring var_names from '{gene_symbol_key}' "
                f"({n_ensembl} Ensembl IDs detected)...")
            adata.var.index = adata.var[gene_symbol_key].astype(str).values
            adata.var_names_make_unique()

    # Filter to query cells only
    if '_dataset' in adata_annotated.obs.columns:
        log(f"  Filtering to query cells only...")
        adata_annotated = adata_annotated[
            adata_annotated.obs['_dataset'] == 'query'
        ].copy()

    # -- WHITELIST merge: keep popv_* (+ two popV internals), drop the rest --
    # Confirmed by Step04_popv_write_diagnostic.py: this removes the
    # all-NaN object column (observation_joinid) and the empty Tabula Sapiens
    # reference columns, leaving only populated, writable popV outputs. All
    # per-method prediction columns are retained for model-bias inspection.
    log(f"  Merging popV annotations (whitelist: popv_* + internals)...")
    popv_columns = [c for c in adata_annotated.obs.columns
                    if c.startswith('popv_') or c in ['_predict_cells', 'over_clustering']]
    popv_obs = adata_annotated.obs[popv_columns]
    log(f"    kept {len(popv_columns)} popV columns "
        f"(dropped {adata_annotated.obs.shape[1] - len(original_obs_columns) - len(popv_columns)} "
        f"reference/internal columns)")

    merged_df = pd.merge(
        adata.obs, popv_obs,
        left_index=True, right_index=True, how='inner'
    )

    n_dropped = adata.n_obs - len(merged_df)
    adata = adata[adata.obs_names.isin(merged_df.index)].copy()
    merged_df = merged_df.loc[adata.obs_names]
    adata.obs = merged_df

    assert all(adata.obs_names == merged_df.index), "Index mismatch after popV merge"

    if n_dropped:
        log(f"  NOTE: {n_dropped} cells dropped (popV excluded for low expression)")
    log(f"  Merged {len(popv_obs.columns)} popV columns, {adata.n_obs:,} cells retained")

    if 'popv_prediction' in adata.obs.columns:
        n_types = adata.obs['popv_prediction'].nunique()
        log(f"  Cell types detected: {n_types}")
        log(f"  Top 10:")
        for ct, count in adata.obs['popv_prediction'].value_counts().head(10).items():
            pct = 100 * count / adata.n_obs
            log(f"    {ct}: {count:,} ({pct:.1f}%)")

    return adata


# =========================================================================
# CELL2LOCATION DECONVOLUTION
# =========================================================================

def _c2l_accelerator():
    """Return 'gpu' if a CUDA device is visible, else 'cpu' (with a warning)."""
    import torch
    if torch.cuda.is_available():
        return "gpu"
    log("  WARNING: no CUDA device visible to torch. cell2location training "
        "on CPU is extremely slow (hours-to-days); resolve the GPU first.")
    return "cpu"


def get_reference_signatures(out_dir, accelerator):
    """
    Stage 4a: per-cell-type expression signatures from the HNSCC single-cell
    reference (cached to CSV; trained only once).
    """
    sig_path = os.path.join(out_dir, "c2l_reference_signatures.csv")
    if os.path.exists(sig_path):
        log(f"  [cache] reference signatures: {sig_path}")
        return pd.read_csv(sig_path, index_col=0)

    from cell2location.models import RegressionModel
    from cell2location.utils.filtering import filter_genes

    log(f"  Loading reference: {C2L_REFERENCE_PATH}")
    ref = sc.read_h5ad(C2L_REFERENCE_PATH)

    # 'subject id' (with a space) -> 'subject_id'
    if 'subject id' in ref.obs.columns and C2L_REF_BATCH_KEY not in ref.obs.columns:
        ref.obs = ref.obs.rename(columns={'subject id': C2L_REF_BATCH_KEY})
    # raw integer counts in X
    if 'counts' in ref.layers:
        ref.X = ref.layers['counts'].copy()
    assert C2L_LABELS_KEY in ref.obs.columns, \
        f"reference missing labels column '{C2L_LABELS_KEY}'"
    assert C2L_REF_BATCH_KEY in ref.obs.columns, \
        f"reference missing batch column '{C2L_REF_BATCH_KEY}'"

    # cell2location-recommended permissive gene filter
    selected = filter_genes(
        ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )
    ref = ref[:, selected].copy()
    log(f"    reference after gene filter: {ref.n_obs:,} x {ref.n_vars:,}")

    RegressionModel.setup_anndata(
        adata=ref, batch_key=C2L_REF_BATCH_KEY, labels_key=C2L_LABELS_KEY,
    )
    mod = RegressionModel(ref)
    log(f"    training RegressionModel (max_epochs={C2L_CFG['ref_epochs']}, "
        f"accelerator={accelerator})...")
    # NOTE (version-sensitive): scvi-tools 1.4 uses accelerator=; older used
    # use_gpu=. If this errors on the device kwarg, that is the line to adjust.
    mod.train(max_epochs=C2L_CFG['ref_epochs'], accelerator=accelerator)
    ref = mod.export_posterior(
        ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500},
    )

    # canonical signature extraction (cell2location tutorial)
    factor_names = ref.uns['mod']['factor_names']
    cols = [f'means_per_cluster_mu_fg_{i}' for i in factor_names]
    if 'means_per_cluster_mu_fg' in ref.varm.keys():
        inf_aver = ref.varm['means_per_cluster_mu_fg'][cols].copy()
    else:
        inf_aver = ref.var[cols].copy()
    inf_aver.columns = factor_names
    inf_aver.to_csv(sig_path)
    log(f"    saved signatures: {sig_path} "
        f"({inf_aver.shape[0]} genes x {inf_aver.shape[1]} cell types)")
    return inf_aver


def run_cell2location(adata, out_dir):
    """
    Stage 4: cell2location deconvolution of the spatial beads.

    Reads RAW counts from adata.layers['counts'] (set in main, before any
    normalization). Writes per-bead abundances + argmax label to adata.obs.
    The spatial mapping is cached per detection_alpha so the hyperparameter
    sweep does not retrain the reference or rerun popV.
    """
    banner("CELL2LOCATION DECONVOLUTION")
    import cell2location
    from cell2location.models import Cell2location

    accelerator = _c2l_accelerator()
    inf_aver = get_reference_signatures(out_dir, accelerator)

    alpha = C2L_CFG['detection_alpha']
    abund_path = os.path.join(out_dir, f"c2l_abundances_alpha{alpha}.csv")

    if os.path.exists(abund_path):
        log(f"  [cache] spatial abundances: {abund_path}")
        abund = pd.read_csv(abund_path, index_col=0)
    else:
        assert 'counts' in adata.layers, \
            "layers['counts'] (raw) required before cell2location"
        sp = ad.AnnData(
            X=adata.layers['counts'].copy(),
            obs=adata.obs.copy(),
            var=adata.var.copy(),
        )
        # genes shared between spatial and reference signatures
        shared = [g for g in sp.var_names if g in set(inf_aver.index)]
        log(f"  genes shared with reference signatures: {len(shared):,}")
        if len(shared) < 200:
            raise RuntimeError(
                f"Only {len(shared)} genes shared between spatial var_names and "
                f"the reference signatures. This usually means the popV "
                f"checkpoint kept Ensembl IDs while the reference uses HUGO "
                f"symbols. Confirm adata.var_names are symbols (e.g. 'CD14'), "
                f"not 'ENSG...'; if Ensembl, delete the popV checkpoint to "
                f"force a fresh run with var_name restoration, or remap here."
            )
        sp = sp[:, shared].copy()
        inf_use = inf_aver.loc[shared].copy()

        Cell2location.setup_anndata(adata=sp, batch_key='puck_id')
        mod = Cell2location(
            sp, cell_state_df=inf_use,
            N_cells_per_location=C2L_CFG['n_cells_per_location'],
            detection_alpha=alpha,
        )
        log(f"  training Cell2location (alpha={alpha}, "
            f"N_cells={C2L_CFG['n_cells_per_location']}, "
            f"max_epochs={C2L_CFG['map_epochs']}, accelerator={accelerator})...")
        mod.train(max_epochs=C2L_CFG['map_epochs'],
                  batch_size=None, train_size=1, accelerator=accelerator)
        sp = mod.export_posterior(
            sp, sample_kwargs={'num_samples': 1000, 'batch_size': sp.n_obs},
        )
        abund = sp.obsm['q05_cell_abundance_w_sf'].copy()
        abund.columns = [c.replace('q05cell_abundance_w_sf_', '')
                         for c in abund.columns]
        abund.to_csv(abund_path)
        log(f"  saved abundances: {abund_path}")

    # align abundances to the working object and write obs columns
    abund = abund.reindex(adata.obs_names)
    for ct in abund.columns:
        adata.obs[f'c2l_{ct}'] = pd.to_numeric(abund[ct], errors='coerce').values
    adata.obs['c2l_argmax'] = abund.astype(float).idxmax(axis=1).values
    adata.obsm['c2l_abundance'] = abund.astype(float).values
    adata.uns['c2l_celltypes'] = list(abund.columns)
    adata.uns['c2l_detection_alpha'] = alpha

    log(f"  cell2location argmax cell types: {adata.obs['c2l_argmax'].nunique()}")
    for ct, n in adata.obs['c2l_argmax'].value_counts().head(10).items():
        log(f"    {ct}: {n:,} ({100 * n / adata.n_obs:.1f}%)")
    return adata


# =========================================================================
# POST-ANNOTATION PROCESSING
# =========================================================================

def post_annotation_processing(adata):
    """
    Normalize, embed, batch-correct, and cluster after popV annotation.

    Order:
      normalize -> HVG (batch-aware) -> PCA on HVGs -> neighbors -> UMAP
      -> [BBKNN if RUN_BBKNN] -> recompute UMAP -> Leiden

    The pre-correction UMAP is stashed in obsm['X_umap_uncorrected'] so the
    before/after batch figure can be drawn.
    """
    banner("POST-ANNOTATION PROCESSING")

    # Raw counts are the canonical source for both the Pearson embedding and
    # cell2location. main() stores them before this runs; only set here if a
    # caller invoked this function standalone without a counts layer.
    if 'counts' not in adata.layers:
        log(f"  Storing raw counts in layers['counts']...")
        adata.layers['counts'] = adata.X.copy()
    else:
        log(f"  Using existing layers['counts'] (raw) as embedding source")

    # --- Embedding via analytic Pearson residuals -----------------------
    # The normalization bake-off (Step04_normalization_bakeoff.py) showed
    # CPM+log1p left a library-size "snake" (depth_umap_corr 0.83) while
    # Pearson residuals removed it (0.16). Residuals are computed from
    # layers['counts'] on a throwaway copy; only X_pca is copied back. X is
    # NOT read here, so the embedding cannot depend on whatever is in X.
    log(f"  Building Pearson-residual embedding (HVG={N_PEARSON_HVG})...")
    emb = adata.copy()
    emb.X = emb.layers['counts'].copy()
    sc.experimental.pp.highly_variable_genes(
        emb, flavor='pearson_residuals', n_top_genes=N_PEARSON_HVG
    )
    emb = emb[:, emb.var['highly_variable']].copy()
    sc.experimental.pp.normalize_pearson_residuals(emb)
    n_comps = min(50, emb.n_vars - 1, emb.n_obs - 1)
    sc.tl.pca(emb, svd_solver='arpack', n_comps=n_comps)
    adata.obsm['X_pca'] = emb.obsm['X_pca']
    adata.uns['pca'] = emb.uns['pca']
    adata.var['highly_variable'] = adata.var_names.isin(emb.var_names)
    log(f"    Pearson HVGs: {emb.n_vars:,}, PCA components: {n_comps}")
    del emb

    # AUDIT: X_pca exists and is Pearson-derived BEFORE log-CPM touches X.
    # Everything downstream (neighbors/UMAP/Leiden) reads use_rep='X_pca'.
    assert 'X_pca' in adata.obsm, "embedding (X_pca) not set before log-CPM"
    log(f"  AUDIT: embedding = Pearson(counts) -> X_pca; log-CPM below is "
        f"for plotting/DE only and never enters PCA/neighbors/UMAP/Leiden")

    # Log-normalized X for visualization / gene plots / .raw (NOT the embedding)
    log(f"  Log-normalizing X for plotting (CPM + log1p)...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    n_pcs = min(N_PCS, n_comps)

    # Pre-correction neighbors + UMAP (stash for before/after figure)
    log(f"  Computing neighbors (pre-correction, n_pcs={n_pcs})...")
    sc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca')
    log(f"  Computing UMAP (pre-correction)...")
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    adata.obsm['X_umap_uncorrected'] = adata.obsm['X_umap'].copy()

    # BBKNN batch correction
    corrected = False
    if RUN_BBKNN:
        log(f"  BBKNN batch correction (batch_key={BBKNN_BATCH_KEY})...")
        try:
            import bbknn  # noqa: F401  (import check; sc.external.pp.bbknn uses it)
            sc.external.pp.bbknn(adata, batch_key=BBKNN_BATCH_KEY, n_pcs=n_pcs)
            log(f"  Recomputing UMAP (batch-corrected)...")
            sc.tl.umap(adata, random_state=RANDOM_SEED)
            corrected = True
            log(f"    BBKNN done; neighbor graph + UMAP are batch-corrected")
        except Exception as e:
            log(f"    BBKNN unavailable/failed ({e}); continuing UNCORRECTED")
            log(f"    (install with: pip install bbknn)")
    else:
        log(f"  RUN_BBKNN=False; skipping batch correction")

    # Leiden clustering
    if corrected and not CLUSTER_AFTER_BBKNN:
        # ClusterCatcher-exact behavior: would have clustered before BBKNN.
        # We deliberately do not support that path cleanly here because the
        # pre-correction graph has been overwritten; flip CLUSTER_AFTER_BBKNN
        # only if you intend to cluster on uncorrected neighbors upstream.
        log(f"  WARNING: CLUSTER_AFTER_BBKNN=False but graph already corrected; "
            f"clustering on corrected graph anyway")
    graph_state = "corrected" if corrected else "uncorrected"
    log(f"  Leiden clustering on {graph_state} graph (resolution={LEIDEN_RESOLUTION})...")
    sc.tl.leiden(adata, key_added='clusters',
                 resolution=LEIDEN_RESOLUTION, random_state=RANDOM_SEED)
    n_clusters = adata.obs['clusters'].nunique()
    log(f"    Found {n_clusters} clusters")

    adata.uns['batch_corrected'] = corrected
    return adata


# =========================================================================
# CONSENSUS ANNOTATION (ClusterCatcher weighted scoring)
# =========================================================================

def assign_cluster_based_annotation(adata):
    """
    Assign final cell type labels using cluster-level weighted scoring.

    Follows ClusterCatcher's assign_cluster_based_annotation():
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
# popV + cell2location RECONCILIATION
# =========================================================================

def build_consensus_annotation(adata, ann_dir):
    """
    Reconcile popV (popv_prediction) with cell2location (c2l_argmax).

    Because the two annotators use different label vocabularies, both are
    mapped through LABEL_HARMONIZE onto a shared coarse vocabulary before
    comparison. Unmapped labels are passed through verbatim (so they will
    simply never 'agree' until added to the map).

    Writes:
      obs['consensus_annotation'] : shared label when both agree, else
                                    'ambiguous:<popV>|<c2l>'
      obs['consensus_agreement']  : 'agree' / 'disagree'
      <ann_dir>/popv_vs_c2l_crosstab.tsv : empirical correspondence table
                                           (use this to complete LABEL_HARMONIZE)
    """
    banner("popV + cell2location RECONCILIATION")

    if 'c2l_argmax' not in adata.obs.columns:
        log(f"  c2l_argmax missing; skipping reconciliation "
            f"(cell2location stage did not run)")
        return adata

    pv = adata.obs['popv_prediction'].astype(str)
    c2 = adata.obs['c2l_argmax'].astype(str)

    # Raw crosstab BEFORE harmonization: this is the diagnostic that tells you
    # which popV string corresponds to which cell2location type empirically.
    ctab = pd.crosstab(pv, c2)
    ctab_path = os.path.join(ann_dir, "popv_vs_c2l_crosstab.tsv")
    ctab.to_csv(ctab_path, sep='\t')
    log(f"  wrote popV-vs-c2l crosstab: {ctab_path}")
    log(f"    (review this, then complete the popV side of LABEL_HARMONIZE)")

    pv_h = pv.map(lambda x: LABEL_HARMONIZE.get(x, x))
    c2_h = c2.map(lambda x: LABEL_HARMONIZE.get(x, x))
    agree = (pv_h.values == c2_h.values)

    adata.obs['consensus_agreement'] = np.where(agree, 'agree', 'disagree')
    adata.obs['consensus_annotation'] = np.where(
        agree, pv_h.values,
        np.char.add(np.char.add('ambiguous:' + pv_h.values.astype(str), '|'),
                    c2_h.values.astype(str))
    )

    n_mapped_pv = pv.isin(LABEL_HARMONIZE).sum()
    log(f"  popV labels covered by LABEL_HARMONIZE: "
        f"{n_mapped_pv:,}/{adata.n_obs:,} beads "
        f"({100 * n_mapped_pv / adata.n_obs:.1f}%)")
    log(f"  popV/cell2location agreement: "
        f"{100 * agree.mean():.1f}% of beads "
        f"(meaningful only once the popV side of the map is filled in)")
    return adata


# =========================================================================
# VISUALIZATION
# =========================================================================

def gen_mpl_labels(adata, groupby, exclude=(), ax=None,
                   adjust_kwargs=None, text_kwargs=None):
    """Non-overlapping labels at category medians (ClusterCatcher pattern)."""
    try:
        from adjustText import adjust_text
    except ImportError:
        log(f"    adjustText not installed; skipping label placement")
        return

    adjust_kwargs = adjust_kwargs or {}
    text_kwargs = text_kwargs or {}

    medians = {}
    for g, idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[idx].obsm["X_umap"], axis=0)

    if ax is None:
        texts = [plt.text(x, y, k, **text_kwargs) for k, (x, y) in medians.items()]
    else:
        texts = [ax.text(x, y, k, **text_kwargs) for k, (x, y) in medians.items()]
    adjust_text(texts, **adjust_kwargs)


def generate_annotation_plots(adata, fig_dir):
    """Generate the full ClusterCatcher-parity annotation plot set + spatial."""
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

    # --- UMAP: final annotation (plain legend) ---
    if 'final_annotation' in adata.obs.columns:
        log(f"  UMAP: final_annotation")
        adata.obs['final_annotation'] = adata.obs['final_annotation'].astype('category')
        fig, ax = plt.subplots(figsize=(14, 10))
        sc.pl.umap(adata, color='final_annotation', ax=ax, show=False,
                   legend_loc='right margin', frameon=False, size=3)
        plt.tight_layout()
        save_fig(fig, 'UMAP_final_annotation', fig_dir)

        # --- UMAP: final annotation with non-overlapping labels ---
        log(f"  UMAP: final_annotation (labeled)")
        plot_umap_final_annotation_labeled(adata, fig_dir)

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

    # --- Before/after batch correction (puck) ---
    if 'X_umap_uncorrected' in adata.obsm:
        log(f"  UMAP: batch correction before/after")
        plot_batch_correction_comparison(adata, fig_dir)

    # --- Stacked bar: cell type composition per cluster ---
    if 'popv_prediction' in adata.obs.columns:
        log(f"  Stacked bar: cluster composition")
        plot_stacked_bar(adata, fig_dir)

    # --- Cell type proportions per puck ---
    if 'final_annotation' in adata.obs.columns:
        log(f"  Stacked bar: cell type proportions per puck")
        plot_cell_type_proportions(adata, fig_dir)

    # --- Spatial plots per puck ---
    annotation_col = 'final_annotation' if 'final_annotation' in adata.obs.columns else 'clusters'
    for puck_name in PUCK_NAMES:
        mask = adata.obs['puck_id'] == puck_name
        if mask.sum() == 0:
            continue
        log(f"  Spatial: {puck_name}")
        plot_spatial_scatter(adata[mask].copy(), puck_name, annotation_col, fig_dir)
        plot_spatial_scatter(adata[mask].copy(), puck_name, 'clusters', fig_dir)


def plot_umap_final_annotation_labeled(adata, fig_dir):
    """final_annotation UMAP with adjustText labels + white-stroke path effects."""
    adata.obs['final_annotation'] = adata.obs['final_annotation'].astype('category')
    effects = [
        pe.withStroke(linewidth=6, foreground="white"),
        pe.Normal(),
    ]
    fig, ax = plt.subplots(figsize=(14, 12))
    sc.pl.umap(adata, color='final_annotation', ax=ax, show=False,
               legend_loc=None, frameon=False, size=3)
    gen_mpl_labels(
        adata, 'final_annotation',
        exclude=("", "None", "nan", "Unknown"),
        ax=ax,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black', lw=0.5)),
        text_kwargs=dict(fontsize=16, fontweight='bold', path_effects=effects),
    )
    fig.tight_layout()
    save_fig(fig, 'UMAP_final_annotation_labeled', fig_dir)


def plot_consensus_umaps(adata, fig_dir):
    """UMAPs for cell2location argmax, the popV/c2l consensus, and agreement."""
    effects = [pe.withStroke(linewidth=6, foreground="white"), pe.Normal()]

    # cell2location argmax (categorical, labeled)
    if 'c2l_argmax' in adata.obs.columns:
        adata.obs['c2l_argmax'] = adata.obs['c2l_argmax'].astype('category')
        fig, ax = plt.subplots(figsize=(14, 12))
        sc.pl.umap(adata, color='c2l_argmax', ax=ax, show=False,
                   legend_loc=None, frameon=False, size=3)
        gen_mpl_labels(
            adata, 'c2l_argmax', exclude=("", "None", "nan"), ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black', lw=0.5)),
            text_kwargs=dict(fontsize=16, fontweight='bold', path_effects=effects),
        )
        fig.tight_layout()
        save_fig(fig, 'UMAP_c2l_argmax_labeled', fig_dir)

    # consensus annotation (categorical, labeled)
    if 'consensus_annotation' in adata.obs.columns:
        adata.obs['consensus_annotation'] = adata.obs['consensus_annotation'].astype('category')
        fig, ax = plt.subplots(figsize=(14, 12))
        sc.pl.umap(adata, color='consensus_annotation', ax=ax, show=False,
                   legend_loc=None, frameon=False, size=3)
        gen_mpl_labels(
            adata, 'consensus_annotation',
            exclude=("", "None", "nan"), ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black', lw=0.5)),
            text_kwargs=dict(fontsize=14, fontweight='bold', path_effects=effects),
        )
        fig.tight_layout()
        save_fig(fig, 'UMAP_consensus_annotation_labeled', fig_dir)

    # agreement map (two fixed hex colors)
    if 'consensus_agreement' in adata.obs.columns:
        adata.obs['consensus_agreement'] = adata.obs['consensus_agreement'].astype('category')
        fig, ax = plt.subplots(figsize=(12, 11))
        sc.pl.umap(
            adata, color='consensus_agreement', ax=ax, show=False,
            frameon=False, size=3,
            palette={'agree': '#4c956c', 'disagree': '#d1495b'},
        )
        fig.tight_layout()
        save_fig(fig, 'UMAP_consensus_agreement', fig_dir)



    pucks = sorted(adata.obs['puck_id'].astype(str).unique())
    cmap = plt.colormaps['turbo']
    colors = {p: to_hex(cmap(i / max(len(pucks) - 1, 1)))
              for i, p in enumerate(pucks)}

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))
    panels = [
        (axes[0], 'X_umap_uncorrected', 'Before (uncorrected)'),
        (axes[1], 'X_umap', 'After BBKNN' if adata.uns.get('batch_corrected') else 'After (no correction)'),
    ]
    for ax, key, title in panels:
        emb = adata.obsm[key]
        for p in pucks:
            m = (adata.obs['puck_id'].astype(str) == p).values
            ax.scatter(emb[m, 0], emb[m, 1], s=2, c=colors[p],
                       label=p, alpha=0.6, rasterized=True)
        ax.set_title(title, fontsize=FONT_SIZE - 6)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                   markerscale=6, fontsize=FONT_SIZE - 16)
    plt.tight_layout()
    save_fig(fig, 'UMAP_batch_correction_puck', fig_dir)


def plot_stacked_bar(adata, fig_dir):
    """Stacked bar of cell type composition per Leiden cluster (hex colors)."""
    cluster_counts = adata.obs['clusters'].value_counts()
    clusters = sorted(
        [c for c in cluster_counts.index if cluster_counts[c] >= 50],
        key=lambda x: int(x) if str(x).isdigit() else x
    )

    fig, ax = plt.subplots(figsize=(max(16, len(clusters) * 1.2), 7))

    all_types = adata.obs['popv_prediction'].astype('category').cat.categories
    cmap = plt.colormaps['turbo']
    colors = {ct: to_hex(cmap(i / max(len(all_types) - 1, 1)))
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


def plot_cell_type_proportions(adata, fig_dir):
    """Stacked bar of final cell type proportions per puck (hex colors)."""
    props = pd.crosstab(
        adata.obs['puck_id'], adata.obs['final_annotation'], normalize='index'
    ) * 100
    all_types = list(props.columns)
    cmap = plt.colormaps['turbo']
    colors = [to_hex(cmap(i / max(len(all_types) - 1, 1)))
              for i in range(len(all_types))]

    fig, ax = plt.subplots(figsize=(12, 8))
    props.plot(kind='bar', stacked=True, ax=ax, color=colors,
               width=0.7, edgecolor='white', linewidth=0.3)
    ax.set_ylabel('Percentage')
    ax.set_xlabel('Puck')
    ax.set_title('Cell Type Proportions per Puck')
    ax.set_ylim(0, 100)
    ax.legend(title='Cell type', bbox_to_anchor=(1.02, 1), loc='upper left',
              fontsize=FONT_SIZE - 16)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    save_fig(fig, 'CellType_Proportions_per_Puck', fig_dir)


def plot_spatial_scatter(adata_puck, puck_name, color_by, fig_dir):
    """Spatial scatter plot for a single puck colored by annotation."""
    categories = sorted(adata_puck.obs[color_by].astype(str).unique())
    cmap = plt.colormaps['turbo']
    colors = {cat: cmap(i / max(len(categories) - 1, 1))
              for i, cat in enumerate(categories)}

    fig, ax = plt.subplots(figsize=(12, 10))
    for cat in categories:
        m = adata_puck.obs[color_by].astype(str) == cat
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
    """Export annotation and cluster-mapping summary tables."""
    banner("EXPORTING SUMMARIES")

    if 'final_annotation' in adata.obs.columns:
        ct_summary = pd.crosstab(
            adata.obs['puck_id'], adata.obs['final_annotation']
        )
        ct_summary['total_cells'] = ct_summary.sum(axis=1)
        ct_path = os.path.join(ann_dir, "annotation_summary.tsv")
        ct_summary.to_csv(ct_path, sep='\t')
        log(f"  Saved: {ct_path}")

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
    banner("STEP 04: CELL TYPE ANNOTATION + CELL2LOCATION (v5)")
    log(f"  Input:  {DIR_03_QC}")
    log(f"  Output: {DIR_04_ANNOTATION}")
    log(f"  BBKNN:  RUN_BBKNN={RUN_BBKNN}, batch_key={BBKNN_BATCH_KEY}, "
        f"cluster_after_bbknn={CLUSTER_AFTER_BBKNN}")
    log(f"  c2l:    alpha={C2L_CFG['detection_alpha']}, "
        f"N_cells={C2L_CFG['n_cells_per_location']}, ref={C2L_REFERENCE_PATH}")

    ann_dir = ensure_dir(DIR_04_ANNOTATION)
    ann_fig_dir = ensure_dir(os.path.join(DIR_04_ANNOTATION, "figures"))
    cache_dir = os.path.join(ann_dir, "popv_cache")
    c2l_dir = ensure_dir(os.path.join(ann_dir, "cell2location"))

    # --- popV annotation (reuse checkpoint if present, else run + save) ---
    popv_path = os.path.join(ann_dir, "all_pucks_popv_annotated.h5ad")
    if os.path.exists(popv_path):
        banner("LOADING popV CHECKPOINT")
        adata = sc.read_h5ad(popv_path)
        log(f"  [cache] {popv_path}")
        log(f"  Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        if 'popv_prediction' not in adata.obs.columns:
            log(f"  FATAL: checkpoint missing popv_prediction; "
                f"delete it to force a fresh popV run.")
            sys.exit(1)
    else:
        banner("LOADING MERGED QC DATA")
        merged_path = os.path.join(DIR_03_QC, "all_pucks_merged_QC.h5ad")
        if not os.path.exists(merged_path):
            log(f"  FATAL: {merged_path} not found. Run Step03 first.")
            sys.exit(1)
        adata = sc.read_h5ad(merged_path)
        log(f"  Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        log(f"  Pucks: {adata.obs['puck_id'].value_counts().to_dict()}")

        adata = run_popv_annotation(adata, cache_dir)
        safe_write(adata, popv_path)
        log(f"  Checkpoint saved: {popv_path}")

    # --- Canonical raw counts (source for BOTH cell2location and Pearson) ---
    # At this point adata.X is still raw integer counts (popV does not modify
    # X in place for the query). Stash before any normalization touches X.
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()
        log(f"  Stored raw counts in layers['counts'] "
            f"(canonical source for c2l + Pearson)")

    # --- cell2location deconvolution (checkpointed) ---
    adata = run_cell2location(adata, c2l_dir)

    # --- Post-annotation processing (Pearson embedding, BBKNN, cluster) ---
    adata = post_annotation_processing(adata)

    # --- popV cluster-based annotation (final_annotation, unchanged) ---
    adata = assign_cluster_based_annotation(adata)

    # --- popV + cell2location reconciliation (consensus_annotation) ---
    adata = build_consensus_annotation(adata, ann_dir)

    # --- Visualization ---
    generate_annotation_plots(adata, ann_fig_dir)
    plot_consensus_umaps(adata, ann_fig_dir)

    # --- Export summaries ---
    export_summaries(adata, ann_dir)

    # --- Save final annotated data (sanitized write) ---
    banner("SAVING FINAL ANNOTATED DATA")
    final_path = os.path.join(ann_dir, "all_pucks_annotated.h5ad")
    safe_write(adata, final_path)
    size_mb = os.path.getsize(final_path) / 1e6
    log(f"  Saved: {final_path} ({size_mb:.1f} MB)")

    # Per-puck annotated h5ad (sanitized writes)
    for puck_name in PUCK_NAMES:
        mask = adata.obs['puck_id'] == puck_name
        if mask.sum() == 0:
            continue
        puck_path = os.path.join(ann_dir, f"{puck_name}_annotated.h5ad")
        safe_write(adata[mask].copy(), puck_path)
        log(f"  Saved: {puck_path}")

    # --- Final summary ---
    banner("STEP 04 SUMMARY")
    log(f"  Total cells:   {adata.n_obs:,}")
    log(f"  Total genes:   {adata.n_vars:,}")
    log(f"  Clusters:      {adata.obs['clusters'].nunique()}")
    log(f"  Batch corrected: {adata.uns.get('batch_corrected', False)}")
    if 'final_annotation' in adata.obs.columns:
        log(f"  Cell types:    {adata.obs['final_annotation'].nunique()}")
    if 'c2l_argmax' in adata.obs.columns:
        log(f"  c2l argmax types: {adata.obs['c2l_argmax'].nunique()} "
            f"(detection_alpha={adata.uns.get('c2l_detection_alpha', '?')})")
    if 'consensus_agreement' in adata.obs.columns:
        frac = (adata.obs['consensus_agreement'] == 'agree').mean()
        log(f"  popV/c2l agreement: {100 * frac:.1f}% "
            f"(after LABEL_HARMONIZE; complete the map for a true number)")
    log(f"  Layers:        {list(adata.layers.keys())}")

    log(f"\n  Cells per puck:")
    for puck_name in PUCK_NAMES:
        n = (adata.obs['puck_id'] == puck_name).sum()
        log(f"    {puck_name}: {n:,}")

    banner("STEP 04 COMPLETE")
    log(f"  Output: {DIR_04_ANNOTATION}/")
    log(f"  Review: {os.path.join(DIR_04_ANNOTATION, 'popv_vs_c2l_crosstab.tsv')}")
    log(f"          -> complete the popV side of LABEL_HARMONIZE, rerun")
    log(f"             (cell2location abundances are cached, so it is fast)")
    log(f"  Sweep:  set C2L_CFG['detection_alpha']=200 and rerun to compare")


if __name__ == "__main__":
    main()

