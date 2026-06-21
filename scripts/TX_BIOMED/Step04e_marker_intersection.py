#!/usr/bin/env python3
"""
Step04e_marker_intersection.py  (diagnostic; compute node, not login)

Three ways to define each compartment's markers, and their overlap, so we can
lock the bake-off panels before writing it:
  curated       : your literature-vetted pairs, pooled to 8 compartments
  reference_sig : Wilcoxon rank_genes_groups on the HNSCC reference
                  (12 final_annotation types -> 8 compartments), adj p < 0.05
  panglaodb     : union of mapped PanglaoDB Hs cell types, present in the panel

Env: sc_pre. Reference is 4.25 GB / 155k cells -> compute node.
Hardcoded TX Biomed paths; does NOT import spatial_config.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
REFERENCE = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_4/00_input/reference_for_cell2location.h5ad"
ANNOTATED = os.path.join(PROJECT_ROOT, "data/outputs/04_annotation/all_pucks_annotated.h5ad")
PG_LOCAL  = os.path.join(PROJECT_ROOT, "data/inputs/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv.gz")
OUTDIR    = os.path.join(PROJECT_ROOT, "data/outputs/04_annotation/bakeoff")
os.makedirs(OUTDIR, exist_ok=True)

ADJ_P = 0.05
TOP_N_REF = 50  # cap reference sig markers per compartment by score rank

REF_TO_COMPARTMENT = {
    "basal cell": "epithelial",
    "CD4-positive, alpha-beta T cell": "T_cell",
    "CD8-positive, alpha-beta T cell": "T_cell",
    "regulatory T cell": "T_cell",
    "B cell": "B_cell",
    "macrophage": "myeloid",
    "myeloid dendritic cell": "myeloid",
    "plasmacytoid dendritic cell": "myeloid",
    "fibroblast": "fibroblast",
    "endothelial cell": "endothelial",
    "smooth muscle cell": "smooth_muscle",
    "mast cell": "mast",
}

# PanglaoDB Hs 'cell type' -> compartment  (CONFIRM / EDIT)
PANGLAO_TO_COMPARTMENT = {
    "Basal cells": "epithelial", "Epithelial cells": "epithelial", "Keratinocytes": "epithelial",
    "T cells": "T_cell", "T cells naive": "T_cell", "T helper cells": "T_cell",
    "T memory cells": "T_cell", "T cytotoxic cells": "T_cell", "T regulatory cells": "T_cell",
    "T follicular helper cells": "T_cell", "Gamma delta T cells": "T_cell",
    "Natural killer T cells": "T_cell",
    "B cells": "B_cell", "B cells memory": "B_cell", "B cells naive": "B_cell",
    "Plasma cells": "B_cell",
    "Macrophages": "myeloid", "Monocytes": "myeloid", "Dendritic cells": "myeloid",
    "Plasmacytoid dendritic cells": "myeloid",
    "Fibroblasts": "fibroblast", "Stromal cells": "fibroblast",
    "Endothelial cells": "endothelial", "Endothelial cells (aorta)": "endothelial",
    "Smooth muscle cells": "smooth_muscle", "Vascular smooth muscle cells": "smooth_muscle",
    "Pericytes": "smooth_muscle",
    "Mast cells": "mast",
    # UNMAPPED on purpose (no reference/c2l counterpart): NK cells, Neutrophils,
    # Basophils, Eosinophils, Megakaryocytes, Myeloid-derived suppressor cells,
    # Myofibroblasts, Myoepithelial cells, tissue-specific macrophages.
}

CURATED = {
    "epithelial":    ["KRT5", "TACSTD2"],
    "T_cell":        ["CD3E", "IL7R", "CD8A", "CCL5", "TIGIT", "CTLA4"],
    "B_cell":        ["CD79A", "MS4A1"],
    "myeloid":       ["CD68", "SPI1", "LAMP3", "CCR7", "IL3RA", "LILRA4"],
    "fibroblast":    ["DCN", "COL1A1"],
    "endothelial":   ["PECAM1", "VWF"],
    "smooth_muscle": ["TAGLN", "RGS5"],
    "mast":          ["TPSAB1", "CPA3"],
}
COMPARTMENTS = list(CURATED.keys())


def log(m): print(m, flush=True)


# spatial panel genes (backed; cheap)
spatial_genes = set(map(str, sc.read_h5ad(ANNOTATED, backed='r').var_names))
log(f"spatial panel genes: {len(spatial_genes):,}")

# PanglaoDB per-compartment sets (Hs, present in panel)
pg = pd.read_table(PG_LOCAL)
pg = pg[pg['species'].astype(str).str.contains('Hs', na=False)]
GENE = 'official gene symbol'
panglao_sets = {c: set() for c in COMPARTMENTS}
for ct, sub in pg.groupby('cell type'):
    comp = PANGLAO_TO_COMPARTMENT.get(ct)
    if comp:
        panglao_sets[comp].update(
            g for g in sub[GENE].dropna().astype(str).unique() if g in spatial_genes
        )

# reference DE at compartment level
log("loading reference (4.25 GB)...")
ref = sc.read_h5ad(REFERENCE)
ref = ref[ref.obs['final_annotation'].isin(REF_TO_COMPARTMENT)].copy()
ref.obs['compartment'] = ref.obs['final_annotation'].map(REF_TO_COMPARTMENT).astype('category')
shared = [g for g in ref.var_names if str(g) in spatial_genes]
ref = ref[:, shared].copy()
ref.X = ref.layers['counts'].copy()
sc.pp.normalize_total(ref, target_sum=1e4)
sc.pp.log1p(ref)
log(f"reference for DE: {ref.n_obs:,} x {ref.n_vars:,}")
sc.tl.rank_genes_groups(ref, groupby='compartment', method='wilcoxon')

ref_sig = {}
for comp in ref.obs['compartment'].cat.categories:
    df = sc.get.rank_genes_groups_df(ref, group=comp)
    df = df[(df['pvals_adj'] < ADJ_P) & (df['logfoldchanges'] > 0)]
    ref_sig[comp] = df.sort_values('scores', ascending=False)['names'].astype(str).head(TOP_N_REF).tolist()

# three-way overlap + proposed panels
rows = []
for comp in COMPARTMENTS:
    cur = set(g for g in CURATED[comp] if g in spatial_genes)
    rs = set(ref_sig.get(comp, []))
    pgs = panglao_sets[comp]
    proposed = cur | (rs & pgs)
    log("\n" + "=" * 70 + f"\n  {comp}\n" + "=" * 70)
    log(f"  curated ({len(cur)}): {sorted(cur)}")
    log(f"  reference_sig adjP<{ADJ_P} top{TOP_N_REF} ({len(rs)}): {sorted(rs)}")
    log(f"  panglaodb ({len(pgs)} genes)")
    log(f"  curated also in panglaodb: {sorted(cur & pgs)}")
    log(f"  reference_sig AND panglaodb ({len(rs & pgs)}): {sorted(rs & pgs)}")
    log(f"  PROPOSED panel = curated + (ref_sig & panglaodb) ({len(proposed)}): {sorted(proposed)}")
    rows.append({
        'compartment': comp, 'n_curated': len(cur), 'n_ref_sig': len(rs),
        'n_panglaodb': len(pgs),
        'curated_in_panglaodb': ";".join(sorted(cur & pgs)),
        'ref_sig_and_panglaodb': ";".join(sorted(rs & pgs)),
        'proposed_panel': ";".join(sorted(proposed)),
    })

pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "marker_panel_intersection.tsv"), sep='\t', index=False)
log(f"\nsaved: {os.path.join(OUTDIR, 'marker_panel_intersection.tsv')}")
