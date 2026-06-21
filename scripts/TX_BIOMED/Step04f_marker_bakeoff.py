#!/usr/bin/env python3
"""
Step04f_marker_bakeoff_v2.py
=========================================================================
(See v1 header for full MARKER PANEL PROVENANCE; unchanged in v2.)

PURPOSE
  Adjudicate popV-vs-cell2location disagreements with an independent marker
  signal, producing `unified_annotation`. popV, cell2location, and consensus
  columns are preserved untouched.

RESOLUTION (changed in v2)
  unified_annotation is now a PER-BEAD call (argmax over the 8 compartment
  marker scores), gated so beads whose top score is not above background
  (<= 0) are labeled 'ambiguous'. The cluster-level margin-weighted vote is
  retained only as a diagnostic (cluster_marker_label / _confidence), NOT as
  the annotation.

CHANGELOG
  v1  cluster-level rollup. REJECTED after diagnostic: sparse beads do not
      cluster into cell types, so all 23 Leiden clusters were majority
      epithelial or myeloid and the minority compartments (T_cell, B_cell,
      endothelial, smooth_muscle, mast) received zero beads. The runner-up
      shares confirmed those populations live as 20-28% minorities inside the
      large clusters, i.e. resolvable per bead but not per cluster.
  v2  unified_annotation = per-bead argmax gated at top_score>0; cluster vote
      demoted to diagnostic; fixed categorical .map crash in the T-cell
      diagnostic (cast to str before mapping).
=========================================================================
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_ROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
DIR_04   = os.path.join(PROJECT_ROOT, "data/outputs/04_annotation")
BAKEOFF  = os.path.join(DIR_04, "bakeoff")
ANNOTATED = os.path.join(DIR_04, "all_pucks_annotated.h5ad")
INTERSECT = os.path.join(BAKEOFF, "marker_panel_intersection.tsv")
OUT_H5AD  = os.path.join(DIR_04, "all_pucks_annotated_unified.h5ad")
os.makedirs(BAKEOFF, exist_ok=True)

DPI = 300
COMPARTMENTS = ["epithelial", "T_cell", "B_cell", "myeloid",
                "fibroblast", "endothelial", "smooth_muscle", "mast"]
PLOT_CATS = COMPARTMENTS + ["ambiguous"]
COLORS = {
    "epithelial": "#4C72B0", "T_cell": "#DD8452", "B_cell": "#55A868",
    "myeloid": "#C44E52", "fibroblast": "#8172B3", "endothelial": "#937860",
    "smooth_muscle": "#DA8BC3", "mast": "#000000", "ambiguous": "#CCCCCC",
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
PLASMA_ADD = ["MZB1", "DERL3", "XBP1", "TNFRSF17", "PRDM1"]


def log(m): print(m, flush=True)


def antigen_presentation(g):
    return g.startswith("HLA-D") or g == "CD74"


log("=" * 70 + "\nSTEP04f v2  marker bake-off\n" + "=" * 70)

inter = pd.read_csv(INTERSECT, sep="\t")
proposed = {r["compartment"]: [g for g in str(r["proposed_panel"]).split(";") if g]
            for _, r in inter.iterrows()}
augmented = {c: list(dict.fromkeys(proposed.get(c, []))) for c in COMPARTMENTS}
augmented["B_cell"] = list(dict.fromkeys(augmented["B_cell"] + PLASMA_ADD))

counts = {}
for c in COMPARTMENTS:
    for g in augmented[c]:
        counts.setdefault(g, set()).add(c)
cross_panel = {g for g, comps in counts.items() if len(comps) >= 2}
ap_block = {g for g in counts if antigen_presentation(g)}
prune = cross_panel | ap_block
log(f"  cross-panel pruned: {sorted(cross_panel)}")
log(f"  antigen-presentation pruned: {sorted(ap_block)}")

log(f"\nloading {ANNOTATED}")
adata = sc.read_h5ad(ANNOTATED)
log(f"  {adata.n_obs:,} beads x {adata.n_vars:,} genes")
spatial_genes = set(map(str, adata.var_names))

FINAL, prov_rows, pruned_rows = {}, [], []
for c in COMPARTMENTS:
    keep = []
    for g in augmented[c]:
        if g in prune:
            reason = ("cross_panel_shared" if g in cross_panel else "") + \
                     ("|antigen_presentation" if g in ap_block else "")
            pruned_rows.append({"gene": g, "compartment": c, "reason": reason.strip("|")})
            continue
        if g not in spatial_genes:
            pruned_rows.append({"gene": g, "compartment": c, "reason": "absent_from_spatial_panel"})
            continue
        keep.append(g)
        src = ("curated" if g in CURATED[c]
               else "plasma_canonical" if (c == "B_cell" and g in PLASMA_ADD)
               else "reference_sig_and_panglaodb")
        prov_rows.append({"gene": g, "compartment": c, "source": src})
    FINAL[c] = keep
    log(f"  {c:14s} n={len(keep):2d}  {keep}")

pd.DataFrame(prov_rows).to_csv(os.path.join(BAKEOFF, "final_marker_panels.tsv"), sep="\t", index=False)
pd.DataFrame(pruned_rows).to_csv(os.path.join(BAKEOFF, "pruned_genes.tsv"), sep="\t", index=False)

log("\nrecomputing CPM+log1p from counts for scoring")
score_ad = adata.copy()
score_ad.X = score_ad.layers["counts"].copy()
score_ad.uns.pop("log1p", None)
sc.pp.normalize_total(score_ad, target_sum=1e4)
sc.pp.log1p(score_ad)
for c in COMPARTMENTS:
    sc.tl.score_genes(score_ad, FINAL[c], score_name=f"_s_{c}", use_raw=False, ctrl_size=50)
    adata.obs[f"marker_score_{c}"] = score_ad.obs[f"_s_{c}"].values
del score_ad

# ---- per-bead call (PRIMARY) -------------------------------------------
S = adata.obs[[f"marker_score_{c}" for c in COMPARTMENTS]].to_numpy()
order = np.argsort(-S, axis=1)
rows = np.arange(S.shape[0])
top_score = S[rows, order[:, 0]]
margin = top_score - S[rows, order[:, 1]]
argmax_comp = np.array(COMPARTMENTS)[order[:, 0]]

adata.obs["marker_argmax"] = pd.Categorical(argmax_comp, categories=COMPARTMENTS)
adata.obs["marker_top_score"] = top_score
adata.obs["marker_margin"] = margin

assigned = np.where(top_score > 0, argmax_comp, "ambiguous")
adata.obs["unified_annotation"] = pd.Categorical(assigned, categories=PLOT_CATS)
adata.obs["unified_confidence"] = margin
adata.obs["unified_basis"] = np.where(top_score > 0, "per_bead_marker_argmax",
                                      "ambiguous_nonpositive_top_score")

# ---- cluster vote (DIAGNOSTIC ONLY) ------------------------------------
clusters = adata.obs["clusters"].astype(str).values
w = np.clip(margin, 0, None)
vote = pd.DataFrame({"cluster": clusters, "comp": argmax_comp, "w": w})
W = (vote.groupby(["cluster", "comp"])["w"].sum()
     .unstack(fill_value=0.0).reindex(columns=COMPARTMENTS, fill_value=0.0))
cl_label = W.idxmax(axis=1)
totals = W.sum(axis=1).replace(0, np.nan)
cl_conf = (W.max(axis=1) / totals).fillna(0.0)
adata.obs["cluster_marker_label"] = pd.Series(clusters, index=adata.obs.index).map(cl_label).astype("category")
adata.obs["cluster_marker_confidence"] = pd.Series(clusters, index=adata.obs.index).map(cl_conf).astype(float)
pd.DataFrame({"cluster": W.index, "n_beads": vote.groupby("cluster").size().reindex(W.index).values,
             "cluster_label": cl_label.values, "confidence": cl_conf.values}).to_csv(
    os.path.join(BAKEOFF, "cluster_assignment.tsv"), sep="\t", index=False)

log("\n=== per-bead marker argmax (ungated) ===")
log(pd.Series(argmax_comp).value_counts().to_string())
log("\n=== unified_annotation (per-bead, gated top_score>0) ===")
log(adata.obs["unified_annotation"].value_counts().to_string())
log(f"  ambiguous (top score <= 0): {(top_score <= 0).sum():,} ({100*(top_score <= 0).mean():.1f}%)")
log("\n=== cluster-level vote (DIAGNOSTIC ONLY, not the annotation) ===")
log(adata.obs["cluster_marker_label"].value_counts().to_string())

# ---- T-cell disagreement resolution ------------------------------------
def parse_consensus(v):
    v = str(v)
    return set(v.split("ambiguous:", 1)[1].split("|")) if v.startswith("ambiguous:") else {v}

con = adata.obs["consensus_annotation"].astype(str).map(parse_consensus)
is_amb = adata.obs["consensus_annotation"].astype(str).str.startswith("ambiguous:")
t_disagree = (is_amb & con.map(lambda s: "T_cell" in s)).to_numpy()
log("\n=== resolution of T-cell disagreement beads ===")
log(f"  T-involving disagreement beads: {int(t_disagree.sum()):,}")
if t_disagree.sum():
    log(adata.obs.loc[t_disagree, "unified_annotation"].value_counts().to_string())
log(f"\n  total unified == T_cell (all beads): "
    f"{int((adata.obs['unified_annotation'] == 'T_cell').sum()):,}")

# ---- figures -----------------------------------------------------------
def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(BAKEOFF, f"Step04f_{stem}.{ext}"), dpi=DPI,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)

if "X_umap" in adata.obsm:
    fig, ax = plt.subplots(figsize=(10, 9))
    for c in PLOT_CATS:
        m = (adata.obs["unified_annotation"] == c).to_numpy()
        ax.scatter(adata.obsm["X_umap"][m, 0], adata.obsm["X_umap"][m, 1],
                   s=2, c=COLORS[c], label=c, linewidths=0)
    ax.axis("off"); ax.legend(markerscale=6, fontsize=14, frameon=False)
    save(fig, "unified_umap")

puck_col = next((c for c in ["puck", "sample", "library", "batch", "sample_id"]
                 if c in adata.obs.columns), None)
if {"x_coord", "y_coord"}.issubset(adata.obs.columns):
    pucks = adata.obs[puck_col].astype(str).unique() if puck_col else ["all"]
    fig, axes = plt.subplots(1, len(pucks), figsize=(7 * len(pucks), 7), squeeze=False)
    for j, pk in enumerate(pucks):
        ax = axes[0, j]
        sub = adata.obs if pk == "all" else adata.obs[adata.obs[puck_col].astype(str) == pk]
        for c in PLOT_CATS:
            m = sub["unified_annotation"] == c
            ax.scatter(sub.loc[m, "x_coord"], sub.loc[m, "y_coord"], s=2, c=COLORS[c], linewidths=0)
        ax.set_title(pk, fontsize=16); ax.axis("off"); ax.set_aspect("equal")
    save(fig, "unified_spatial")

mean_score = adata.obs.groupby("clusters")[[f"marker_score_{c}" for c in COMPARTMENTS]].mean()
mean_score.columns = COMPARTMENTS
fig, ax = plt.subplots(figsize=(10, max(6, 0.4 * len(mean_score))))
im = ax.imshow(mean_score.to_numpy(), aspect="auto", cmap="magma")
ax.set_xticks(range(len(COMPARTMENTS))); ax.set_xticklabels(COMPARTMENTS, rotation=45, ha="right")
ax.set_yticks(range(len(mean_score))); ax.set_yticklabels(mean_score.index)
fig.colorbar(im, ax=ax, label="mean marker score")
save(fig, "cluster_score_heatmap")

fig, ax = plt.subplots(figsize=(9, 6))
vc = adata.obs["unified_annotation"].value_counts().reindex(PLOT_CATS).fillna(0)
ax.bar(range(len(PLOT_CATS)), vc.values, color=[COLORS[c] for c in PLOT_CATS])
ax.set_xticks(range(len(PLOT_CATS))); ax.set_xticklabels(PLOT_CATS, rotation=45, ha="right")
ax.set_ylabel("beads")
save(fig, "unified_composition")

log(f"\nwriting {OUT_H5AD}")
adata.write_h5ad(OUT_H5AD)
log("DONE.")
