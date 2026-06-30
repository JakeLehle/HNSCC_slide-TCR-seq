#!/usr/bin/env python3
"""
Diagnostic_marker_DEG_and_Tcell_mask.py
=========================================================================
READ-ONLY DIAGNOSTIC. Does NOT modify the unified annotation.

Reads:
  data/outputs/04_annotation/all_pucks_annotated_unified.h5ad   (locked)
  data/outputs/04_annotation/bakeoff/final_marker_panels.tsv    (04f panels)
Writes only to:
  data/outputs/04_annotation/diagnostics/                       (TSVs + figures)

PURPOSE
  Gather, in one pass, the evidence needed to decide what (if anything) to
  change in Step04f before any rerun. The annotation->validation loop has a
  feedback edge (findings here may change the 04f marker panel), so this is a
  throwaway diagnostic, not the eventual Step05. Once we read the output we
  fold all 04f edits in together and rerun 04f once.

WHAT IT PRODUCES
  A. Compartment DEG (rank_genes_groups, Wilcoxon, on log1p-CPM X), each gene
     flagged panel vs novel. The novel genes are the publication cross-check.
  B. clusters x unified_annotation cross-tab (row-normalized), confirming the
     minority compartments live as within-cluster minorities (the reason the
     annotation is per-bead, not per-cluster).
  C. T-cell mask QC, built for a clear judgment:
       C1 mask size per puck
       C2 CD3 detection rate in T beads vs background (with fold enrichment)
       C3 NK-marker enrichment in T beads vs background (the leakage tell)
       C4 CD3+/NK- vs CD3-/NK+ breakdown among the T calls
       C5 PTPRC-drop flip test: recompute the per-bead argmax with PTPRC
          removed from the T panel and report how many T calls survive and
          where the flips go. This is the single number for whether PTPRC is
          inflating the mask.
       C6 margin distribution for T beads (confident vs barely over the gate)
       C7 popV / cell2location T-lineage concordance (triangulation)
       C8 spatial homotypic enrichment per puck (coherent foci vs scatter),
          which also previews the TLS check (task 4.3)
  D. A DECISION GUIDE printed at the end.

CONVENTIONS
  Self-contained (mirrors 04f; no spatial_config import). matplotlib Agg.
  Hex colors. PDF + PNG at 300 DPI. No em dashes. Scoring for the flip test
  replicates 04f exactly (CPM target_sum 1e4 + log1p from layers['counts'],
  score_genes use_raw=False, ctrl_size=50) so the baseline argmax can be
  checked against the stored unified_annotation before any flips are read.

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Slide-TCR-seq (Sophia Liu / Ragon collaboration)
=========================================================================
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --------------------------------------------------------------------- paths
PROJECT_ROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
DIR_04   = os.path.join(PROJECT_ROOT, "data/outputs/04_annotation")
BAKEOFF  = os.path.join(DIR_04, "bakeoff")
DIAG     = os.path.join(DIR_04, "diagnostics")
UNIFIED  = os.path.join(DIR_04, "all_pucks_annotated_unified.h5ad")
PANELS   = os.path.join(BAKEOFF, "final_marker_panels.tsv")
os.makedirs(DIAG, exist_ok=True)

# ----------------------------------------------------------------- constants
DPI = 300
TOPN_DEG = 50          # rows kept per compartment in the DEG table
TOPN_NOVEL_PRINT = 12  # novel genes printed to console per compartment
K_SPATIAL = 15         # neighbors for spatial homotypic enrichment
N_PERM = 200           # label permutations for the spatial null

COMPARTMENTS = ["epithelial", "T_cell", "B_cell", "myeloid",
                "fibroblast", "endothelial", "smooth_muscle", "mast"]
COLORS = {
    "epithelial": "#4C72B0", "T_cell": "#DD8452", "B_cell": "#55A868",
    "myeloid": "#C44E52", "fibroblast": "#8172B3", "endothelial": "#937860",
    "smooth_muscle": "#DA8BC3", "mast": "#000000", "ambiguous": "#CCCCCC",
}

# leakage probes
CD3_GENES = ["CD3D", "CD3E", "CD3G", "CD247"]
NK_GENES  = ["NKG7", "GNLY", "KLRD1", "KLRF1", "KLRB1", "NCAM1", "FCGR3A", "NCR1"]
PTPRC     = "PTPRC"

# figure fonts (paper convention is 28 to 34; gene-dense dotplot ticks are
# reduced for legibility since this is an internal diagnostic, not a panel)
plt.rcParams.update({"font.size": 28, "axes.titlesize": 32,
                     "axes.labelsize": 30, "savefig.facecolor": "white"})


def log(m):
    print(m, flush=True)


def banner(t):
    log("\n" + "=" * 74 + f"\n{t}\n" + "=" * 74)


def save(fig, stem):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(DIAG, f"diag_{stem}.{ext}"),
                    dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# ===================================================================== load
banner("LOAD")
log(f"reading {UNIFIED}")
adata = sc.read_h5ad(UNIFIED)
log(f"  {adata.n_obs:,} beads x {adata.n_vars:,} genes")
log(f"  obs has: {[c for c in ['unified_annotation','clusters','puck_id','consensus_annotation','c2l_argmax','popv_prediction'] if c in adata.obs.columns]}")

panels = pd.read_csv(PANELS, sep="\t")
FINAL = {c: panels.loc[panels["compartment"] == c, "gene"].tolist() for c in COMPARTMENTS}
for c in COMPARTMENTS:
    log(f"  panel {c:14s} n={len(FINAL[c]):2d}")

# resolve the puck column
puck_col = next((c for c in ["puck_id", "puck", "sample", "library", "batch", "sample_id"]
                 if c in adata.obs.columns), None)
log(f"  puck column: {puck_col}")
have_coords = {"x_coord", "y_coord"}.issubset(adata.obs.columns)

var_names = pd.Index(adata.var_names.astype(str))
def present(genes):
    return [g for g in genes if g in var_names]


# ======================================================== A. compartment DEG
banner("A. COMPARTMENT DEG (panel vs novel)")
de = adata[adata.obs["unified_annotation"].isin(COMPARTMENTS)].copy()
de.obs["unified_annotation"] = de.obs["unified_annotation"].astype(str).astype("category")
log(f"  beads in DEG (ambiguous excluded): {de.n_obs:,}")
log("  rank_genes_groups: wilcoxon, use_raw=False, pts=True (on log1p-CPM X)")
sc.tl.rank_genes_groups(de, "unified_annotation", method="wilcoxon",
                        use_raw=False, pts=True)
res = de.uns["rank_genes_groups"]
groups = list(res["names"].dtype.names)
pts = res.get("pts")
pts_rest = res.get("pts_rest")

rows = []
for g in groups:
    panel_set = set(FINAL.get(g, []))
    for rank in range(min(TOPN_DEG, len(res["names"][g]))):
        gn = str(res["names"][g][rank])
        pin = float(pts.loc[gn, g]) if (pts is not None and gn in pts.index) else np.nan
        pre = float(pts_rest.loc[gn, g]) if (pts_rest is not None and gn in pts_rest.index) else np.nan
        rows.append({
            "compartment": g, "rank": rank + 1, "gene": gn,
            "logFC": float(res["logfoldchanges"][g][rank]),
            "pval_adj": float(res["pvals_adj"][g][rank]),
            "score": float(res["scores"][g][rank]),
            "pct_in": pin, "pct_rest": pre,
            "in_panel": gn in panel_set,
        })
deg = pd.DataFrame(rows)
deg.to_csv(os.path.join(DIAG, "compartment_DEG.tsv"), sep="\t", index=False)
log(f"  wrote compartment_DEG.tsv ({len(deg)} rows)")

log("\n  top NOVEL (non-panel) markers per compartment "
    "(rank, gene, logFC, pct_in, pct_rest):")
for g in COMPARTMENTS:
    sub = deg[(deg["compartment"] == g) & (~deg["in_panel"])].head(TOPN_NOVEL_PRINT)
    log(f"  --- {g} ---")
    for _, r in sub.iterrows():
        log(f"      {int(r['rank']):>3}  {r['gene']:<12} "
            f"lfc={r['logFC']:+.2f}  in={r['pct_in']:.2f}  rest={r['pct_rest']:.2f}")

# compact diagnostic dotplot: top 4 markers per compartment, dot size = pct_in,
# color = column-z of mean log1p-CPM across compartments
top4 = {g: deg[deg["compartment"] == g].head(4)["gene"].tolist() for g in COMPARTMENTS}
gene_order, gene_home = [], []
for g in COMPARTMENTS:
    for gn in top4[g]:
        if gn not in gene_order:
            gene_order.append(gn); gene_home.append(g)
present_order = [gn for gn in gene_order if gn in var_names]
if present_order:
    Xsub = de[:, present_order].X
    Xsub = Xsub.toarray() if sparse.issparse(Xsub) else np.asarray(Xsub)
    csub = de[:, present_order].layers["counts"]
    csub = csub.toarray() if sparse.issparse(csub) else np.asarray(csub)
    lab = de.obs["unified_annotation"].to_numpy()
    mean_expr = np.zeros((len(COMPARTMENTS), len(present_order)))
    frac = np.zeros_like(mean_expr)
    for i, comp in enumerate(COMPARTMENTS):
        m = lab == comp
        mean_expr[i] = Xsub[m].mean(axis=0)
        frac[i] = (csub[m] > 0).mean(axis=0)
    z = (mean_expr - mean_expr.mean(0, keepdims=True)) / (mean_expr.std(0, keepdims=True) + 1e-9)
    fig, ax = plt.subplots(figsize=(max(16, 0.7 * len(present_order)), 9))
    for i in range(len(COMPARTMENTS)):
        for j in range(len(present_order)):
            ax.scatter(j, i, s=20 + 600 * frac[i, j], c=[plt.cm.magma((z[i, j] + 3) / 6)],
                       edgecolors="#333333", linewidths=0.4)
    ax.set_yticks(range(len(COMPARTMENTS))); ax.set_yticklabels(COMPARTMENTS, fontsize=22)
    ax.set_xticks(range(len(present_order)))
    ax.set_xticklabels(present_order, rotation=90, fontsize=14)  # diagnostic ticks
    ax.set_title("Top markers per compartment (size=pct expressing, color=z mean expr)",
                 fontsize=22)
    ax.set_xlim(-1, len(present_order)); ax.set_ylim(-1, len(COMPARTMENTS))
    save(fig, "compartment_DEG_dotplot")
    log("  wrote diag_compartment_DEG_dotplot.{pdf,png}")
del de


# ===================================================== B. cluster cross-tab
banner("B. clusters x unified_annotation cross-tab (row-normalized)")
ct = pd.crosstab(adata.obs["clusters"], adata.obs["unified_annotation"])
ct = ct.reindex(columns=[c for c in (COMPARTMENTS + ["ambiguous"]) if c in ct.columns], fill_value=0)
ct_row = ct.div(ct.sum(1), axis=0)
ct_row.to_csv(os.path.join(DIAG, "cluster_x_unified_rownorm.tsv"), sep="\t")
log("  per-cluster composition (fraction; T_cell column highlighted in log):")
for cl in ct_row.index:
    tfrac = ct_row.loc[cl, "T_cell"] if "T_cell" in ct_row.columns else np.nan
    log(f"    cluster {str(cl):>3}  n={int(ct.loc[cl].sum()):>6}  T_cell_frac={tfrac:.2f}")
fig, ax = plt.subplots(figsize=(10, max(7, 0.42 * len(ct_row))))
im = ax.imshow(ct_row.to_numpy(), aspect="auto", cmap="magma", vmin=0, vmax=1)
ax.set_xticks(range(ct_row.shape[1])); ax.set_xticklabels(ct_row.columns, rotation=45, ha="right", fontsize=20)
ax.set_yticks(range(ct_row.shape[0])); ax.set_yticklabels(ct_row.index, fontsize=16)
fig.colorbar(im, ax=ax, label="fraction of cluster")
ax.set_title("cluster composition", fontsize=24)
save(fig, "cluster_x_unified_heatmap")
log("  wrote cluster_x_unified_rownorm.tsv + heatmap")


# =============================== rebuild scores for C4/C5 (replicate 04f) ===
banner("rebuild marker scores from counts (replicate 04f for the flip test)")
score_ad = adata.copy()
score_ad.X = score_ad.layers["counts"].copy()
score_ad.uns.pop("log1p", None)
sc.pp.normalize_total(score_ad, target_sum=1e4)
sc.pp.log1p(score_ad)

def score_panel(genes, name):
    gp = present(genes)
    sc.tl.score_genes(score_ad, gp, score_name=name, use_raw=False, ctrl_size=50)
    return score_ad.obs[name].to_numpy()

base = {c: score_panel(FINAL[c], f"_b_{c}") for c in COMPARTMENTS}
S_base = np.column_stack([base[c] for c in COMPARTMENTS])
order = np.argsort(-S_base, axis=1)
top_base = S_base[np.arange(len(S_base)), order[:, 0]]
arg_base = np.array(COMPARTMENTS)[order[:, 0]]
label_base = np.where(top_base > 0, arg_base, "ambiguous")

stored = adata.obs["unified_annotation"].astype(str).to_numpy()
agree = (label_base == stored).mean()
log(f"  baseline argmax vs stored unified_annotation agreement: {agree*100:.2f}%")
if agree < 0.99:
    log("  WARNING: baseline does not reproduce stored labels at >=99%. "
        "Read the flip test against BASELINE (recomputed) labels, not stored.")

# NK signature score (composite leakage probe)
nk_score = score_panel(NK_GENES, "_nk_score")


# ======================================================= C. T-cell mask QC
banner("C. T-CELL MASK QC")
t_mask = stored == "T_cell"
n_t = int(t_mask.sum())
log(f"C1  T_cell beads (stored): {n_t:,}")
if puck_col:
    for pk, n in adata.obs.loc[t_mask, puck_col].value_counts().sort_index().items():
        log(f"      {pk}: {int(n):,}")

# dense probe matrix (counts) for CD3 + NK
probe_genes = present(CD3_GENES + NK_GENES)
missing = [g for g in (CD3_GENES + NK_GENES) if g not in var_names]
if missing:
    log(f"  probes absent from panel (skipped): {missing}")
cnt = adata[:, probe_genes].layers["counts"]
cnt = cnt.toarray() if sparse.issparse(cnt) else np.asarray(cnt)
det = cnt > 0
gidx = {g: i for i, g in enumerate(probe_genes)}
bg = ~t_mask  # background = every non-T bead

def rate(gene, mask):
    return det[mask, gidx[gene]].mean() if gene in gidx else np.nan

log("\nC2  CD3 detection (counts>0): rate_in_T  rate_bg  fold")
cd3_rows = []
for g in CD3_GENES:
    if g not in gidx:
        continue
    rt, rb = rate(g, t_mask), rate(g, bg)
    fold = rt / rb if rb > 0 else np.inf
    cd3_rows.append((g, rt, rb, fold))
    log(f"      {g:<8} {rt:.3f}     {rb:.3f}    {fold:.1f}x")

log("\nC3  NK-marker detection: rate_in_T  rate_bg  fold  (leakage tell)")
nk_rows = []
for g in NK_GENES:
    if g not in gidx:
        continue
    rt, rb = rate(g, t_mask), rate(g, bg)
    fold = rt / rb if rb > 0 else np.inf
    nk_rows.append((g, rt, rb, fold))
    log(f"      {g:<8} {rt:.3f}     {rb:.3f}    {fold:.1f}x")

idx_cd3 = [gidx[g] for g in CD3_GENES if g in gidx]
idx_nk = [gidx[g] for g in NK_GENES if g in gidx]
anyCD3 = det[:, idx_cd3].any(axis=1) if idx_cd3 else np.zeros(adata.n_obs, bool)
anyNK = det[:, idx_nk].any(axis=1) if idx_nk else np.zeros(adata.n_obs, bool)
tt = t_mask
b_cd3pos_nkneg = int((tt & anyCD3 & ~anyNK).sum())
b_cd3neg_nkpos = int((tt & ~anyCD3 & anyNK).sum())
b_both = int((tt & anyCD3 & anyNK).sum())
b_neither = int((tt & ~anyCD3 & ~anyNK).sum())
log("\nC4  among T_cell beads (detection-based; dropout inflates 'neither'):")
log(f"      CD3+ NK-  : {b_cd3pos_nkneg:>6}  ({100*b_cd3pos_nkneg/n_t:.1f}%)  bona fide T")
log(f"      CD3- NK+  : {b_cd3neg_nkpos:>6}  ({100*b_cd3neg_nkpos/n_t:.1f}%)  likely NK leakage")
log(f"      CD3+ NK+  : {b_both:>6}  ({100*b_both/n_t:.1f}%)")
log(f"      neither   : {b_neither:>6}  ({100*b_neither/n_t:.1f}%)  dropout-dominated")
log(f"      NK signature score in T beads: median={np.median(nk_score[tt]):+.3f} "
    f"vs background median={np.median(nk_score[bg]):+.3f}")

# C5 PTPRC-drop flip test (only the T score changes; others held at baseline)
banner("C5. PTPRC-DROP FLIP TEST")
t_noptprc_genes = [g for g in FINAL["T_cell"] if g != PTPRC]
if PTPRC not in FINAL["T_cell"]:
    log(f"  note: {PTPRC} is not in the realized T panel; flip test is a no-op")
t_np = score_panel(t_noptprc_genes, "_t_noptprc")
S_np = S_base.copy()
S_np[:, COMPARTMENTS.index("T_cell")] = t_np
order_np = np.argsort(-S_np, axis=1)
top_np = S_np[np.arange(len(S_np)), order_np[:, 0]]
arg_np = np.array(COMPARTMENTS)[order_np[:, 0]]
label_np = np.where(top_np > 0, arg_np, "ambiguous")

base_t = label_base == "T_cell"
n_base_t = int(base_t.sum())
stay = int((base_t & (label_np == "T_cell")).sum())
flipped = label_np[base_t & (label_np != "T_cell")]
log(f"  baseline T_cell beads: {n_base_t:,}")
log(f"  remain T_cell after PTPRC removal: {stay:,} ({100*stay/max(n_base_t,1):.1f}%)")
log(f"  flipped away: {n_base_t - stay:,} ({100*(n_base_t-stay)/max(n_base_t,1):.1f}%)")
if len(flipped):
    log("  flips go to:")
    for k, v in pd.Series(flipped).value_counts().items():
        log(f"      {k:<14} {int(v):,}")
flip_tbl = pd.Series(label_np[base_t]).value_counts().rename_axis("new_label").reset_index(name="n_beads")
flip_tbl.to_csv(os.path.join(DIAG, "ptprc_flip_table.tsv"), sep="\t", index=False)

# C6 margin distribution for T beads
banner("C6. MARGIN (confidence) DISTRIBUTION FOR T BEADS")
if "unified_confidence" in adata.obs.columns:
    mg = adata.obs.loc[t_mask, "unified_confidence"].to_numpy()
    qs = np.percentile(mg, [10, 25, 50, 75, 90])
    log(f"  margin percentiles (10/25/50/75/90): "
        f"{qs[0]:.3f} {qs[1]:.3f} {qs[2]:.3f} {qs[3]:.3f} {qs[4]:.3f}")
    log(f"  T beads with margin < 0.02 (barely over gate): "
        f"{int((mg < 0.02).sum()):,} ({100*(mg<0.02).mean():.1f}%)")
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.hist(mg, bins=60, color=COLORS["T_cell"])
    ax.set_xlabel("margin (top1 - top2)"); ax.set_ylabel("T beads")
    ax.set_title("T-cell call confidence", fontsize=24)
    save(fig, "tcell_margin_hist")

# C7 popV / c2l T-lineage concordance
banner("C7. popV / cell2location T-LINEAGE CONCORDANCE (T beads)")
def is_T_label(s):
    s = str(s).lower()
    return ("t cell" in s) or ("cd4" in s and "t" in s) or ("cd8" in s and "t" in s) or ("regulatory t" in s)

if "popv_prediction" in adata.obs.columns:
    popv_T = adata.obs["popv_prediction"].map(is_T_label).to_numpy()
    log(f"  popV called T-lineage on T beads: "
        f"{100*popv_T[t_mask].mean():.1f}%")
if "c2l_argmax" in adata.obs.columns:
    c2l_T = adata.obs["c2l_argmax"].map(is_T_label).to_numpy()
    log(f"  c2l argmax called T on T beads:    {100*c2l_T[t_mask].mean():.1f}%")
if "popv_prediction" in adata.obs.columns and "c2l_argmax" in adata.obs.columns:
    either = (popv_T | c2l_T)[t_mask].mean()
    both = (popv_T & c2l_T)[t_mask].mean()
    log(f"  either popV or c2l agrees:         {100*either:.1f}%")
    log(f"  both agree:                        {100*both:.1f}%")

# C8 spatial homotypic enrichment per puck
banner("C8. SPATIAL HOMOTYPIC ENRICHMENT OF T BEADS (per puck)")
spatial_rows = []
if have_coords and puck_col:
    for pk in sorted(adata.obs[puck_col].astype(str).unique()):
        m = (adata.obs[puck_col].astype(str) == pk).to_numpy()
        coords = adata.obs.loc[m, ["x_coord", "y_coord"]].to_numpy()
        is_t = (stored[m] == "T_cell")
        gfrac = is_t.mean()
        if is_t.sum() < 5 or len(coords) <= K_SPATIAL:
            log(f"  {pk}: too few T beads for enrichment"); continue
        nn = NearestNeighbors(n_neighbors=K_SPATIAL + 1).fit(coords)
        _, nbr = nn.kneighbors(coords)
        nbr = nbr[:, 1:]  # drop self
        homot = is_t[nbr].mean(axis=1)            # per-bead neighbor T fraction
        obs_stat = homot[is_t].mean()             # mean over T beads
        rng = np.random.default_rng(42)
        perm = np.empty(N_PERM)
        for i in range(N_PERM):
            sh = rng.permutation(is_t)
            perm[i] = sh[nbr].mean(axis=1)[sh].mean()
        z = (obs_stat - perm.mean()) / (perm.std() + 1e-9)
        emp_p = (np.sum(perm >= obs_stat) + 1) / (N_PERM + 1)
        enr = obs_stat / gfrac if gfrac > 0 else np.nan
        log(f"  {pk}: global T frac={gfrac:.3f}  neighbor T frac (T beads)={obs_stat:.3f}  "
            f"enrichment={enr:.2f}x  z={z:.1f}  emp_p={emp_p:.3f}")
        spatial_rows.append({"puck": pk, "global_T_frac": gfrac,
                             "neighbor_T_frac": obs_stat, "enrichment": enr,
                             "z": z, "emp_p": emp_p})
        # T-only spatial scatter for eyeballing coherence
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(coords[~is_t, 0], coords[~is_t, 1], s=1, c="#E5E5E5", linewidths=0)
        ax.scatter(coords[is_t, 0], coords[is_t, 1], s=4, c=COLORS["T_cell"], linewidths=0)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{pk}  T beads", fontsize=30)
        save(fig, f"tcell_spatial_{pk}")
    if spatial_rows:
        pd.DataFrame(spatial_rows).to_csv(
            os.path.join(DIAG, "tcell_spatial_enrichment.tsv"), sep="\t", index=False)
else:
    log("  coordinates or puck column unavailable; skipped")

# CD3 vs NK enrichment bar (one figure summarizing C2/C3)
if cd3_rows or nk_rows:
    labels = [r[0] for r in cd3_rows] + [r[0] for r in nk_rows]
    folds = [r[3] for r in cd3_rows] + [r[3] for r in nk_rows]
    cols = [COLORS["T_cell"]] * len(cd3_rows) + ["#777777"] * len(nk_rows)
    folds = [min(f, 50) for f in folds]  # cap inf/extreme for display
    fig, ax = plt.subplots(figsize=(max(10, 0.8 * len(labels)), 7))
    ax.bar(range(len(labels)), folds, color=cols)
    ax.axhline(1.0, color="#000000", lw=1, ls="--")
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=90, fontsize=18)
    ax.set_ylabel("fold detection vs background")
    ax.set_title("CD3 (orange) vs NK (grey) in T beads", fontsize=22)
    save(fig, "cd3_vs_nk_enrichment")


# ============================================================ D. decision guide
banner("D. DECISION GUIDE (read these to set the 04f edits)")
log(f"""  1) PTPRC: if C5 'remain T_cell' is high (say >90%), PTPRC is not
     inflating the mask, keep it. If a large share flips to myeloid/B/
     ambiguous, drop PTPRC from the 04f T panel.
  2) NK leakage: if C3 NK folds and the C4 'CD3- NK+' share are high, and
     C7 concordance is low, the T mask is pulling in NK. Options for 04f:
     add an NK-exclusion rule, or a CD3-confirmation gate on T calls.
  3) Spatial (C8): enrichment >1 with z well above 0 suggests T beads form
     coherent foci (real structure / candidate TLS), supporting the mask;
     near 1 suggests scatter and weaker support.
  4) DEG novel genes (section A): scan the T_cell novel list. CD3/cytotoxic
     genes corroborate; a wall of NKG7/GNLY/KLRD1 corroborates leakage.
  All outputs in: {DIAG}
""")
log("DONE (read-only; unified object untouched).")
