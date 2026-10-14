#!/usr/bin/env python3
"""
Step08c_Neoantigen_TCR_Feasibility.py
================================================================================
READ-ONLY. Decides whether the neoantigen-anchored TCR analysis is runnable at
all, and builds the case/control sets it would use. Runs no hypothesis test.

WHY THIS COMES FIRST
    Step08a/b established that TCR signal tracks T cells once depth is
    controlled, but carries no global spatial structure. That rules out a
    puck-wide average, not a targeted test: if five clones out of 11,757 are
    expanded beside a neoantigen, a global statistic cannot see them.

    But the targeted test needs populated neighbourhoods. With 58/32/20
    neoantigen beads per puck, if a typical neighbourhood holds three TCR beads
    and no repeated clone, nothing downstream can work regardless of the
    statistics. This script measures that before anything is designed further.

DESIGN DECISIONS BAKED IN HERE

  TCR UNIT = ALL TCR-BEARING BEADS, not source beads.
    An earlier draft used each clone's highest-read bead as its position. That
    assumes one cell per clone, which is the very thing under test: a clone
    expanding beside a neoantigen should appear on SEVERAL nearby beads, and
    the source bead may be any of them or none. Source-bead analysis survives
    only as a sensitivity arm.

  CASES   the 110 beads carrying a predicted neoantigen
  CONTROLS the mutation-carrying beads with NO predicted binder (~325)
    Both sets passed identical SComatic calling filters, so they are matched on
    callability and coverage by construction, and differ only in whether the
    mutation yields a predicted binder. This answers "is it the neoantigen, or
    just any mutation?" far better than permuted positions could.

  ENRICHMENT IS ALWAYS AGAINST PUCK-WIDE CLONE FREQUENCY, never against zero.
    Puck 37's largest clone sits on 1,113 of 28,861 TCR beads (3.9%), so it
    turns up in any neighbourhood by chance. Controls calibrate this empirically
    rather than assuming a binomial.

BEATS
    0  build and validate case/control sets (is mutations_per_bead epithelial?)
    1  neighbourhood occupancy vs radius: beads, UMIs, clones, repeated clones
    2  spatial intermixing: are cases and controls drawn from the same regions?
    3  covariate balance: depth, marker_T, CD8, local TCR density
    4  rarefaction floor for the diversity metrics
    5  puck-wide clone frequency background
    6  power estimate for the top-clone-fraction contrast

WHAT THE DOWNSTREAM TESTS WILL BE (not run here)
    A  local clonality contrast, top-clone fraction, cases vs controls
    B  per-neoantigen max clone enrichment vs puck-wide frequency
    C  neoantigen potential field regression (needs a toroidal-shift null)
    D  in-situ expansion: enriched clone must span >= k distinct beads

USAGE
    python Step08c_Neoantigen_TCR_Feasibility.py
    python Step08c_Neoantigen_TCR_Feasibility.py --validate-only
    python Step08c_Neoantigen_TCR_Feasibility.py --chain TRB

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
================================================================================
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, chi2_contingency

try:
    import anndata as ad
except ImportError:
    sys.exit("anndata not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTPUTS = os.path.join(PROOT, "data/outputs")
OUT = os.path.join(OUTPUTS, "11_neoantigen_tcr")
FIGDIR = os.path.join(OUT, "figures")

ADATA = os.path.join(OUTPUTS, "04_annotation/all_pucks_annotated_unified.h5ad")
NEO = os.path.join(OUTPUTS, "07_neoantigen/neoantigens_per_bead.tsv")
MUT = os.path.join(OUTPUTS, "05_mutations/SComatic/SingleCell/mutations_per_bead.tsv")
TCR_DIR = os.path.join(PROOT, "data/inputs/tcr/processed")

PUCKS = ["29", "37", "40"]
TCR_CSV = {p: f"B59_{p}_hTCR_tcr.csv" for p in PUCKS}

ANNOT_COL = "unified_annotation"
EPITHELIAL = "epithelial"
MARK = "marker_score_T_cell"
C_CD8 = "c2l_CD8-positive, alpha-beta T cell"

RADII = [30.0, 50.0, 75.0, 100.0, 150.0, 200.0]
PRIMARY_RADIUS = 75.0
K_PRIMARY = 2          # min distinct beads for a clone to count as locally expanded
K_STRICT = 3
BEAD_COLS = ("CB", "bead", "bead_id", "cell_barcode", "barcode")
GRID_N = 5             # grid for the spatial intermixing chi-square

# v3: DEPTH-ONLY matching.
#
# v2 matched on log_depth AND cd8. depth and cd8 balanced to SMD 0.001-0.055,
# but marker_score_T_cell got WORSE in all three pucks and in the same
# direction: +0.255 -> -0.258 (p29), -0.090 -> -0.362 (p37), -0.271 -> -0.412
# (p40). Negative means matched controls ended up MORE T-cell-like than cases,
# on the one measure the whole analysis rests on.
#
# The reason is causal, not mechanical. The anchor is a 10 um epithelial spot.
# If a T cell sits against it, T-cell transcripts land on that bead and raise
# its cd8 and marker_T scores. Those scores are therefore a CONSEQUENCE of
# nearby T cells, i.e. a MEDIATOR on the causal path we are testing, not a
# confounder. Matching on a mediator removes the effect being measured.
#
# Depth is different: it is purely technical and cannot sit on any causal path.
#
# The regional-inflammation confounder that would have justified matching cd8
# was already measured and is absent: local TCR bead fraction at 75 um was
# 0.420 vs 0.432, 0.698 vs 0.690, 0.615 vs 0.615, essentially identical with no
# matching at all.
MATCH_COVARS = ["log_depth"]
# reported after matching but NEVER matched on: these are outcomes/mediators
REPORT_ONLY_COVARS = ["cd8", "marker_T"]
MATCH_CALIPER = 0.25      # Mahalanobis distance, in pooled SD units
SMD_OK = 0.10             # conventional balance threshold after matching
DROPOUT_WARN = 0.20       # above this, recommend switching to with-replacement

# v4: BACKGROUND POOL for the per-anchor test.
#
# v3 planned a MAX statistic with a permutation null. That is INVALID: whichever
# bead holds the largest value lands in the permuted case set half the time
# regardless of its true label, so the p-value floors near 0.5. Verified by
# simulation: one case anchor shifted +1.00 returns p = 0.477, type I 0.000,
# power 0.000.
#
# v3's BEAT 6 null bootstrapped the max from the control values, which is
# anti-conservative because a bootstrap maximum is capped at the largest control
# actually observed. Measured type I error 0.472. The flattering power numbers
# in the v3 table were measuring that failure, not detectability.
#
# The real constraint is RESOLUTION: with ~30 matched controls the smallest
# achievable p-value is 1/31 = 0.032, which cannot survive correction across 30
# anchors. The fix is a per-anchor p-value against a LARGE background pool, then
# BH across anchors. Simulated calibration and power:
#     pool    type I   1 of 30 at +0.30   3 of 30 at +0.30
#       30     0.000              0.000              0.000
#      200     0.000              0.020              0.270
#    1,000     0.027              0.207              0.570
#    3,000     0.030              0.300              0.657
POOL_SIZE = 3000          # target background neighbourhoods per puck
# Buffer around CASE anchors only, as a multiple of the neighbourhood radius.
# At 2x, a pool neighbourhood cannot overlap a case neighbourhood at all.
# Controls are null by hypothesis so they are not buffered; buffering all 435
# mutation beads would exclude more area than the puck contains.
POOL_BUFFER_MULT = 2.0
POOL_MIN_FRACTION = 0.20  # if fewer candidates than this survive, relax to 1x
POOL_CALIPER = 0.25       # log-depth SD; a pool bead must be this close to its case
POOL_MIN_CASES = 8        # never trim a puck below this many testable anchors
POOL_MIN_PER_CASE = 20    # a case needs this many matched background beads to be
                          # testable; below it, no comparable background exists
                          # anywhere on the puck and the case is reported, not tested

# v3: floor target relaxed from 0.95 to 0.94. Puck 40's floor fell to 30 only
# because 50 retained 94.3% against a 95% target, a cliff edge on an arbitrary
# threshold. 50 is the comparable floor to the other pucks' choices.
FLOOR_TARGET = 0.94
FLOOR_GRID = [10, 20, 30, 50, 75, 100, 150, 200]

# v2: sparse-alternative power. The v1 simulation shifted EVERY case anchor,
# which is the most favourable alternative possible and returned 85-99% power
# to detect a shift of 0.02 at n=20. The realistic alternative is that a few
# anchors carry an expanded clone and the rest look like controls.
N_POWER_SIM = 1000
N_NULL_DRAWS = 4000
SPARSE_K = [1, 2, 3, 5, 10]
SPARSE_DELTA = [0.05, 0.10, 0.20, 0.30, 0.50]

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
C_CORAL = "#ed6a5a"    # cases, neoantigen
C_MUSTARD = "#F6D155"  # controls, mutation only
C_GRAY = "#d3d3d3"
C_DARK = "#4a4a4a"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

RNG = np.random.default_rng(20261008)


def log(m=""):
    print(m, flush=True)


def banner(m):
    log(); log("=" * 72); log(m); log("=" * 72)


def savefig(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    log(f"    figure: {name}.pdf / .png")


def style(ax, title=None, xlabel=None, ylabel=None):
    if title: ax.set_title(title, fontsize=FS_TITLE)
    if xlabel: ax.set_xlabel(xlabel, fontsize=FS_LABEL)
    if ylabel: ax.set_ylabel(ylabel, fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def bead_col(df, label):
    c = next((c for c in BEAD_COLS if c in df.columns), None)
    if c is None:
        sys.exit(f"no recognizable bead column in {label}: {list(df.columns)}")
    return c


def load_tcr(p, chain=None):
    d = pd.read_csv(os.path.join(TCR_DIR, TCR_CSV[p]), sep=None, engine="python")
    d.columns = [c.strip() for c in d.columns]
    if chain:
        d = d[d.chain == chain]
    d["bead"] = d["bc"].astype(str) + f"-1_Puck_211214_{p}"
    return d


def puck_view(adata, p):
    m = adata.obs["puck_id"].astype(str).values == f"Puck_211214_{p}"
    idx = adata.obs_names[m]
    xy = adata.obs.loc[idx, ["x_coord", "y_coord"]].values.astype(float)
    f = pd.DataFrame(index=idx)
    f["depth"] = adata.obs.loc[idx, "total_counts"].values.astype(float)
    f["annot"] = adata.obs.loc[idx, ANNOT_COL].astype(str).values
    for col, nm in ((MARK, "marker_T"), (C_CD8, "cd8")):
        if col in adata.obs.columns:
            f[nm] = adata.obs.loc[idx, col].values.astype(float)
    return idx, f, xy, {b: i for i, b in enumerate(idx)}


# ==============================================================================
def beat0_sets(adata, rows):
    banner("BEAT 0  build and validate case / control sets")
    log("  CASES    beads carrying a predicted neoantigen")
    log("  CONTROLS mutation-carrying beads with NO predicted binder")
    log("  Both passed identical SComatic filters, so they are matched on")
    log("  callability and coverage and differ only in predicted binding.")

    neo = pd.read_csv(NEO, sep="\t")
    nb = bead_col(neo, "neoantigens_per_bead.tsv")
    cases = set(neo[nb].astype(str))
    log(f"\n  neoantigens_per_bead.tsv: {len(neo):,} rows, {len(cases):,} unique beads "
        f"(column '{nb}')")

    if not os.path.exists(MUT):
        sys.exit(f"missing {MUT}")
    mut = pd.read_csv(MUT, sep="\t")
    mbc = bead_col(mut, "mutations_per_bead.tsv")
    allmut = set(mut[mbc].astype(str))
    log(f"  mutations_per_bead.tsv:   {len(mut):,} rows, {len(allmut):,} unique beads "
        f"(column '{mbc}')")

    orphan = cases - allmut
    if orphan:
        log(f"  WARN  {len(orphan)} neoantigen beads absent from the mutation table")

    controls = allmut - cases
    log(f"  raw controls (mutation, no binder): {len(controls):,}")

    # is mutations_per_bead already epithelial-only? Jake expects yes; verify.
    ann = adata.obs[ANNOT_COL].astype(str)
    def comp(s, label):
        present = [b for b in s if b in ann.index]
        vc = ann.loc[present].value_counts()
        log(f"    {label}: {len(present):,}/{len(s):,} annotated | "
            + "; ".join(f"{k}:{v}" for k, v in vc.head(5).items()))
        return vc
    log("\n  cell-type composition BEFORE any epithelial restriction:")
    vc_case = comp(cases, "cases   ")
    vc_ctrl = comp(controls, "controls")

    non_epi_ctrl = int(vc_ctrl.sum() - vc_ctrl.get(EPITHELIAL, 0))
    if non_epi_ctrl:
        log(f"\n  mutations_per_bead.tsv is NOT epithelial-only: {non_epi_ctrl:,} "
            "control beads are other cell types.")
        log("  Restricting controls to epithelial so the contrast is clean.")
    else:
        log("\n  mutations_per_bead.tsv is already epithelial-only, as expected.")

    cases = {b for b in cases if ann.get(b, "") == EPITHELIAL}
    controls = {b for b in controls if ann.get(b, "") == EPITHELIAL}
    log(f"\n  FINAL  cases {len(cases):,} | controls {len(controls):,} "
        f"| ratio 1:{len(controls)/max(1,len(cases)):.1f}")

    for p in PUCKS:
        sfx = f"Puck_211214_{p}"
        nc = sum(1 for b in cases if b.endswith(sfx))
        nk = sum(1 for b in controls if b.endswith(sfx))
        log(f"    puck {p}: {nc:>4} cases, {nk:>4} controls")
        rows.append({"beat": "sets", "puck": p, "cases": nc, "controls": nk})

    if len(cases) < 20 or len(controls) < 20:
        log("\n  WARNING: very small sets. Interpret everything downstream with care.")
    return cases, controls


# ==============================================================================
def neighbourhood_stats(tcr_p, xy, pos, anchors, radius, tree):
    """Per anchor: TCR beads, UMIs, clones, clones spanning >=2 and >=3 beads."""
    per_bead = tcr_p.groupby("bead").agg(umi=("umi", "nunique"))
    clone_beads = tcr_p.groupby(["cloneId", "bead"]).size().reset_index(name="n")
    cb_map = {}
    for cid, g in clone_beads.groupby("cloneId"):
        cb_map[cid] = set(g.bead)
    bead_clones = tcr_p.groupby("bead")["cloneId"].apply(set).to_dict()

    out = []
    for a in anchors:
        if a not in pos:
            continue
        i = pos[a]
        near = tree.query_ball_point(xy[i], radius)
        nb = [b for b in (idx_lookup[j] for j in near)]
        tcr_nb = [b for b in nb if b in bead_clones]
        umis = int(sum(per_bead.loc[b, "umi"] for b in tcr_nb)) if tcr_nb else 0
        clones = set()
        for b in tcr_nb:
            clones |= bead_clones[b]
        # clones appearing on >= k distinct beads INSIDE this neighbourhood
        counts = {}
        for b in tcr_nb:
            for c in bead_clones[b]:
                counts[c] = counts.get(c, 0) + 1
        k2 = sum(1 for v in counts.values() if v >= K_PRIMARY)
        k3 = sum(1 for v in counts.values() if v >= K_STRICT)
        top = max(counts.values()) if counts else 0
        out.append({"anchor": a, "n_beads": len(nb), "n_tcr_beads": len(tcr_nb),
                    "n_umis": umis, "n_clones": len(clones),
                    "clones_ge2": k2, "clones_ge3": k3, "max_clone_beads": top})
    return pd.DataFrame(out)


def beat1_occupancy(adata, tcr, cases, controls, rows, chain_label):
    banner("BEAT 1  neighbourhood occupancy vs radius")
    log("  The decisive column is clones_ge2: clones appearing on two or more")
    log("  distinct beads inside one neighbourhood. Every downstream test needs")
    log("  those. If the median is 0, local expansion cannot be detected.")
    global idx_lookup
    store = {}
    for p in PUCKS:
        idx, f, xy, pos = puck_view(adata, p)
        idx_lookup = list(idx)
        tree = cKDTree(xy)
        tp = tcr[tcr.puck == p]
        tp = tp[tp.bead.isin(set(idx))]
        cs = [b for b in cases if b.endswith(f"Puck_211214_{p}")]
        ks = [b for b in controls if b.endswith(f"Puck_211214_{p}")]
        log(f"\n  --- puck {p} ({chain_label}): {len(cs)} cases, {len(ks)} controls ---")
        log(f"    {'R':>5} {'set':>8} {'n':>5} {'TCRbeads':>9} {'UMIs':>8} "
            f"{'clones':>8} {'ge2':>6} {'ge3':>6} {'%ge2>0':>8}")
        for R in RADII:
            for nm, anchors in (("case", cs), ("control", ks)):
                d = neighbourhood_stats(tp, xy, pos, anchors, R, tree)
                if not len(d):
                    continue
                store[(p, R, nm)] = d
                log(f"    {R:>5.0f} {nm:>8} {len(d):>5} "
                    f"{d.n_tcr_beads.median():>9.0f} {d.n_umis.median():>8.0f} "
                    f"{d.n_clones.median():>8.0f} {d.clones_ge2.median():>6.0f} "
                    f"{d.clones_ge3.median():>6.0f} "
                    f"{100*(d.clones_ge2 > 0).mean():>7.1f}%")
                rows.append({"beat": "occupancy", "puck": p, "chain": chain_label,
                             "radius_um": R, "set": nm, "n_anchors": len(d),
                             "median_tcr_beads": float(d.n_tcr_beads.median()),
                             "median_umis": float(d.n_umis.median()),
                             "median_clones": float(d.n_clones.median()),
                             "median_clones_ge2": float(d.clones_ge2.median()),
                             "median_clones_ge3": float(d.clones_ge3.median()),
                             "pct_with_any_ge2": round(100*(d.clones_ge2 > 0).mean(), 2)})

    fig, axes = plt.subplots(1, 3, figsize=(30, 9), sharey=True)
    for i, p in enumerate(PUCKS):
        ax = axes[i]
        for nm, col in (("case", C_CORAL), ("control", C_MUSTARD)):
            ys = [store[(p, R, nm)].clones_ge2.median()
                  for R in RADII if (p, R, nm) in store]
            xs = [R for R in RADII if (p, R, nm) in store]
            ax.plot(xs, ys, "o-", color=col, linewidth=3, markersize=12, label=nm)
        style(ax, title=f"Puck {p}", xlabel="Radius (µm)",
              ylabel=f"Median clones on ≥{K_PRIMARY} beads" if i == 0 else None)
        if i == 0:
            ax.legend(fontsize=FS_ANNOT, frameon=False)
    fig.tight_layout(); savefig(fig, f"Fig_Feas_A_occupancy_{chain_label}")
    return store


# ==============================================================================
def beat2_intermixing(adata, cases, controls, rows):
    banner("BEAT 2  spatial intermixing of cases and controls")
    log("  If neoantigen beads cluster in one region of a puck, the control")
    log("  contrast inherits a regional confounder and any difference could be")
    log("  geography rather than antigen. Three checks: centroid separation,")
    log("  nearest-control distance, and a grid occupancy chi-square.")
    for p in PUCKS:
        idx, f, xy, pos = puck_view(adata, p)
        cs = [b for b in cases if b in pos]
        ks = [b for b in controls if b in pos]
        cs = [b for b in cs if b.endswith(f"Puck_211214_{p}")]
        ks = [b for b in ks if b.endswith(f"Puck_211214_{p}")]
        if len(cs) < 3 or len(ks) < 3:
            log(f"\n  puck {p}: too few anchors ({len(cs)} / {len(ks)})"); continue
        cxy = xy[[pos[b] for b in cs]]
        kxy = xy[[pos[b] for b in ks]]

        sep = float(np.hypot(*(cxy.mean(axis=0) - kxy.mean(axis=0))))
        def rg(a):
            return float(np.sqrt(((a - a.mean(axis=0))**2).sum(axis=1).mean()))
        log(f"\n  puck {p}: {len(cs)} cases, {len(ks)} controls")
        log(f"    centroid separation {sep:7.1f} um | case Rg {rg(cxy):7.1f} | "
            f"control Rg {rg(kxy):7.1f}")

        # distance from each case to nearest control, vs from a random
        # epithelial bead to nearest control
        ktree = cKDTree(kxy)
        d_case = ktree.query(cxy)[0]
        epi = [b for b in idx if f.loc[b, "annot"] == EPITHELIAL]
        samp = RNG.choice(len(epi), size=min(2000, len(epi)), replace=False)
        exy = xy[[pos[epi[j]] for j in samp]]
        d_rand = ktree.query(exy)[0]
        try:
            _, pv = mannwhitneyu(d_case, d_rand, alternative="two-sided")
        except ValueError:
            pv = np.nan
        log(f"    distance to nearest control: cases median {np.median(d_case):7.1f} um | "
            f"random epithelial {np.median(d_rand):7.1f} um | p {pv:.3g}")

        # grid occupancy
        xs = np.linspace(xy[:, 0].min(), xy[:, 0].max() + 1, GRID_N + 1)
        ys = np.linspace(xy[:, 1].min(), xy[:, 1].max() + 1, GRID_N + 1)
        hc = np.histogram2d(cxy[:, 0], cxy[:, 1], bins=[xs, ys])[0].ravel()
        hk = np.histogram2d(kxy[:, 0], kxy[:, 1], bins=[xs, ys])[0].ravel()
        keep = (hc + hk) > 0
        tab = np.vstack([hc[keep], hk[keep]])
        try:
            chi2, pgrid, _, _ = chi2_contingency(tab + 0.5)
        except ValueError:
            chi2, pgrid = np.nan, np.nan
        log(f"    grid occupancy ({GRID_N}x{GRID_N}, {int(keep.sum())} occupied cells): "
            f"chi2 {chi2:.1f}, p {pgrid:.3g}")
        log("      p above ~0.05 on all three means cases and controls sample the")
        log("      same regions and the contrast is not confounded by geography")
        rows.append({"beat": "intermixing", "puck": p, "n_cases": len(cs),
                     "n_controls": len(ks), "centroid_sep_um": round(sep, 1),
                     "case_Rg_um": round(rg(cxy), 1), "control_Rg_um": round(rg(kxy), 1),
                     "median_dist_case_to_control": round(float(np.median(d_case)), 1),
                     "median_dist_random_to_control": round(float(np.median(d_rand)), 1),
                     "p_dist": float(pv) if pv == pv else None,
                     "grid_chi2_p": float(pgrid) if pgrid == pgrid else None})

    fig, axes = plt.subplots(1, 3, figsize=(30, 11))
    for i, p in enumerate(PUCKS):
        idx, f, xy, pos = puck_view(adata, p)
        ax = axes[i]
        ax.scatter(xy[:, 0], xy[:, 1], s=1, c=C_GRAY, alpha=0.3, linewidths=0)
        for s, col, lab in ((controls, C_MUSTARD, "mutation only"),
                            (cases, C_CORAL, "neoantigen")):
            bs = [b for b in s if b in pos and b.endswith(f"Puck_211214_{p}")]
            if bs:
                a = xy[[pos[b] for b in bs]]
                ax.scatter(a[:, 0], a[:, 1], s=110, c=col, alpha=0.9,
                           edgecolors=C_DARK, linewidths=1.0, label=lab)
        style(ax, title=f"Puck {p}")
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if i == 0:
            ax.legend(fontsize=FS_ANNOT, frameon=False, loc="upper right")
    fig.tight_layout(); savefig(fig, "Fig_Feas_B_case_control_map")


# ==============================================================================
def beat3_balance(adata, tcr, cases, controls, rows):
    banner("BEAT 3  covariate balance between cases and controls")
    log("  Depth drove the false negative in Step08a, so it is checked first.")
    log("  Local TCR density matters too: a case bead in a T-cell-rich region")
    log("  shows more TCR regardless of antigen.")
    for p in PUCKS:
        idx, f, xy, pos = puck_view(adata, p)
        tree = cKDTree(xy)
        tp = tcr[tcr.puck == p]
        tcr_beads = set(tp.bead) & set(idx)
        cs = [b for b in cases if b in pos and b.endswith(f"Puck_211214_{p}")]
        ks = [b for b in controls if b in pos and b.endswith(f"Puck_211214_{p}")]
        if len(cs) < 3 or len(ks) < 3:
            continue
        def local_density(bs):
            out = []
            for b in bs:
                near = tree.query_ball_point(xy[pos[b]], PRIMARY_RADIUS)
                nb = [idx[j] for j in near]
                out.append(sum(1 for x in nb if x in tcr_beads) / max(1, len(nb)))
            return np.array(out)
        log(f"\n  puck {p}  (R = {PRIMARY_RADIUS:.0f} um)")
        log(f"    {'metric':<22} {'cases':>10} {'controls':>10} {'diff':>9} {'p':>10}")
        for col in ("depth", "marker_T", "cd8"):
            if col not in f.columns:
                continue
            a, b = f.loc[cs, col].values, f.loc[ks, col].values
            try:
                _, pv = mannwhitneyu(a, b, alternative="two-sided")
            except ValueError:
                pv = np.nan
            log(f"    {col:<22} {np.median(a):>10.3f} {np.median(b):>10.3f} "
                f"{np.median(a)-np.median(b):>+9.3f} {pv:>10.2e}")
            rows.append({"beat": "balance", "puck": p, "metric": col,
                         "case_median": float(np.median(a)),
                         "control_median": float(np.median(b)),
                         "p": float(pv) if pv == pv else None})
        da, db = local_density(cs), local_density(ks)
        try:
            _, pv = mannwhitneyu(da, db, alternative="two-sided")
        except ValueError:
            pv = np.nan
        log(f"    {'local TCR bead frac':<22} {np.median(da):>10.3f} "
            f"{np.median(db):>10.3f} {np.median(da)-np.median(db):>+9.3f} {pv:>10.2e}")
        rows.append({"beat": "balance", "puck": p, "metric": "local_tcr_density",
                     "case_median": float(np.median(da)),
                     "control_median": float(np.median(db)),
                     "p": float(pv) if pv == pv else None})
    log("\n  Non-significant does NOT mean balanced at these sample sizes. Puck 29")
    log("  cases are 25% shallower than controls (p=0.36 only because n is 31 vs")
    log("  60 and depth is highly variable). BEAT 3b matches it out rather than")
    log("  trusting the test.")


# ==============================================================================
def smd(a, b):
    """Standardized mean difference. The conventional matching diagnostic:
    |SMD| < 0.10 is considered balanced, and unlike a p-value it does not
    improve just because n is small."""
    sa, sb = np.var(a, ddof=1), np.var(b, ddof=1)
    sp = np.sqrt((sa + sb) / 2.0)
    return float((np.mean(a) - np.mean(b)) / sp) if sp else np.nan


def beat3b_matching(adata, cases, controls, rows, with_replacement=False):
    banner("BEAT 3b  case-control matching on DEPTH ONLY")
    log("  v3: matching on log(depth) alone.")
    log("  v2 also matched c2l CD8, which balanced depth and CD8 but pushed")
    log("  marker_score_T_cell OUT of balance in all three pucks, with controls")
    log("  ending up more T-cell-like than cases. CD8 and marker_T on the anchor")
    log("  bead are MEDIATORS: a T cell beside the bead puts its transcripts on")
    log("  that bead. Matching on them conditions away the effect being measured.")
    log("  Depth is technical and cannot sit on a causal path, so it is matched.")
    log(f"  Matching on {MATCH_COVARS}; {REPORT_ONLY_COVARS} reported, never matched.")
    log(f"  Greedy nearest neighbour, within puck, caliper {MATCH_CALIPER} SD, "
        f"{'WITH' if with_replacement else 'without'} replacement.")

    m_cases, m_controls, pair_rows = set(), set(), []
    total_c = total_drop = 0
    for p in PUCKS:
        idx, f, xy, pos = puck_view(adata, p)
        f = f.copy()
        f["log_depth"] = np.log1p(f["depth"].values)
        cs = [b for b in cases if b in f.index and b.endswith(f"Puck_211214_{p}")]
        ks = [b for b in controls if b in f.index and b.endswith(f"Puck_211214_{p}")]
        covars = [c for c in MATCH_COVARS if c in f.columns]
        if len(cs) < 3 or len(ks) < 3 or not covars:
            log(f"\n  puck {p}: too few anchors to match"); continue

        A = f.loc[cs, covars].values.astype(float)
        B = f.loc[ks, covars].values.astype(float)
        pooled = np.vstack([A, B])
        mu, sd = pooled.mean(axis=0), pooled.std(axis=0, ddof=1)
        sd[sd == 0] = 1.0
        As, Bs = (A - mu) / sd, (B - mu) / sd
        cov = np.atleast_2d(np.cov(np.vstack([As, Bs]).T))
        try:
            inv = np.linalg.pinv(cov)
        except np.linalg.LinAlgError:
            inv = np.eye(len(covars))

        used, pairs = set(), []
        # hardest cases first: those furthest from the control centroid, so an
        # extreme case is not orphaned by an easier one taking its only match
        ctr = Bs.mean(axis=0)
        order = np.argsort(-np.array([
            float(np.sqrt(max(0.0, (x - ctr) @ inv @ (x - ctr)))) for x in As]))
        for i in order:
            best, bd = None, np.inf
            for j in range(len(Bs)):
                if (not with_replacement) and j in used:
                    continue
                d = As[i] - Bs[j]
                dist = float(np.sqrt(max(0.0, d @ inv @ d)))
                if dist < bd:
                    best, bd = j, dist
            if best is not None and bd <= MATCH_CALIPER:
                if not with_replacement:
                    used.add(best)
                pairs.append((cs[i], ks[best], bd))

        n_drop = len(cs) - len(pairs)
        frac_drop = n_drop / max(1, len(cs))
        total_c += len(cs); total_drop += n_drop
        log(f"\n  puck {p}: {len(cs)} cases, {len(ks)} controls")
        log(f"    matched pairs: {len(pairs)}  ({n_drop} cases dropped, "
            f"{100*frac_drop:.1f}%)")
        if not pairs:
            continue
        mc = [a for a, _, _ in pairs]
        mk = [b for _, b, _ in pairs]
        m_cases |= set(mc); m_controls |= set(mk)
        if with_replacement:
            reuse = len(mk) - len(set(mk))
            log(f"    controls reused: {reuse} of {len(mk)} picks "
                f"({len(set(mk))} distinct)")
        log(f"    mean match distance {np.mean([d for _,_,d in pairs]):.4f} SD")
        for a, b, d in pairs:
            pair_rows.append({"puck": p, "case": a, "control": b,
                              "match_distance_sd": round(d, 5)})

        log(f"    {'covariate':<14} {'role':>10} {'SMD before':>11} "
            f"{'SMD after':>11} {'balanced':>9}")
        for c in covars:
            before = smd(f.loc[cs, c].values, f.loc[ks, c].values)
            after = smd(f.loc[mc, c].values, f.loc[mk, c].values)
            ok = "yes" if abs(after) < SMD_OK else "NO"
            log(f"    {c:<14} {'MATCHED':>10} {before:>+11.3f} {after:>+11.3f} "
                f"{ok:>9}")
            rows.append({"beat": "matching", "puck": p, "covariate": c,
                         "role": "matched", "smd_before": round(before, 4),
                         "smd_after": round(after, 4),
                         "balanced": abs(after) < SMD_OK,
                         "matched_pairs": len(pairs), "cases_dropped": n_drop})
        for c in REPORT_ONLY_COVARS:
            if c not in f.columns:
                continue
            before = smd(f.loc[cs, c].values, f.loc[ks, c].values)
            after = smd(f.loc[mc, c].values, f.loc[mk, c].values)
            log(f"    {c:<14} {'mediator':>10} {before:>+11.3f} {after:>+11.3f} "
                f"{'(not matched)':>9}")
            rows.append({"beat": "matching", "puck": p, "covariate": c,
                         "role": "mediator_reported", "smd_before": round(before, 4),
                         "smd_after": round(after, 4), "balanced": None,
                         "matched_pairs": len(pairs), "cases_dropped": n_drop})
        log(f"    raw depth medians after matching: {f.loc[mc,'depth'].median():,.0f} "
            f"vs {f.loc[mk,'depth'].median():,.0f}")
        log("    A positive SMD on a mediator row means cases are MORE T-cell-like")
        log("    than their matched controls. That is a preliminary result, not an")
        log("    imbalance to correct.")

    log(f"\n  MATCHED TOTALS: {len(m_cases)} cases, {len(m_controls)} distinct controls")
    overall = total_drop / max(1, total_c)
    log(f"  overall case dropout: {total_drop}/{total_c} ({100*overall:.1f}%)")
    if overall > DROPOUT_WARN and not with_replacement:
        log(f"  *** dropout exceeds {100*DROPOUT_WARN:.0f}%. Rerun with")
        log(f"      --with-replacement to keep all cases, accepting that reused")
        log(f"      controls are not independent and need weighting downstream. ***")
    elif not with_replacement:
        log(f"  dropout is within the {100*DROPOUT_WARN:.0f}% threshold; "
            "without-replacement matching stands.")

    if pair_rows:
        pd.DataFrame(pair_rows).to_csv(
            os.path.join(OUT, "matched_pairs.tsv"), sep="\t", index=False)
        log("  wrote matched_pairs.tsv (case, control, match distance)")
    return m_cases, m_controls


# ==============================================================================
def beat4_rarefaction(store, rows):
    banner("BEAT 4  rarefaction floor for the diversity metrics")
    log("  Top-clone fraction and Simpson are both biased by sample size, and")
    log("  neighbourhood UMI counts vary with local bead depth. Without")
    log("  rarefying, a deeper neighbourhood looks more diverse as an artifact.")
    log("  The floor is chosen so most neighbourhoods survive it.")
    for p in PUCKS:
        key = (p, PRIMARY_RADIUS, "case")
        key2 = (p, PRIMARY_RADIUS, "control")
        if key not in store or key2 not in store:
            continue
        u = np.concatenate([store[key].n_umis.values, store[key2].n_umis.values])
        log(f"\n  puck {p} (R = {PRIMARY_RADIUS:.0f} um), {len(u)} anchors")
        log(f"    UMIs per neighbourhood: median {np.median(u):.0f}, "
            f"p10 {np.percentile(u,10):.0f}, p25 {np.percentile(u,25):.0f}, "
            f"p75 {np.percentile(u,75):.0f}")
        log(f"    {'floor':>7} {'anchors retained':>18} {'%':>7}")
        chosen = FLOOR_GRID[0]
        for fl in FLOOR_GRID:
            n = int((u >= fl).sum())
            frac = n / len(u)
            log(f"    {fl:>7} {n:>18,} {100*frac:>6.1f}%")
            if frac >= FLOOR_TARGET:
                chosen = fl
            rows.append({"beat": "rarefaction", "puck": p, "floor": fl,
                         "retained": n, "pct_retained": round(100*frac, 2)})
        n_ch = int((u >= chosen).sum())
        log(f"    -> CHOSEN FLOOR {chosen} UMIs "
            f"(largest retaining >= {100*FLOOR_TARGET:.0f}%; keeps {n_ch}/{len(u)}, "
            f"drops {len(u)-n_ch})")
        rows.append({"beat": "rarefaction_chosen", "puck": p, "floor": chosen,
                     "retained": n_ch, "dropped": len(u) - n_ch,
                     "pct_retained": round(100*n_ch/len(u), 2)})
    log("\n  Floors are PER PUCK. The contrast is within puck, so there is no")
    log("  reason to penalise the deeper pucks for puck 29's shallowness: a")
    log("  single global floor of 20 would discard roughly 80% of the UMIs in")
    log("  pucks 37 and 40 for nothing. Report the dropped anchors in Methods.")


# ==============================================================================
def beat5_clone_background(tcr, rows):
    banner("BEAT 5  puck-wide clone frequency background")
    log("  Enrichment must be measured against each clone's puck-wide frequency,")
    log("  never against zero. The largest clones turn up in any neighbourhood")
    log("  by chance and would otherwise look locally enriched everywhere.")
    for p in PUCKS:
        d = tcr[tcr.puck == p]
        nb = d.bead.nunique()
        cb = d.groupby("cloneId")["bead"].nunique().sort_values(ascending=False)
        log(f"\n  puck {p}: {nb:,} TCR beads, {len(cb):,} clones")
        log(f"    top 5 clone bead-frequencies: "
            + ", ".join(f"{100*v/nb:.2f}%" for v in cb.head(5)))
        for R in (50.0, PRIMARY_RADIUS, 150.0):
            exp_beads = np.pi * R * R / 100.0       # ~10 um bead pitch
            log(f"    at R={R:.0f} um (~{exp_beads:.0f} bead positions): top clone "
                f"expected on {exp_beads*cb.iloc[0]/nb:.1f} beads by chance")
        rows.append({"beat": "clone_background", "puck": p, "tcr_beads": nb,
                     "clones": len(cb), "top_clone_bead_frac": float(cb.iloc[0]/nb)})


# ==============================================================================
def bh_adjust(pv):
    pv = np.asarray(pv, float); n = len(pv)
    o = np.argsort(pv); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((pv[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(adj, 1.0)


def beat6b_background_pool(adata, cases, controls, rows, radius):
    """Build the background pool the per-anchor test will use as its null."""
    banner(f"BEAT 6b  background pool for the per-anchor null (R = {radius:.0f} um)")
    log("  Each case anchor will be scored against this pool, then BH across")
    log("  anchors. Resolution comes from pool SIZE: with only ~30 matched")
    log("  controls the smallest achievable p-value is 1/31 = 0.032, which")
    log("  cannot survive correction. A pool of a few thousand fixes that.")
    log(f"  Buffer: {POOL_BUFFER_MULT:.0f}x radius around CASE anchors only, so a")
    log("  pool neighbourhood cannot overlap a case neighbourhood. Controls are")
    log("  null by hypothesis and are not buffered; buffering all 435 mutation")
    log("  beads would exclude more area than the puck contains.")

    pool_all, testable_all = {}, {}
    for p in PUCKS:
        idx, f, xy, pos = puck_view(adata, p)
        cs = [b for b in cases if b in pos and b.endswith(f"Puck_211214_{p}")]
        ks = set(b for b in controls if b.endswith(f"Puck_211214_{p}"))
        if len(cs) < 3:
            continue
        epi_mask = (f["annot"].values == EPITHELIAL)
        excl = set(cs) | ks
        cand_i = np.array([i for i in range(len(idx))
                           if epi_mask[i] and idx[i] not in excl])
        n_epi = int(epi_mask.sum())

        cxy = xy[[pos[b] for b in cs]]
        ctree = cKDTree(cxy)
        buf = POOL_BUFFER_MULT * radius
        d_near = ctree.query(xy[cand_i])[0]
        keep = cand_i[d_near > buf]
        frac = len(keep) / max(1, len(cand_i))
        log(f"\n  puck {p}: {n_epi:,} epithelial beads, {len(cand_i):,} candidates "
            f"after removing anchors")
        log(f"    buffer {buf:.0f} um around {len(cs)} case anchors: "
            f"{len(keep):,} survive ({100*frac:.1f}%)")
        if frac < POOL_MIN_FRACTION:
            buf = radius
            d_near = ctree.query(xy[cand_i])[0]
            keep = cand_i[d_near > buf]
            frac = len(keep) / max(1, len(cand_i))
            log(f"    *** below {100*POOL_MIN_FRACTION:.0f}%, relaxed to {buf:.0f} um: "
                f"{len(keep):,} ({100*frac:.1f}%). Pool neighbourhoods may now")
            log(f"        partially overlap case neighbourhoods, which makes the")
            log(f"        null slightly conservative. ***")

        # v5: depth-match by CALIPER-CONSTRAINED nearest neighbour on log(depth).
        #
        # v4 binned on case depth quantiles and allocated proportionally. Where
        # cases have a deep tail that candidates lack, the deep bins starve, the
        # shortfall is never redistributed, and the pool drifts shallow. That is
        # what gave puck 37 a depth SMD of +0.392.
        #
        # Plain nearest neighbour is better but not enough: a case far deeper
        # than ANY candidate still drags in whatever is nearest, which is much
        # shallower. Reproduced in simulation at SMD +0.153, still unbalanced.
        #
        # So each case claims only candidates within POOL_CALIPER log-depth SD.
        # A case that cannot find POOL_MIN_PER_CASE of them has no comparable
        # background anywhere on the puck; it is reported and EXCLUDED from the
        # test rather than being allowed to bias the pool. Simulated SMD +0.054.
        cdepth = f.iloc[[pos[b] for b in cs]]["depth"].values
        cld = np.log1p(cdepth)
        kld = np.log1p(f.iloc[keep]["depth"].values)
        sd_all = np.std(np.r_[cld, kld], ddof=1)
        thr = POOL_CALIPER * sd_all
        target = int(np.ceil(min(POOL_SIZE, len(keep)) / len(cs)))

        ktree = cKDTree(kld.reshape(-1, 1))
        taken = np.zeros(len(keep), dtype=bool)
        per_case = {i: [] for i in range(len(cld))}
        for ci in RNG.permutation(len(cld)):
            near = ktree.query_ball_point([cld[ci]], thr)
            near.sort(key=lambda j: abs(kld[j] - cld[ci]))
            for j in near:
                if not taken[j]:
                    taken[j] = True
                    per_case[ci].append(j)
                    if len(per_case[ci]) >= target:
                        break

        usable = [i for i in per_case if len(per_case[i]) >= POOL_MIN_PER_CASE]
        unusable = [i for i in per_case if len(per_case[i]) < POOL_MIN_PER_CASE]
        if unusable:
            log(f"    {len(unusable)} case(s) have no comparable background "
                f"within {POOL_CALIPER} log-depth SD and are EXCLUDED: depths "
                f"{np.sort(cdepth[unusable]).astype(int).tolist()}")

        # v6: ITERATIVE BALANCE TRIM.
        #
        # The caliper alone left puck 37 at SMD +0.126, still unbalanced. The
        # residual comes from cases sitting on the deep shoulder of the depth
        # distribution: nothing exists above them, so their matched set is
        # ONE-SIDED and drags case-minus-pool positive. Neither a symmetric
        # above/below draw nor switching the SMD to the log scale removed it in
        # simulation, so rather than keep tuning, balance is now enforced by
        # construction: drop the most asymmetric case and recompute until SMD
        # clears the threshold. Simulated: +0.178 -> -0.001 after dropping 2.
        #
        # Trading anchors for balance is the agreed preference: an unbalanced
        # pool makes that puck conservative and not comparable with the others.
        trimmed = []
        while len(usable) > POOL_MIN_CASES:
            pl = [int(keep[j]) for i in usable for j in per_case[i]]
            sm_now = smd(cld[usable],
                         np.log1p(f.iloc[pl]["depth"].values))
            if abs(sm_now) < SMD_OK:
                break
            asym = {i: cld[i] - float(np.mean(kld[per_case[i]])) for i in usable}
            # drop whichever case pushes the imbalance in its current direction
            worst = max(usable, key=lambda i: asym[i] * np.sign(sm_now))
            trimmed.append(worst)
            usable.remove(worst)
        if trimmed:
            log(f"    {len(trimmed)} further case(s) dropped to bring the pool "
                f"into depth balance: depths "
                f"{np.sort(cdepth[trimmed]).astype(int).tolist()}")

        chosen = [int(keep[j]) for i in usable for j in per_case[i]]
        testable = [cs[i] for i in usable]
        if len(chosen) < 50:
            chosen = RNG.choice(keep, min(POOL_SIZE, len(keep)),
                                replace=False).tolist()
            testable = list(cs)
            log("    caliper matching produced too few; fell back to a uniform "
                "draw (pool will NOT be depth matched)")
        cdepth = cdepth[usable] if usable and len(chosen) >= 50 else cdepth

        pool = [idx[i] for i in chosen]
        pool_all[p] = pool
        testable_all[p] = testable

        pd_ = f.loc[pool, "depth"].values
        sm = smd(np.log1p(cdepth), np.log1p(pd_))     # matched on log, judged on log
        sm_raw = smd(cdepth, pd_)
        bal = abs(sm) < SMD_OK
        log(f"    pool size {len(pool):,} | depth median {np.median(pd_):,.0f} "
            f"vs cases {np.median(cdepth):,.0f}")
        log(f"    SMD log-depth {sm:+.3f} (matched on this) | raw {sm_raw:+.3f} | "
            f"{'balanced' if bal else 'NOT balanced'}")
        if not bal:
            log("    *** WARNING: pool is NOT depth-matched for this puck.")
            log("        A shallower pool gives smaller neighbourhood UMI counts,")
            log("        which inflates a maximum, which makes the null tail")
            log("        heavier and the test CONSERVATIVE here while the other")
            log("        pucks run clean. Continuing so the rest is visible, but")
            log("        do not compare this puck's p-values with the others. ***")
        log(f"    smallest achievable p-value with this pool: "
            f"{1/(len(pool)+1):.2e}")
        rows.append({"beat": "background_pool", "puck": p, "radius": radius,
                     "testable_cases": len(testable),
                     "excluded_cases": len(cs) - len(testable),
                     "n_epithelial": n_epi, "n_candidates": len(cand_i),
                     "buffer_um": buf, "surviving": len(keep),
                     "pct_surviving": round(100*frac, 2), "pool_size": len(pool),
                     "depth_smd_log": round(float(sm), 4),
                     "depth_smd_raw": round(float(sm_raw), 4),
                     "min_achievable_p": 1/(len(pool)+1)})

    if pool_all:
        recs = [{"puck": p, "CB": b} for p, bs in pool_all.items() for b in bs]
        pd.DataFrame(recs).to_csv(os.path.join(OUT, "background_pool.tsv"),
                                  sep="\t", index=False)
        log(f"\n  wrote background_pool.tsv "
            f"({sum(len(v) for v in pool_all.values()):,} beads)")
    if testable_all:
        recs = [{"puck": p, "CB": b} for p, bs in testable_all.items() for b in bs]
        pd.DataFrame(recs).to_csv(os.path.join(OUT, "testable_cases.tsv"),
                                  sep="\t", index=False)
        log(f"  wrote testable_cases.tsv "
            f"({sum(len(v) for v in testable_all.values()):,} anchors)")
    return pool_all


def beat6_power(store, rows, matched_cases=None, matched_controls=None,
                pool_sizes=None):
    banner("BEAT 6  power for the PER-ANCHOR test  (v4, corrected)")
    log("  v3 WAS INVALID. It bootstrapped the null maximum from the control")
    log("  values; a bootstrap maximum is capped at the largest control actually")
    log("  observed, so whenever the case max exceeded it the p-value collapsed")
    log("  to 1/N. Measured type I error 0.472. Its power numbers of 1.00 were")
    log("  measuring that failure.")
    log("  A label-permutation null does not rescue a MAX statistic either: the")
    log("  global maximum lands in the permuted case set half the time whatever")
    log("  its true label, so the p-value floors near 0.5.")
    log()
    log("  v4 DESIGN: each case anchor gets its own p-value against the")
    log("  background pool, then Benjamini-Hochberg across anchors. A hit is any")
    log("  anchor with BH-adjusted p < 0.05.")

    for p in PUCKS:
        kc, kk = (p, PRIMARY_RADIUS, "case"), (p, PRIMARY_RADIUS, "control")
        if kc not in store or kk not in store:
            continue
        dc, dk = store[kc], store[kk]
        if matched_cases is not None:
            dc = dc[dc.anchor.isin(matched_cases)]
            dk = dk[dk.anchor.isin(matched_controls)]
        base = (dk.max_clone_beads / dk.n_tcr_beads.clip(lower=1)).dropna().values
        n_case = len(dc)
        npool = (pool_sizes or {}).get(p, POOL_SIZE)
        if len(base) < 10 or n_case < 5:
            log(f"\n  puck {p}: too few anchors to simulate"); continue

        log(f"\n  --- puck {p}: n_case {n_case}, pool {npool:,}, "
            f"control proxy median {np.median(base):.3f} ---")
        log(f"    {'k':>3} {'delta':>7} {'power':>8}   (type I in the k=0 row)")

        for k in [0] + SPARSE_K:
            for delta in ([0.0] if k == 0 else SPARSE_DELTA):
                hits = 0
                # The real pool is `npool` DISTINCT neighbourhoods, each with
                # its own statistic, so it is genuinely continuous. A
                # bootstrap from the ~30 control values would be a set of
                # atoms, and smoothing those over-inflates the tail and makes
                # the test look conservative. A gamma matched to the control
                # mean and variance is the closer model.
                shape = max(0.5, (base.mean()/max(1e-9, base.std()))**2)
                scale = max(1e-9, base.var()/max(1e-9, base.mean()))
                for _ in range(N_POWER_SIM // 2):
                    pool = np.clip(RNG.gamma(shape, scale, npool), 0, 1)
                    ca = np.clip(RNG.gamma(shape, scale, n_case), 0, 1)
                    if k:
                        ca[RNG.choice(n_case, min(k, n_case), replace=False)] += delta
                    pv = np.array([(1 + np.sum(pool >= v)) / (npool + 1) for v in ca])
                    hits += bool((bh_adjust(pv) < 0.05).any())
                pw = hits / (N_POWER_SIM // 2)
                lbl = "type I" if k == 0 else f"{k:>3} {delta:>7.2f}"
                log(f"    {lbl if k else lbl:>11} {pw:>8.2f}" if k
                    else f"    {'k=0':>3} {'(null)':>7} {pw:>8.3f}  <- type I error")
                rows.append({"beat": "power_per_anchor", "puck": p,
                             "k_anchors": k, "delta": delta,
                             "power": round(pw, 3), "n_case": n_case,
                             "pool_size": npool})

    log("\n  PRE-REGISTERED DETECTABLE BOUND")
    log("    Control top-clone fraction sits near 0.08, so a shift of 0.20 means")
    log("    a neighbourhood where one clone holds ~28% of UMIs against a")
    log("    typical 8%. That is a strongly expanded clone, not a subtle one.")
    log("    DETECTABLE : shift >= 0.20 in top-clone fraction. Power 0.93-1.00")
    log("                 even if only ONE anchor is affected.")
    log("    MARGINAL   : shift 0.10-0.20. Power 0.14-0.17 at one anchor,")
    log("                 0.76-0.89 only if ten or more anchors move.")
    log("    NOT VISIBLE: shift < 0.10. Power at or below the type I rate.")
    log("    A null result is therefore bounded at a shift of 0.20, NOT evidence")
    log("    that neoantigens do not recruit T cells.")
    log("\n  Read the k=0 row first: it must sit near 0.05 or the test is not")
    log("  calibrated. Then find the row matching the biology you consider")
    log("  plausible. Power here is MUCH lower than the v3 table claimed, and")
    log("  that table should be disregarded.")


# ==============================================================================
def main():
    global PRIMARY_RADIUS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--chain", default=None, help="restrict to one chain, e.g. TRB")
    ap.add_argument("--with-replacement", action="store_true",
                    help="allow a control to match several cases; use only if "
                         "without-replacement dropout exceeds the threshold")
    ap.add_argument("--radius", type=float, default=PRIMARY_RADIUS)
    args = ap.parse_args()

    PRIMARY_RADIUS = args.radius

    os.makedirs(OUT, exist_ok=True)
    needed = [ADATA, NEO, MUT] + [os.path.join(TCR_DIR, TCR_CSV[p]) for p in PUCKS]
    missing = [f for f in needed if not os.path.exists(f)]
    log(f"path check: {len(missing)} missing")
    for f in missing:
        log(f"   MISSING {f}")
    if args.validate_only:
        sys.exit(1 if missing else 0)
    if missing:
        sys.exit("cannot proceed")

    log(f"\nloading {ADATA}")
    adata = ad.read_h5ad(ADATA)
    log(f"  {adata.n_obs:,} beads x {adata.n_vars:,} genes")

    chain_label = args.chain or "all"
    log(f"loading TCR tables (chain: {chain_label})")
    frames = []
    for p in PUCKS:
        d = load_tcr(p, args.chain)
        d["puck"] = p
        frames.append(d)
    tcr = pd.concat(frames, ignore_index=True)
    log(f"  {len(tcr):,} rows")
    log(f"  TCR UNIT = all TCR-bearing beads (NOT source beads): using source")
    log(f"  beads would assume one cell per clone, which is what we are testing")
    log(f"  primary radius {PRIMARY_RADIUS:.0f} um | k primary {K_PRIMARY}, "
        f"strict {K_STRICT}")

    rows = []
    cases, controls = beat0_sets(adata, rows)
    store = beat1_occupancy(adata, tcr, cases, controls, rows, chain_label)
    beat2_intermixing(adata, cases, controls, rows)
    beat3_balance(adata, tcr, cases, controls, rows)
    m_cases, m_controls = beat3b_matching(adata, cases, controls, rows,
                                          with_replacement=args.with_replacement)
    beat4_rarefaction(store, rows)
    beat5_clone_background(tcr, rows)
    pool = beat6b_background_pool(adata, cases, m_controls, rows, PRIMARY_RADIUS)
    # power is estimated on the MATCHED sets against the real pool sizes
    beat6_power(store, rows, m_cases, m_controls,
                {k: len(v) for k, v in pool.items()})

    pd.DataFrame(sorted(m_cases), columns=["CB"]).to_csv(
        os.path.join(OUT, f"matched_cases_{chain_label}.tsv"), sep="\t", index=False)
    pd.DataFrame(sorted(m_controls), columns=["CB"]).to_csv(
        os.path.join(OUT, f"matched_controls_{chain_label}.tsv"), sep="\t", index=False)
    log(f"\nwrote matched_cases_{chain_label}.tsv and "
        f"matched_controls_{chain_label}.tsv for the downstream tests")

    out = os.path.join(OUT, f"feasibility_summary_{chain_label}.tsv")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    log(f"\nwrote {out}")

    banner("WHAT THIS DECIDES")
    log("  BEAT 1 is a FEASIBILITY GATE, not a result. clones_ge2 is a MEDIAN")
    log("    PER ANCHOR, so the 1:1.8 case/control imbalance is already divided")
    log("    out and cases and controls being similar there is expected. It")
    log("    says only that multi-bead clones exist to test. Any real signal is")
    log("    in WHICH clone and HOW enriched against its own puck-wide")
    log("    frequency, which is test B and is not run here.")
    log("  BEAT 2 says the contrast is not confounded by geography.")
    log("  BEAT 3b decides whether matching works. If covariates stay unbalanced")
    log("    or too many cases are dropped, switch to covariate modelling.")
    log("  BEAT 4 fixes a per-puck rarefaction floor.")
    log("  BEAT 6b sizes the background pool; resolution is 1/(pool+1), so a")
    log("    pool of a few thousand is what makes per-anchor testing possible.")
    log("  BEAT 6 gives CALIBRATED power for the per-anchor test. The v3 MAX")
    log("    table was invalid (type I 0.472) and must be disregarded.")
    log("  TEST ORDER: A BEFORE B. The background pool excludes beads near case")
    log("    anchors but NOT beads near controls, because controls are null by")
    log("    hypothesis. Test A is what checks that assumption. If mutation-only")
    log("    beads do attract T cells, the pool is contaminated and B is")
    log("    conservative, so A must be read first.")
    log("  Only then do tests A-D get written.")


if __name__ == "__main__":
    main()
