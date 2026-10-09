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
N_POWER_SIM = 2000

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
    log("\n  Imbalance here does not invalidate the design, but any imbalanced")
    log("  covariate must enter the downstream model or the matching.")


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
        for fl in (10, 20, 30, 50, 75, 100):
            n = int((u >= fl).sum())
            log(f"    {fl:>7} {n:>18,} {100*n/len(u):>6.1f}%")
            rows.append({"beat": "rarefaction", "puck": p, "floor": fl,
                         "retained": n, "pct_retained": round(100*n/len(u), 2)})
    log("\n  Pick the largest floor that retains most anchors; report how many")
    log("  fall below it rather than silently dropping them.")


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
def beat6_power(store, rows):
    banner("BEAT 6  power for the top-clone-fraction contrast")
    log("  Using the OBSERVED control distribution as the null, what shift in")
    log("  top-clone fraction would we detect at 80% power with these n?")
    log("  Simulation only; no test is run on the real cases.")
    for p in PUCKS:
        kc = (p, PRIMARY_RADIUS, "case")
        kk = (p, PRIMARY_RADIUS, "control")
        if kc not in store or kk not in store:
            continue
        dc, dk = store[kc], store[kk]
        # proxy for top-clone fraction from occupancy: max clone beads / tcr beads
        base = (dk.max_clone_beads / dk.n_tcr_beads.clip(lower=1)).dropna().values
        n_case, n_ctrl = len(dc), len(dk)
        if len(base) < 10 or n_case < 5:
            log(f"\n  puck {p}: too few anchors to simulate"); continue
        log(f"\n  puck {p}: n_case {n_case}, n_control {n_ctrl}, "
            f"control top-clone proxy median {np.median(base):.3f}")
        log(f"    {'shift':>8} {'power':>8}")
        detect = None
        for shift in (0.02, 0.05, 0.10, 0.15, 0.20, 0.30):
            hits = 0
            for _ in range(N_POWER_SIM // 4):
                a = RNG.choice(base, n_case, replace=True) + shift
                b = RNG.choice(base, n_ctrl, replace=True)
                try:
                    _, pv = mannwhitneyu(a, b, alternative="greater")
                except ValueError:
                    continue
                hits += (pv < 0.05)
            pw = hits / (N_POWER_SIM // 4)
            log(f"    {shift:>8.2f} {pw:>8.2f}")
            if detect is None and pw >= 0.80:
                detect = shift
            rows.append({"beat": "power", "puck": p, "shift": shift,
                         "power": round(pw, 3), "n_case": n_case, "n_control": n_ctrl})
        log(f"    smallest detectable shift at 80% power: "
            + (f"{detect:.2f}" if detect else "> 0.30, underpowered"))
    log("\n  This uses a proxy (max clone beads / TCR beads) rather than the real")
    log("  rarefied top-clone fraction, so treat it as an order-of-magnitude")
    log("  guide to whether the design is worth running at all.")


# ==============================================================================
def main():
    global PRIMARY_RADIUS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--chain", default=None, help="restrict to one chain, e.g. TRB")
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
    beat4_rarefaction(store, rows)
    beat5_clone_background(tcr, rows)
    beat6_power(store, rows)

    out = os.path.join(OUT, f"feasibility_summary_{chain_label}.tsv")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    log(f"\nwrote {out}")

    banner("WHAT THIS DECIDES")
    log("  BEAT 1 clones_ge2 is the gate. If the median is 0 at every radius,")
    log("    local expansion cannot be detected and the design stops here.")
    log("  BEAT 2 tells us whether the case/control contrast is confounded by")
    log("    geography. All three p-values above ~0.05 means it is not.")
    log("  BEAT 3 names any covariate that must enter the model.")
    log("  BEAT 4 fixes the rarefaction floor.")
    log("  BEAT 6 says whether the effect size we could detect is plausible.")
    log("  Only then do tests A-D get written.")


if __name__ == "__main__":
    main()
