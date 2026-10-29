#!/usr/bin/env python3
"""
Step08e_CD8_Phase0_Diagnostic.py
================================================================================
READ-ONLY (one exception: it backs up the Step08c design artifacts before
anything else). Runs NO hypothesis test. Its job is to decide whether a
CD8-restricted version of the Step08d tests can run at all, and on what terms.

WHY THIS EXISTS
    Step08d tests A-E use no cell-type information. The statistic is the
    rarefied top-clone fraction over TRB UMIs in a neighbourhood. TRB is
    expressed by CD8, CD4 and regulatory T cells alike, so CD4-lineage UMIs
    enter every neighbourhood and dilute the clonality signal. MHC class I
    presentation is read by CD8, so a CD8-focused statistic is the biologically
    matched one.

    BUT the dilution is SYMMETRIC. Test A pairs each anchor with its own
    depth-matched mutation-only control; tests B and D score each anchor
    against a background pool on the same puck. CD4 UMIs land in both arms.
    The expected effect of removing them is attenuation relief, not a signal
    appearing from nothing. Read any improvement in that light.

THE CONSTRAINT THAT DRIVES EVERY DESIGN HERE
    TCR-bearing beads: 18,927 / 28,861 / 23,816.
    T_cell-annotated beads, all three pucks combined: 5,184.
    Most TCR UMIs sit on beads that are not called T cells, which is what
    Step08b's diffusion kernel and the 71.7-75.0% error-derived barcode rate
    predict. A HARD subset of UMIs to CD8-called beads therefore discards well
    over 90% of the pool, and the Step08c rarefaction floors (20 / 75 / 50 UMIs
    at R = 100 um) would almost never be met. That is the v2 floor bug in a new
    costume: floors calibrated on one UMI pool, applied to a much smaller one.

    So the candidate designs are CONTINUOUS WEIGHTS, not subsets, and this
    script measures the weighted effective sample size against the real floors
    before any test is written.

THE THREE DESIGNS THIS SCRIPT SIZES
    1  CD8-weighted clonality. Keep every TRB UMI, weight it by its bead's CD8
       share of the T lineage. Tests A, B, D become weighted versions.
    2  CD8 abundance as the OUTCOME. Never run before. Local c2l CD8 around
       each anchor against its matched control. No rarefaction floor, no TCR
       sparsity. CD8 recruitment without clonal expansion is a real and
       reportable result, so CD4 and Treg are carried as separate outcomes
       rather than folded into a denominator.
    3  Subtype-stratified anchors. Split anchors by local CD8 dominance and run
       the existing tests within stratum. Weakest, because it halves an already
       small anchor set; sized here so that is a number rather than a guess.

BEATS
    0  back up the Step08c design artifacts (they are NOT chain-suffixed and
       every 08c run overwrites them)
    1  c2l T-lineage column inventory and total T-lineage mass per bead
    2  collinearity of c2l CD8 / CD4 / Treg. If these track each other the
       split carries no information and designs 1 and 3 are dead on arrival
    3  marker-gene cross-check of the c2l split (CD8A, CD8B, FOXP3, GZMB, ...)
    4  the weight itself: distribution, and whether it is a depth readout
    5  weighted UMI accounting per neighbourhood vs the REAL Step08c floors,
       swept across radius
    6  projected pair attrition under each weight definition (design 1 go/no-go)
    7  design 2 preview: local CD8 / CD4 / Treg, cases vs matched controls
    8  design 3 feasibility: anchors per stratum vs the Step08c minimum
    9  verdict

WHAT THIS SCRIPT IS NOT
    It does not fix the HLA problem. The anchors are still the 110 beads whose
    neoantigens were predicted against a fixed ten-allele panel rather than
    each patient's genotype. A CD8-restricted test on a partly non-presentable
    anchor set carries the same interpretive ceiling. This is a sensitivity
    analysis to have in hand before the remap, not a rescue of the 08d result.

ON BEAT 7 AND DOUBLE DIPPING
    Beat 7 is effectively design 2's test statistic. Looking at it here and
    then "testing" it in Step08f is not independent. So the primary radius is
    PRE-SPECIFIED HERE as 100 um, inherited from Step08d v4, before the numbers
    are seen. Every other radius in the sweep is reported as exploratory and
    must be labelled that way downstream.

USAGE
    python Step08e_CD8_Phase0_Diagnostic.py
    python Step08e_CD8_Phase0_Diagnostic.py --chain TRB
    python Step08e_CD8_Phase0_Diagnostic.py --radii 30 50 75 100 150 200
    python Step08e_CD8_Phase0_Diagnostic.py --skip-markers
    python Step08e_CD8_Phase0_Diagnostic.py --validate-only

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
================================================================================
"""

import argparse
import os
import shutil
import sys
from datetime import datetime

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, mannwhitneyu, wilcoxon, binomtest
from scipy import sparse as sp

try:
    import anndata as ad
except ImportError:
    sys.exit("anndata not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTPUTS = os.path.join(PROOT, "data/outputs")
FEAS = os.path.join(OUTPUTS, "11_neoantigen_tcr")
TESTS08D = os.path.join(OUTPUTS, "12_neoantigen_tcr_tests")
OUT = os.path.join(OUTPUTS, "13_cd8_phase0")
FIGDIR = os.path.join(OUT, "figures")

ADATA = os.path.join(OUTPUTS, "04_annotation/all_pucks_annotated_unified.h5ad")
NEO = os.path.join(OUTPUTS, "07_neoantigen/neoantigens_per_bead.tsv")
PAIRS = os.path.join(FEAS, "matched_pairs.tsv")
TESTABLE = os.path.join(FEAS, "testable_cases.tsv")
POOL = os.path.join(FEAS, "background_pool.tsv")
TCR_DIR = os.path.join(PROOT, "data/inputs/tcr/processed")

PUCKS = ["29", "37", "40"]
TCR_CSV = {p: f"B59_{p}_hTCR_tcr.csv" for p in PUCKS}

ANNOT_COL = "unified_annotation"
BASAL_T = "T_cell"
EPITHELIAL = "epithelial"
MARK = "marker_score_T_cell"

# The three c2l T-lineage columns. Note the comma inside the CD4 and CD8 names.
C_CD8 = "c2l_CD8-positive, alpha-beta T cell"
C_CD4 = "c2l_CD4-positive, alpha-beta T cell"
C_TREG = "c2l_regulatory T cell"

# Marker genes for the beat 3 cross-check of the c2l split.
#   CD4 THE GENE IS NOT A CD4 T-CELL MARKER on these beads: monocytes and
#   macrophages express it, and myeloid is the second largest compartment
#   (25,062 beads). It is reported for completeness and must not be used to
#   define a CD4 fraction.
MARKER_GENES = {
    "CD8A": "CD8 lineage",
    "CD8B": "CD8 lineage",
    "GZMB": "cytotoxic",
    "PRF1": "cytotoxic",
    "GZMK": "cytotoxic",
    "FOXP3": "Treg",
    "IL2RA": "Treg / activated",
    "CD4": "NOT CD8-exclusive, myeloid expressed",
    "IL7R": "CD4 memory / also non-T",
    "CD3E": "pan-T",
    "TRBC2": "pan-T, TRB constant",
    "PTPRC": "pan-immune",
}

# v4 of Step08d dropped 75 um and settled on 100. Jake asked for a sweep here
# because the weighted UMI count per neighbourhood is the binding constraint
# and it scales with area.
RADII = [30.0, 50.0, 75.0, 100.0, 150.0, 200.0]
PRIMARY_RADIUS = 100.0        # PRE-SPECIFIED for design 2, see header
DEFAULT_CHAIN = "TRB"

# Beads whose total T-lineage abundance is near zero have an unstable CD8
# share: a ratio of two tiny numbers. The floor is a percentile of the
# T-lineage mass observed on T_cell-ANNOTATED beads in that puck, so it is
# calibrated to real T-cell beads rather than to an arbitrary constant.
MASS_FLOOR_PCT = 25.0
MASS_FLOOR_GRID = [0.0, 25.0, 50.0]   # sensitivity, reported in beat 5

# Step08c's floor for a testable anchor. Reused here so design 3's strata are
# judged against the same bar the original design used.
POOL_MIN_CASES = 8

EPS = 1e-12

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
C_CORAL = "#ed6a5a"     # cases, neoantigen
C_MUSTARD = "#F6D155"   # controls, mutation only
C_GRAY = "#d3d3d3"
C_DARK = "#4a4a4a"
C_BLUE = "#7eb0d5"      # CD8
C_GREEN = "#b2e061"     # CD4
C_SALMON = "#fd7f6f"    # Treg
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

RNG = np.random.default_rng(20261027)


# ==============================================================================
# helpers
# ==============================================================================
def log(m=""):
    print(m, flush=True)


def banner(m):
    log()
    log("=" * 78)
    log(m)
    log("=" * 78)


def savefig(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), dpi=DPI,
                    bbox_inches="tight")
    plt.close(fig)
    log(f"  figure: {name}.pdf / .png")


def style(ax, title=None, xlabel=None, ylabel=None):
    if title:
        ax.set_title(title, fontsize=FS_TITLE)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=FS_LABEL)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=FS_LABEL)
    ax.tick_params(labelsize=FS_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def smd(a, b):
    """Standardized mean difference, the Step08c balance diagnostic."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    spool = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2.0)
    return float((a.mean() - b.mean()) / spool) if spool > 0 else np.nan


def kish_ess(w):
    """Kish effective sample size, (sum w)^2 / sum w^2.

    This is how many effective draws a weighted rarefaction gets out of the
    neighbourhood, so it is the quantity the Step08c floor applies to: a
    weighted draw of `floor` UMIs with probability proportional to w needs
    ess >= floor.

    IT IS SCALE INVARIANT, and that matters. 200 UMIs all at w = 0.10 give
    ess = 200, not 20, because uniform weights are no weighting at all. That
    is arithmetically right and substantively a warning: a neighbourhood with
    a large ess and a tiny sum(w) is one where the CD8 weight is flat, so the
    weighted test is sampling the same CD4-contaminated pool Step08d used and
    has simply relabelled it. BEAT 5 therefore reports ess AND sum(w), and the
    two answer different questions:
        ess     can a weighted test RUN here (floor criterion)
        sum(w)  is there any CD8 contrast here to run it ON
    """
    w = np.asarray(w, float)
    w = w[np.isfinite(w) & (w > 0)]
    if not len(w):
        return 0.0
    s = w.sum()
    return float(s * s / np.square(w).sum())


# ==============================================================================
# loading
# ==============================================================================
def validate_paths(require_markers=True):
    ok = True
    need = [("adata", ADATA), ("neoantigens", NEO), ("matched_pairs", PAIRS),
            ("testable_cases", TESTABLE), ("background_pool", POOL)]
    for p in PUCKS:
        need.append((f"tcr {p}", os.path.join(TCR_DIR, TCR_CSV[p])))
    for name, path in need:
        hit = os.path.exists(path)
        log(f"  {'OK ' if hit else '***'} {name:<22} {path}")
        ok &= hit
    return ok


def read_floors(chain):
    """Floors come from the Step08c summary for the SAME CHAIN. Never hardcoded.

    v2 of Step08d hardcoded 20/75/50, which Step08c had chosen on chain=all
    (269,594 rows) while 08d ran chain=TRB (198,037 rows). The same
    neighbourhoods held ~26% fewer UMIs than the floors were sized for and pair
    attrition at R=75 hit 8 of 18 in puck 37. Hard stop if the file is absent.
    """
    path = os.path.join(FEAS, f"feasibility_summary_{chain}.tsv")
    if not os.path.exists(path):
        sys.exit(f"FATAL: {path} not found. Run Step08c --chain {chain} first. "
                 f"Floors are never hardcoded.")
    d = pd.read_csv(path, sep="\t")
    rows = d[d.beat == "rarefaction_chosen"]
    if not len(rows):
        sys.exit(f"FATAL: no 'rarefaction_chosen' rows in {path}.")
    floors = {str(r.puck): int(r.floor) for r in rows.itertuples()}
    log(f"  floors read from feasibility_summary_{chain}.tsv: {floors}")
    return floors


def load_tcr(chain):
    frames = []
    for p in PUCKS:
        d = pd.read_csv(os.path.join(TCR_DIR, TCR_CSV[p]))
        d["puck"] = p
        d["bead"] = d["bc"].astype(str) + f"-1_Puck_211214_{p}"
        frames.append(d)
    t = pd.concat(frames, ignore_index=True)
    log(f"  TCR rows, all chains: {len(t):,}")
    if chain and chain != "all":
        keep = [c.strip() for c in chain.split(",")]
        t = t[t["chain"].isin(keep)].copy()
        log(f"  after chain filter {keep}: {len(t):,} rows "
            f"({100*len(t)/max(1,len(frames[0])+len(frames[1])+len(frames[2])):.1f}% of all)")
    for p in PUCKS:
        d = t[t.puck == p]
        log(f"    puck {p}: {len(d):,} rows, {d.bead.nunique():,} beads, "
            f"{d.cloneId.nunique():,} clones")
    return t


def load_bead_frame(adata, p, mass_floor_pct):
    """Per-bead frame for one puck, with every candidate CD8 weight attached."""
    m = adata.obs["puck_id"].astype(str).values == f"Puck_211214_{p}"
    idx = adata.obs_names[m]
    o = adata.obs.loc[idx]

    f = pd.DataFrame(index=idx)
    f["depth"] = o["total_counts"].values.astype(float)
    f["annot"] = o[ANNOT_COL].astype(str).values
    f["is_T"] = (f["annot"] == BASAL_T).astype(int)
    f["marker_T"] = (o[MARK].values.astype(float)
                     if MARK in o.columns else np.nan)
    f["cd8"] = o[C_CD8].values.astype(float) if C_CD8 in o.columns else np.nan
    f["cd4"] = o[C_CD4].values.astype(float) if C_CD4 in o.columns else np.nan
    f["treg"] = o[C_TREG].values.astype(float) if C_TREG in o.columns else np.nan

    f["t_mass"] = f[["cd8", "cd4", "treg"]].sum(axis=1, min_count=1)

    # mass floor, calibrated on T_cell-annotated beads in THIS puck
    tm = f.loc[f.is_T == 1, "t_mass"].dropna().values
    floor = float(np.percentile(tm, mass_floor_pct)) if len(tm) else 0.0

    ok = f["t_mass"].values > max(floor, EPS)

    # --- candidate weights -----------------------------------------------
    # w_share   CD8 share of T lineage. Scale free, so a shallow bead with a
    #           little of everything is treated the same as a deep one. Set to
    #           0 below the mass floor rather than NaN: a bead with no T-lineage
    #           signal contributes nothing, which is the right behaviour, not
    #           a missing value that silently drops its UMIs.
    share = np.where(ok, f["cd8"].values / np.maximum(f["t_mass"].values, EPS), 0.0)
    f["w_share"] = share

    # w_abs     absolute CD8 abundance, scaled to the puck 99th percentile and
    #           clipped. Keeps "how much CD8 is here" rather than "what fraction
    #           of T lineage is CD8", at the cost of inheriting any depth
    #           dependence in c2l output. Beat 4 measures that dependence.
    q99 = np.nanpercentile(f["cd8"].values, 99) if f["cd8"].notna().any() else np.nan
    f["w_abs"] = (np.clip(f["cd8"].values / q99, 0, 1)
                  if np.isfinite(q99) and q99 > 0 else 0.0)

    # w_argmax  hard call: CD8 is the largest of the three AND mass clears the
    #           floor. Included so the cost of a hard call is visible as a
    #           number next to the continuous alternatives, not argued about.
    tri = f[["cd8", "cd4", "treg"]].values
    arg_cd8 = np.zeros(len(f), bool)
    good = np.isfinite(tri).all(axis=1)
    arg_cd8[good] = tri[good].argmax(axis=1) == 0
    f["w_argmax"] = (arg_cd8 & ok).astype(float)

    # w_flat    no weighting at all. The Step08d status quo, carried through
    #           every table so every CD8 number has its own control.
    f["w_flat"] = 1.0

    xy = o[["x_coord", "y_coord"]].values.astype(float)
    pos = {b: i for i, b in enumerate(idx)}
    return f, xy, pos, floor


WEIGHTS = [("w_flat", "unweighted (08d status quo)"),
           ("w_share", "CD8 share of T lineage"),
           ("w_abs", "CD8 abundance, puck-scaled"),
           ("w_argmax", "hard call, CD8 is argmax")]


# ==============================================================================
# BEAT 0
# ==============================================================================
def beat0_backup(do_backup, rows):
    banner("BEAT 0  back up the Step08c design artifacts")
    log("  matched_pairs.tsv, testable_cases.tsv and background_pool.tsv are")
    log("  NOT chain-suffixed and are overwritten by every Step08c run. Any")
    log("  recalibration for a weighted pool would destroy the sets the 08d")
    log("  results were computed on. Copy first, argue later.")
    if not do_backup:
        log("  --no-backup given, skipping")
        return
    stamp = datetime.now().strftime("%Y%m%d")
    dest = os.path.join(FEAS, f"BACKUP_pre_cd8_{stamp}")
    if os.path.isdir(dest) and os.listdir(dest):
        log(f"  backup already present, leaving it alone: {dest}")
        rows.append({"beat": "backup", "path": dest, "action": "already_present"})
        return
    os.makedirs(dest, exist_ok=True)
    n = 0
    for src in (PAIRS, TESTABLE, POOL):
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dest, os.path.basename(src)))
            log(f"    copied {os.path.basename(src)}")
            n += 1
    # the chain-suffixed summaries too, cheap and they pin the floors
    for fn in sorted(os.listdir(FEAS)):
        if fn.startswith("feasibility_summary_") and fn.endswith(".tsv"):
            shutil.copy2(os.path.join(FEAS, fn), os.path.join(dest, fn))
            log(f"    copied {fn}")
            n += 1
    # and the 08d result tables, so the comparison arm is frozen
    if os.path.isdir(TESTS08D):
        for fn in sorted(os.listdir(TESTS08D)):
            if fn.endswith(".tsv"):
                shutil.copy2(os.path.join(TESTS08D, fn),
                             os.path.join(dest, f"08d_{fn}"))
                log(f"    copied 08d_{fn}")
                n += 1
    log(f"  {n} files -> {dest}")
    rows.append({"beat": "backup", "path": dest, "n_files": n, "action": "created"})


# ==============================================================================
# BEAT 1
# ==============================================================================
def beat1_inventory(adata, mass_floor_pct, rows):
    banner("BEAT 1  c2l T-lineage column inventory")
    cols = adata.obs.columns
    for name, c in (("CD8", C_CD8), ("CD4", C_CD4), ("Treg", C_TREG)):
        log(f"  {name:<5} {'FOUND' if c in cols else '*** MISSING ***'}  {c}")
    missing = [c for c in (C_CD8, C_CD4, C_TREG) if c not in cols]
    if missing:
        log()
        log("  Available c2l columns, for repair:")
        for c in sorted([c for c in cols if c.startswith("c2l_")]):
            log(f"    {c}")
        sys.exit("FATAL: a required c2l T-lineage column is absent. "
                 "Designs 1 and 3 cannot be defined without it.")

    log()
    log("  Per puck: T-lineage abundance, all beads and T_cell-annotated beads.")
    for p in PUCKS:
        f, _, _, floor = load_bead_frame(adata, p, mass_floor_pct)
        nT = int(f.is_T.sum())
        log(f"\n  --- puck {p}  ({len(f):,} beads, {nT:,} annotated T_cell) ---")
        log(f"    {'measure':<10} {'median all':>12} {'median T':>12} "
            f"{'frac zero':>11} {'frac zero T':>12}")
        for lab, c in (("CD8", "cd8"), ("CD4", "cd4"), ("Treg", "treg"),
                       ("T mass", "t_mass")):
            v = f[c].values
            vt = f.loc[f.is_T == 1, c].values
            log(f"    {lab:<10} {np.nanmedian(v):>12.4f} {np.nanmedian(vt):>12.4f} "
                f"{np.nanmean(v <= EPS):>11.3f} {np.nanmean(vt <= EPS):>12.3f}")
            rows.append({"beat": "c2l_inventory", "puck": p, "measure": lab,
                         "median_all": float(np.nanmedian(v)),
                         "median_Tcell": float(np.nanmedian(vt)),
                         "frac_zero_all": float(np.nanmean(v <= EPS)),
                         "frac_zero_Tcell": float(np.nanmean(vt <= EPS))})
        log(f"    T-lineage mass floor at p{mass_floor_pct:.0f} of T beads: {floor:.4f}")
        log(f"    beads clearing that floor: "
            f"{100*np.mean(f.t_mass.values > max(floor, EPS)):.1f}% of all, "
            f"{100*np.mean(f.loc[f.is_T==1,'t_mass'].values > max(floor, EPS)):.1f}% of T beads")
        rows.append({"beat": "mass_floor", "puck": p,
                     "floor_pct": mass_floor_pct, "floor_value": floor,
                     "frac_beads_above": float(np.mean(f.t_mass.values > max(floor, EPS)))})

    log()
    log("  A high 'frac zero' on T_cell beads is the warning sign. cell2location")
    log("  ran with N_cells_per_location = 1 on beads where T cells are 24-26%")
    log("  shallower than average; if CD8 is zero on most T beads the split")
    log("  cannot carry a weight and design 2 is the only survivor.")


# ==============================================================================
# BEAT 2
# ==============================================================================
def beat2_collinearity(adata, mass_floor_pct, tcr, rows):
    banner("BEAT 2  collinearity of c2l CD8 / CD4 / Treg")
    log("  THE GO/NO-GO FOR DESIGNS 1 AND 3. If cell2location cannot separate")
    log("  the three T subtypes on these beads, the three abundances move")
    log("  together, the CD8 share is a constant plus noise, and weighting by")
    log("  it reproduces the 08d result with extra machinery. Spearman is")
    log("  reported on three bead sets because the answer can differ: the")
    log("  deconvolution may separate them where there is real T signal and")
    log("  not elsewhere.")

    fig, axes = plt.subplots(1, len(PUCKS), figsize=(13 * len(PUCKS), 12))
    for i, p in enumerate(PUCKS):
        f, _, _, floor = load_bead_frame(adata, p, mass_floor_pct)
        tb = set(tcr.loc[tcr.puck == p, "bead"])
        sets = (("all beads", np.ones(len(f), bool)),
                ("T_cell beads", f.is_T.values == 1),
                ("TCR-bearing", np.array([b in tb for b in f.index])),
                ("above mass floor", f.t_mass.values > max(floor, EPS)))
        log(f"\n  --- puck {p} ---")
        log(f"    {'bead set':<20} {'n':>8} {'rho CD8:CD4':>13} "
            f"{'rho CD8:Treg':>14} {'rho CD4:Treg':>14}")
        for name, m in sets:
            sub = f.loc[m]
            if len(sub) < 50:
                log(f"    {name:<20} {len(sub):>8,}  too few beads")
                continue
            r = {}
            for a, b, k in (("cd8", "cd4", "cd8_cd4"), ("cd8", "treg", "cd8_treg"),
                            ("cd4", "treg", "cd4_treg")):
                v1, v2 = sub[a].values, sub[b].values
                ok = np.isfinite(v1) & np.isfinite(v2)
                r[k] = float(spearmanr(v1[ok], v2[ok]).statistic) if ok.sum() > 10 else np.nan
            log(f"    {name:<20} {len(sub):>8,} {r['cd8_cd4']:>13.3f} "
                f"{r['cd8_treg']:>14.3f} {r['cd4_treg']:>14.3f}")
            rows.append({"beat": "collinearity", "puck": p, "bead_set": name,
                         "n": int(len(sub)), **{f"rho_{k}": v for k, v in r.items()}})

        ax = axes[i] if len(PUCKS) > 1 else axes
        m = f.is_T.values == 1
        ax.scatter(f.loc[m, "cd4"].values, f.loc[m, "cd8"].values, s=40,
                   c=C_BLUE, alpha=0.45, linewidths=0, label="T_cell beads")
        lim = np.nanpercentile(np.r_[f.loc[m, "cd4"].values,
                                     f.loc[m, "cd8"].values], 99)
        if np.isfinite(lim) and lim > 0:
            ax.plot([0, lim], [0, lim], color=C_DARK, lw=2, ls="--")
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
        style(ax, title=f"Puck {p}",
              xlabel="c2l CD4 abundance",
              ylabel="c2l CD8 abundance" if i == 0 else None)
        if i == 0:
            ax.legend(fontsize=FS_ANNOT, frameon=False)
    fig.tight_layout()
    savefig(fig, "Fig_08e_A_cd8_vs_cd4_Tbeads")

    log()
    log("  Reading guide. |rho| above ~0.8 on T_cell beads means the two")
    log("  abundances are one measurement wearing two labels: stop, and run")
    log("  design 2 only. Between 0.4 and 0.8 the weight carries some signal")
    log("  and the continuous designs are worth running with the caveat")
    log("  reported. Below 0.4 the split is informative.")


# ==============================================================================
# BEAT 3
# ==============================================================================
def load_marker_matrix(genes):
    """Pull a few gene columns out of the annotated object without loading X.

    X in all_pucks_annotated_unified.h5ad is normalized, not raw counts, so
    only DETECTION (value > 0) is used below. Detection is invariant to the
    normalization; magnitudes are not, and are deliberately not interpreted.
    """
    try:
        a = ad.read_h5ad(ADATA, backed="r")
        present = [g for g in genes if g in a.var_names]
        if not present:
            return None, []
        sub = a[:, present].to_memory()
    except Exception as e:
        log(f"  backed gene slice failed ({e}); falling back to a full read")
        a = ad.read_h5ad(ADATA)
        present = [g for g in genes if g in a.var_names]
        if not present:
            return None, []
        sub = a[:, present]
    X = sub.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    return pd.DataFrame(X, index=sub.obs_names, columns=present), present


def beat3_markers(adata, mass_floor_pct, tcr, rows, skip):
    banner("BEAT 3  marker-gene cross-check of the c2l split")
    if skip:
        log("  --skip-markers given, skipping")
        return
    log("  cell2location is a model. Before a weight is built on its CD8/CD4")
    log("  split, that split should agree with the genes. The test: is the c2l")
    log("  CD8 share higher on beads where CD8A or CD8B is detected than on")
    log("  beads where it is not, among T_cell beads at matched depth?")
    log()
    log("  CD4 THE GENE is reported but is NOT usable as a CD4 T-cell marker")
    log("  here: monocytes and macrophages express it and myeloid is 25,062")
    log("  beads, the second largest compartment.")

    M, present = load_marker_matrix(list(MARKER_GENES))
    if M is None:
        log("  none of the marker genes are in var_names, skipping")
        return
    absent = [g for g in MARKER_GENES if g not in present]
    if absent:
        log(f"  not in var_names, skipped: {', '.join(absent)}")

    for p in PUCKS:
        f, _, _, _ = load_bead_frame(adata, p, mass_floor_pct)
        tb = set(tcr.loc[tcr.puck == p, "bead"])
        sub = M.reindex(f.index)
        is_t = f.is_T.values == 1
        is_tcr = np.array([b in tb for b in f.index])
        log(f"\n  --- puck {p}  detection rate (fraction of beads with > 0) ---")
        log(f"    {'gene':<8} {'role':<34} {'all':>8} {'T_cell':>9} {'TCR+':>8}")
        for g in present:
            v = sub[g].values > 0
            log(f"    {g:<8} {MARKER_GENES[g][:33]:<34} {v.mean():>8.3f} "
                f"{v[is_t].mean():>9.3f} {v[is_tcr].mean():>8.3f}")
            rows.append({"beat": "marker_detection", "puck": p, "gene": g,
                         "frac_all": float(v.mean()),
                         "frac_Tcell": float(v[is_t].mean()),
                         "frac_TCRpos": float(v[is_tcr].mean())})

        # the agreement test
        cd8_genes = [g for g in ("CD8A", "CD8B") if g in present]
        if cd8_genes and is_t.sum() > 50:
            pos = (sub[cd8_genes].values > 0).any(axis=1) & is_t
            neg = (~(sub[cd8_genes].values > 0).any(axis=1)) & is_t
            if pos.sum() >= 10 and neg.sum() >= 10:
                a, b = f.loc[pos, "w_share"].values, f.loc[neg, "w_share"].values
                try:
                    _, pv = mannwhitneyu(a, b, alternative="greater")
                except ValueError:
                    pv = np.nan
                log(f"\n    c2l CD8 share on T beads: "
                    f"CD8A/B detected {np.median(a):.3f} (n={pos.sum():,}) vs "
                    f"not detected {np.median(b):.3f} (n={neg.sum():,})")
                log(f"    one-sided MWU p = {pv:.3g} | SMD {smd(a, b):+.3f}")
                log(f"    median depth {np.median(f.loc[pos,'depth']):,.0f} vs "
                    f"{np.median(f.loc[neg,'depth']):,.0f}  <- if these differ")
                log(f"    materially the comparison is partly a depth readout")
                rows.append({"beat": "marker_agreement", "puck": p,
                             "n_cd8_detected": int(pos.sum()),
                             "n_cd8_not": int(neg.sum()),
                             "median_share_detected": float(np.median(a)),
                             "median_share_not": float(np.median(b)),
                             "smd": smd(a, b),
                             "p_mwu": float(pv) if pv == pv else None,
                             "median_depth_detected": float(np.median(f.loc[pos, "depth"])),
                             "median_depth_not": float(np.median(f.loc[neg, "depth"]))})
            else:
                log("\n    too few beads on one side for the agreement test")

    log()
    log("  If the c2l share does not separate CD8A/B-positive from negative T")
    log("  beads, the deconvolution is not recovering the subtype and the")
    log("  weight is noise, whatever beat 2 said about collinearity.")


# ==============================================================================
# BEAT 4
# ==============================================================================
def beat4_weight(adata, mass_floor_pct, rows):
    banner("BEAT 4  the weight itself: shape, and is it a depth readout?")
    log("  A weight that tracks sequencing depth would push the clonality")
    log("  statistic toward deep beads, and deep beads carry more UMIs, which")
    log("  is the exact confound Step08a BEAT 2 was rewritten to remove. If")
    log("  w_share rises monotonically across depth deciles, weighted tests")
    log("  need depth in the model rather than assumed away.")

    fig, axes = plt.subplots(1, len(PUCKS), figsize=(13 * len(PUCKS), 11))
    for i, p in enumerate(PUCKS):
        f, _, _, _ = load_bead_frame(adata, p, mass_floor_pct)
        log(f"\n  --- puck {p} ---")
        for col, lab in WEIGHTS:
            v = f[col].values
            log(f"    {lab:<34} mean {np.nanmean(v):.4f} | "
                f"median {np.nanmedian(v):.4f} | "
                f"frac zero {np.nanmean(v <= EPS):.3f} | "
                f"frac > 0.5 {np.nanmean(v > 0.5):.3f}")
            rows.append({"beat": "weight_shape", "puck": p, "weight": col,
                         "mean": float(np.nanmean(v)),
                         "median": float(np.nanmedian(v)),
                         "frac_zero": float(np.nanmean(v <= EPS)),
                         "frac_gt_half": float(np.nanmean(v > 0.5))})

        f["dep_dec"] = pd.qcut(f["depth"], 10, labels=False, duplicates="drop")
        log(f"\n    w_share by depth decile (T_cell beads only)")
        log(f"      {'decile':>7} {'median depth':>13} {'n':>7} {'mean w_share':>14}")
        tt = f[f.is_T == 1]
        dd_mean = []
        for dd, g in tt.groupby("dep_dec"):
            if len(g) < 10:
                continue
            mw = float(np.nanmean(g["w_share"].values))
            dd_mean.append((int(dd), mw))
            log(f"      {int(dd):>7} {g.depth.median():>13,.0f} {len(g):>7,} {mw:>14.4f}")
        if len(dd_mean) >= 3:
            d_, w_ = np.array(dd_mean).T
            rho = float(spearmanr(d_, w_).statistic)
            log(f"      Spearman(decile, mean w_share) = {rho:+.3f}")
            log("      |rho| near 1 means the weight is partly a depth measure.")
            rows.append({"beat": "weight_depth", "puck": p,
                         "spearman_decile_wshare": rho})

        ax = axes[i] if len(PUCKS) > 1 else axes
        ax.hist(f.loc[f.is_T == 1, "w_share"].values, bins=40, color=C_BLUE,
                alpha=0.85, label="T_cell beads")
        ax.hist(f.loc[f.is_T == 0, "w_share"].values, bins=40, color=C_GRAY,
                alpha=0.55, label="other beads", density=False)
        ax.set_yscale("log")
        style(ax, title=f"Puck {p}", xlabel="CD8 share of T lineage",
              ylabel="beads (log)" if i == 0 else None)
        if i == 0:
            ax.legend(fontsize=FS_ANNOT, frameon=False)
    fig.tight_layout()
    savefig(fig, "Fig_08e_B_weight_distribution")


# ==============================================================================
# BEAT 5
# ==============================================================================
def neighbourhood_rows(tree, xy, pos, idx, bead_rows, anchor, R):
    """Row indices of the TCR table that fall inside R um of `anchor`."""
    if anchor not in pos:
        return None
    near = tree.query_ball_point(xy[pos[anchor]], R)
    out = []
    for j in near:
        r = bead_rows.get(idx[j])
        if r is not None:
            out.append(r)
    return np.concatenate(out) if out else np.array([], int)


def build_bead_rows(tp, idx_set):
    """bead -> array of TCR row positions, restricted to annotated beads."""
    d = tp[tp.bead.isin(idx_set)]
    g = {}
    for bead, sub in d.groupby("bead", sort=False):
        g[bead] = sub.index.values
    return g, d


def beat5_accounting(adata, tcr, floors, radii, mass_floor_pct, cases_by_puck,
                     controls_by_puck, rows):
    banner("BEAT 5  weighted UMI accounting vs the real Step08c floors")
    log("  The binding constraint for design 1. For each anchor neighbourhood:")
    log("    n_raw      UMIs present, what Step08d counted")
    log("    sum(w)     CD8-equivalent UMIs, the weighted total")
    log("    ess        Kish effective sample size, (sum w)^2 / sum w^2")
    log()
    log("  THE FLOOR APPLIES TO ESS, not to sum(w). A weighted rarefaction")
    log("  drawing `floor` UMIs with probability proportional to w needs")
    log("  ess >= floor.")
    log()
    log("  ESS IS SCALE INVARIANT, so read it next to sum(w), never alone.")
    log("  200 UMIs all at w = 0.10 give ess = 200 and sum(w) = 20: the test")
    log("  can run, but the weights are flat, so it is sampling the same")
    log("  CD4-contaminated pool Step08d used under a new name. A high ESS with")
    log("  a low sum(w) is a BEAT 2 collinearity failure surfacing locally, not")
    log("  a pass. Both columns are reported for exactly that reason.")
    log(f"  Floors in use (from Step08c, same chain): {floors}")

    acc = []
    for p in PUCKS:
        f, xy, pos, _ = load_bead_frame(adata, p, mass_floor_pct)
        idx = f.index.values
        tree = cKDTree(xy)
        tp = tcr[tcr.puck == p].reset_index(drop=True)
        bead_rows, _ = build_bead_rows(tp, set(idx))
        wmap = {c: f[c].to_dict() for c, _ in WEIGHTS}
        row_w = {c: np.array([wmap[c].get(b, 0.0) for b in tp["bead"].values])
                 for c, _ in WEIGHTS}
        floor = floors.get(p)
        anchors = list(cases_by_puck.get(p, [])) + list(controls_by_puck.get(p, []))
        kind = (["case"] * len(cases_by_puck.get(p, []))
                + ["control"] * len(controls_by_puck.get(p, [])))
        log(f"\n  --- puck {p}  (floor {floor} UMIs, {len(anchors)} anchors) ---")
        log(f"    {'R um':>6} {'weight':<30} {'median n_raw':>13} "
            f"{'med sum(w)':>11} {'median ESS':>12} {'% >= floor':>11}")
        for R in radii:
            per_anchor = {c: [] for c, _ in WEIGHTS}
            nraw = []
            for a in anchors:
                ri = neighbourhood_rows(tree, xy, pos, idx, bead_rows, a, R)
                if ri is None:
                    nraw.append(np.nan)
                    for c, _ in WEIGHTS:
                        per_anchor[c].append(np.nan)
                    continue
                nraw.append(len(ri))
                for c, _ in WEIGHTS:
                    w = row_w[c][ri] if len(ri) else np.array([])
                    per_anchor[c].append((kish_ess(w), float(np.nansum(w))))
            nraw = np.array(nraw, float)
            for c, lab in WEIGHTS:
                arr = np.array(per_anchor[c], float)
                e, sw = arr[:, 0], arr[:, 1]
                frac = float(np.nanmean(e >= floor)) if floor else np.nan
                log(f"    {R:>6.0f} {lab:<30} {np.nanmedian(nraw):>13.0f} "
                    f"{np.nanmedian(sw):>11.1f} {np.nanmedian(e):>12.1f} "
                    f"{100*frac:>10.1f}%")
                rows.append({"beat": "umi_accounting", "puck": p, "radius": R,
                             "weight": c, "floor": floor,
                             "median_n_raw": float(np.nanmedian(nraw)),
                             "median_sum_w": float(np.nanmedian(sw)),
                             "median_ess": float(np.nanmedian(e)),
                             "frac_above_floor": frac})
                for a, k, e1, s1, n1 in zip(anchors, kind, e, sw, nraw):
                    acc.append({"puck": p, "radius": R, "weight": c,
                                "anchor": a, "kind": k, "n_raw": n1,
                                "sum_w": s1, "ess": e1,
                                "above_floor": bool(e1 >= floor) if floor else None})

        # mass-floor sensitivity, primary radius only, w_share only
        log(f"\n    mass-floor sensitivity at R = {PRIMARY_RADIUS:.0f} um, w_share")
        for mp in MASS_FLOOR_GRID:
            f2, xy2, pos2, fl2 = load_bead_frame(adata, p, mp)
            w2 = f2["w_share"].to_dict()
            rw = np.array([w2.get(b, 0.0) for b in tp["bead"].values])
            e = []
            for a in anchors:
                ri = neighbourhood_rows(tree, xy2, pos2, f2.index.values,
                                        bead_rows, a, PRIMARY_RADIUS)
                e.append(kish_ess(rw[ri]) if ri is not None and len(ri) else 0.0)  # ESS
            e = np.array(e, float)
            log(f"      p{mp:>5.0f}  mass floor {fl2:>8.4f}  median ESS "
                f"{np.nanmedian(e):>8.1f}  % >= floor "
                f"{100*np.nanmean(e >= floor):>6.1f}%")
            rows.append({"beat": "mass_floor_sensitivity", "puck": p,
                         "floor_pct": mp, "mass_floor": fl2,
                         "median_ess": float(np.nanmedian(e)),
                         "frac_above_floor": float(np.nanmean(e >= floor))})

    pd.DataFrame(acc).to_csv(os.path.join(OUT, "weight_accounting.tsv"),
                             sep="\t", index=False)
    log(f"\n  per-anchor detail: {os.path.join(OUT, 'weight_accounting.tsv')}")


# ==============================================================================
# BEAT 6
# ==============================================================================
def beat6_attrition(adata, tcr, floors, radii, mass_floor_pct, pairs, rows):
    banner("BEAT 6  projected pair attrition  (DESIGN 1 GO / NO-GO)")
    log("  Test A is paired. A pair survives only if BOTH members clear the")
    log("  floor. Step08d v4 dropped 75 um precisely because it lost 26 of 30")
    log("  and 16 of 18 pairs that way. The same question, asked of the")
    log("  weighted pool, before anything is written.")

    out = []
    for p in PUCKS:
        f, xy, pos, _ = load_bead_frame(adata, p, mass_floor_pct)
        idx = f.index.values
        tree = cKDTree(xy)
        tp = tcr[tcr.puck == p].reset_index(drop=True)
        bead_rows, _ = build_bead_rows(tp, set(idx))
        row_w = {c: np.array([f[c].to_dict().get(b, 0.0)
                              for b in tp["bead"].values])
                 for c, _ in WEIGHTS}
        floor = floors.get(p)
        pp = [(a, b) for a, b in pairs.get(p, [])]
        log(f"\n  --- puck {p}  ({len(pp)} matched pairs, floor {floor}) ---")
        log(f"    {'R um':>6} " + " ".join(f"{c.replace('w_',''):>12}" for c, _ in WEIGHTS))
        for R in radii:
            line = {}
            for c, _ in WEIGHTS:
                keep = 0
                for a, b in pp:
                    ra = neighbourhood_rows(tree, xy, pos, idx, bead_rows, a, R)
                    rb = neighbourhood_rows(tree, xy, pos, idx, bead_rows, b, R)
                    if ra is None or rb is None:
                        continue
                    ea = kish_ess(row_w[c][ra]) if len(ra) else 0.0
                    eb = kish_ess(row_w[c][rb]) if len(rb) else 0.0
                    if floor and ea >= floor and eb >= floor:
                        keep += 1
                line[c] = keep
                rows.append({"beat": "pair_attrition", "puck": p, "radius": R,
                             "weight": c, "n_pairs": len(pp), "n_retained": keep,
                             "frac_retained": keep / max(1, len(pp))})
                out.append({"puck": p, "radius": R, "weight": c,
                            "n_pairs": len(pp), "n_retained": keep})
            log(f"    {R:>6.0f} " + " ".join(f"{line[c]:>7}/{len(pp):<4}"
                                             for c, _ in WEIGHTS))

    df = pd.DataFrame(out)
    df.to_csv(os.path.join(OUT, "pair_attrition.tsv"), sep="\t", index=False)

    fig, axes = plt.subplots(1, len(PUCKS), figsize=(13 * len(PUCKS), 11),
                             sharey=True)
    colmap = {"w_flat": C_GRAY, "w_share": C_BLUE,
              "w_abs": C_CORAL, "w_argmax": C_MUSTARD}
    for i, p in enumerate(PUCKS):
        ax = axes[i] if len(PUCKS) > 1 else axes
        d = df[df.puck == p]
        for c, lab in WEIGHTS:
            s = d[d.weight == c].sort_values("radius")
            if not len(s):
                continue
            ax.plot(s.radius, s.n_retained / s.n_pairs.clip(lower=1), "-o",
                    lw=3, ms=14, color=colmap[c], label=lab)
        ax.axhline(0.5, color=C_DARK, ls="--", lw=2)
        ax.set_ylim(0, 1.02)
        style(ax, title=f"Puck {p}", xlabel="radius (µm)",
              ylabel="fraction of pairs retained" if i == 0 else None)
        if i == 0:
            ax.legend(fontsize=FS_ANNOT, frameon=False, loc="lower right")
    fig.tight_layout()
    savefig(fig, "Fig_08e_C_pair_retention_vs_radius")

    log()
    log("  Dashed line is 50% retention. Below it, a weighted test A is being")
    log("  run on a different and smaller anchor set than Step08d used, and")
    log("  the two are no longer comparable even informally.")


# ==============================================================================
# BEAT 7
# ==============================================================================
def beat7_design2(adata, radii, mass_floor_pct, pairs, rows):
    banner("BEAT 7  DESIGN 2 PREVIEW: local CD8 / CD4 / Treg, cases vs controls")
    log("  No TCR data, no rarefaction floor, no dependence on the error-")
    log("  derived barcodes. The outcome is the MEAN c2l abundance over beads")
    log("  within R of the anchor. Mean rather than sum, because bead density")
    log("  varies across the puck and a sum would partly measure that.")
    log()
    log(f"  PRIMARY RADIUS IS PRE-SPECIFIED AT {PRIMARY_RADIUS:.0f} um, inherited")
    log("  from Step08d v4. Every other radius here is EXPLORATORY and must be")
    log("  labelled that way downstream. This beat is effectively design 2's")
    log("  statistic, so pre-specifying the radius now is what keeps Step08f")
    log("  from being a second look at the same numbers.")
    log()
    log("  ANCHOR BEAD EXCLUDED from the neighbourhood by default. A T cell")
    log("  touching a 10 um epithelial anchor puts its transcripts on that")
    log("  bead, so the anchor's own CD8 is a direct mediator, which is why")
    log("  Step08c v3 removed cd8 from the matching covariates. Included-anchor")
    log("  values are reported alongside so the size of that effect is visible.")
    log()
    log("  Treg is carried as its OWN outcome, not folded into a denominator.")
    log("  Treg accumulation near neoantigens is a result, not contamination.")

    out = []
    for p in PUCKS:
        f, xy, pos, _ = load_bead_frame(adata, p, mass_floor_pct)
        idx = f.index.values
        tree = cKDTree(xy)
        pp = pairs.get(p, [])
        if len(pp) < 5:
            log(f"\n  --- puck {p}: only {len(pp)} pairs, skipping ---")
            continue
        vals = {c: f[c].values for c in ("cd8", "cd4", "treg", "t_mass", "depth")}

        def local(anchor, R, exclude_self):
            if anchor not in pos:
                return None
            j = pos[anchor]
            near = tree.query_ball_point(xy[j], R)
            if exclude_self:
                near = [k for k in near if k != j]
            if not near:
                return None
            near = np.asarray(near)
            o = {k: float(np.nanmean(v[near])) for k, v in vals.items()}
            o["n_beads"] = len(near)
            tm = o["t_mass"]
            o["cd8_share"] = o["cd8"] / tm if tm > EPS else np.nan
            return o

        for R in radii:
            for excl in (True, False):
                ca, kb = [], []
                for a, b in pp:
                    va, vb = local(a, R, excl), local(b, R, excl)
                    if va is None or vb is None:
                        continue
                    ca.append(va)
                    kb.append(vb)
                if len(ca) < 5:
                    continue
                ca = pd.DataFrame(ca)
                kb = pd.DataFrame(kb)
                tag = "anchor excluded" if excl else "anchor included"
                if R == PRIMARY_RADIUS and excl:
                    log(f"\n  --- puck {p}  R = {R:.0f} um  [{tag}]  "
                        f"*** PRE-SPECIFIED PRIMARY *** ---")
                else:
                    log(f"\n  --- puck {p}  R = {R:.0f} um  [{tag}]  exploratory ---")
                log(f"    {len(ca)} complete pairs | median beads in "
                    f"neighbourhood {np.median(ca.n_beads):.0f}")
                log(f"    {'outcome':<14} {'case med':>10} {'ctrl med':>10} "
                    f"{'diff':>9} {'case>ctrl':>11} {'Wilcoxon':>11} {'sign':>9}")
                for metric in ("cd8", "cd4", "treg", "cd8_share", "depth"):
                    a = ca[metric].values.astype(float)
                    b = kb[metric].values.astype(float)
                    ok = np.isfinite(a) & np.isfinite(b)
                    a, b = a[ok], b[ok]
                    if len(a) < 5:
                        continue
                    try:
                        _, pw = wilcoxon(a, b, alternative="greater")
                    except ValueError:
                        pw = np.nan
                    npos = int((a > b).sum())
                    try:
                        ps = binomtest(npos, len(a), 0.5,
                                       alternative="greater").pvalue
                    except Exception:
                        ps = np.nan
                    log(f"    {metric:<14} {np.median(a):>10.4f} "
                        f"{np.median(b):>10.4f} {np.median(a-b):>+9.4f} "
                        f"{npos:>6}/{len(a):<4} {pw:>11.4g} {ps:>9.4g}")
                    rec = {"beat": "design2_preview", "puck": p, "radius": R,
                           "anchor_excluded": excl, "metric": metric,
                           "n_pairs": len(a),
                           "median_case": float(np.median(a)),
                           "median_control": float(np.median(b)),
                           "median_diff": float(np.median(a - b)),
                           "n_case_gt": npos,
                           "p_wilcoxon": float(pw) if pw == pw else None,
                           "p_sign": float(ps) if ps == ps else None,
                           "prespecified_primary": bool(R == PRIMARY_RADIUS and excl)}
                    rows.append(rec)
                    out.append(rec)
                log("    `depth` is a NEGATIVE CONTROL. If local depth separates")
                log("    cases from controls, the c2l differences above are")
                log("    partly a coverage readout and design 2 needs depth in")
                log("    the model rather than matched at the anchor only.")

    if out:
        d = pd.DataFrame(out)
        d.to_csv(os.path.join(OUT, "design2_preview.tsv"), sep="\t", index=False)

        prim = d[(d.prespecified_primary) & (d.metric.isin(["cd8", "cd4", "treg"]))]
        if len(prim):
            fig, axes = plt.subplots(1, len(PUCKS), figsize=(13 * len(PUCKS), 11),
                                     sharey=False)
            for i, p in enumerate(PUCKS):
                ax = axes[i] if len(PUCKS) > 1 else axes
                s = prim[prim.puck == p]
                if not len(s):
                    style(ax, title=f"Puck {p} (n/a)")
                    continue
                xs = np.arange(len(s))
                ax.bar(xs - 0.2, s.median_control.values, width=0.4,
                       color=C_MUSTARD, label="mutation only")
                ax.bar(xs + 0.2, s.median_case.values, width=0.4,
                       color=C_CORAL, label="neoantigen")
                ax.set_xticks(xs)
                ax.set_xticklabels([m.upper() for m in s.metric.values])
                style(ax, title=f"Puck {p}  (R = {PRIMARY_RADIUS:.0f} µm)",
                      ylabel="median local c2l abundance" if i == 0 else None)
                if i == 0:
                    ax.legend(fontsize=FS_ANNOT, frameon=False)
            fig.tight_layout()
            savefig(fig, "Fig_08e_D_design2_preview_primary_radius")


# ==============================================================================
# BEAT 8
# ==============================================================================
def beat8_design3(adata, radii, mass_floor_pct, cases_by_puck, rows):
    banner("BEAT 8  DESIGN 3 feasibility: anchors per CD8-dominance stratum")
    log("  Split anchors by whether the local T lineage is CD8-dominated, then")
    log("  run the existing tests within stratum. The split is at the puck")
    log("  median of local CD8 share among anchors, so the two strata are")
    log("  balanced by construction and neither is defined by a threshold")
    log("  chosen after seeing an outcome.")
    log(f"  Step08c required {POOL_MIN_CASES} testable anchors per puck. A")
    log("  stratum below that cannot support a test.")

    out = []
    for p in PUCKS:
        f, xy, pos, _ = load_bead_frame(adata, p, mass_floor_pct)
        tree = cKDTree(xy)
        cd8 = f["cd8"].values
        tmass = f["t_mass"].values
        anchors = [a for a in cases_by_puck.get(p, []) if a in pos]
        log(f"\n  --- puck {p}  ({len(anchors)} anchors) ---")
        log(f"    {'R um':>6} {'median share':>14} {'CD8-high':>10} "
            f"{'CD8-low':>9} {'testable?':>11}")
        for R in radii:
            shares = []
            for a in anchors:
                near = tree.query_ball_point(xy[pos[a]], R)
                near = [k for k in near if k != pos[a]]
                if not near:
                    shares.append(np.nan)
                    continue
                near = np.asarray(near)
                tm = float(np.nanmean(tmass[near]))
                shares.append(float(np.nanmean(cd8[near])) / tm if tm > EPS else np.nan)
            shares = np.array(shares, float)
            good = np.isfinite(shares)
            if good.sum() < 4:
                log(f"    {R:>6.0f}  too few anchors with a defined share")
                continue
            med = float(np.nanmedian(shares))
            hi = int((good & (shares > med)).sum())
            lo = int((good & (shares <= med)).sum())
            ok = (hi >= POOL_MIN_CASES) and (lo >= POOL_MIN_CASES)
            log(f"    {R:>6.0f} {med:>14.4f} {hi:>10} {lo:>9} "
                f"{'yes' if ok else 'NO':>11}")
            rows.append({"beat": "design3_strata", "puck": p, "radius": R,
                         "median_local_cd8_share": med, "n_cd8_high": hi,
                         "n_cd8_low": lo, "testable": bool(ok)})
            out.append({"puck": p, "radius": R, "median_share": med,
                        "n_high": hi, "n_low": lo, "testable": ok})
    if out:
        pd.DataFrame(out).to_csv(os.path.join(OUT, "design3_strata.tsv"),
                                 sep="\t", index=False)
    log()
    log("  With 31 / 17 / 53 testable anchors to start from, puck 37 is the one")
    log("  at risk, and puck 37 is where the only nominal test A result lives.")
    log("  If 37 cannot be stratified, design 3 answers a question about the")
    log("  other two pucks only, and should be reported that way.")


# ==============================================================================
# BEAT 9
# ==============================================================================
def beat9_verdict(rows):
    banner("BEAT 9  verdict")
    d = pd.DataFrame(rows)

    def grab(beat):
        return d[d.beat == beat] if "beat" in d.columns else pd.DataFrame()

    co = grab("collinearity")
    if len(co):
        t = co[co.bead_set == "T_cell beads"]
        if len(t):
            worst = float(np.nanmax(np.abs(t["rho_cd8_cd4"].values)))
            log(f"  DESIGN 1 and 3 gate: max |rho(CD8, CD4)| on T_cell beads "
                f"= {worst:.3f}")
            if worst >= 0.8:
                log("    -> the c2l split is one measurement wearing two labels.")
                log("       Do not build a weight on it. Run design 2 only.")
            elif worst >= 0.4:
                log("    -> partial separation. Weighted designs are runnable")
                log("       with the collinearity reported as a limitation.")
            else:
                log("    -> the split carries information. Proceed.")

    at = grab("pair_attrition")
    if len(at):
        s = at[(at.weight == "w_share") & (at.radius == PRIMARY_RADIUS)]
        if len(s):
            log(f"\n  DESIGN 1 gate: pairs retained at R = {PRIMARY_RADIUS:.0f} um "
                f"under w_share")
            for _, r in s.iterrows():
                log(f"    puck {r.puck}: {int(r.n_retained)}/{int(r.n_pairs)} "
                    f"({100*r.frac_retained:.0f}%)")
            best = at[at.weight == "w_share"].groupby("radius")["frac_retained"].mean()
            if len(best):
                log(f"    best radius for retention under w_share: "
                    f"{best.idxmax():.0f} um at {100*best.max():.0f}% mean")

    d2 = grab("design2_preview")
    if len(d2):
        prim = d2[(d2.prespecified_primary == True)]
        if len(prim):
            log(f"\n  DESIGN 2, pre-specified primary (R = {PRIMARY_RADIUS:.0f} um, "
                f"anchor excluded):")
            for _, r in prim.iterrows():
                log(f"    puck {r.puck} {r.metric:<10} diff {r.median_diff:+.4f} "
                    f"| Wilcoxon p {r.p_wilcoxon}")
            log("    These are THREE PUCKS = THREE PATIENTS. Combine evidence,")
            log("    never data. No correction is applied here because this is")
            log("    a diagnostic; Step08f applies it across the family.")

    d3 = grab("design3_strata")
    if len(d3):
        s = d3[d3.radius == PRIMARY_RADIUS]
        if len(s):
            log(f"\n  DESIGN 3 gate at R = {PRIMARY_RADIUS:.0f} um:")
            for _, r in s.iterrows():
                log(f"    puck {r.puck}: {int(r.n_cd8_high)} high / "
                    f"{int(r.n_cd8_low)} low -> "
                    f"{'testable' if r.testable else 'NOT testable'}")

    log()
    log("  STANDING LIMIT, unchanged by anything above: the anchors are still")
    log("  the 110 beads whose neoantigens were predicted against a fixed")
    log("  ten-allele panel rather than patient HLA. A share of them are")
    log("  peptides the patient cannot present. Whatever this analysis shows,")
    log("  it is bounded by that, and the fix is HLA typing from the uBAMs.")


# ==============================================================================
def main():
    global PRIMARY_RADIUS
    ap = argparse.ArgumentParser(
        description="Phase 0 diagnostic for the CD8-restricted reanalysis.")
    ap.add_argument("--chain", default=DEFAULT_CHAIN,
                    help="TRB (default), all, or a comma list")
    ap.add_argument("--radii", type=float, nargs="+", default=RADII)
    ap.add_argument("--primary-radius", type=float, default=PRIMARY_RADIUS)
    ap.add_argument("--mass-floor-pct", type=float, default=MASS_FLOOR_PCT,
                    help="percentile of T-lineage mass among T_cell beads "
                         "below which the CD8 share is treated as zero")
    ap.add_argument("--skip-markers", action="store_true",
                    help="skip BEAT 3, which needs gene columns from X")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()

    PRIMARY_RADIUS = args.primary_radius

    banner("Step08e_CD8_Phase0_Diagnostic")
    log(f"  started   {datetime.now():%Y-%m-%d %H:%M:%S}")
    log(f"  chain     {args.chain}")
    log(f"  radii     {args.radii}")
    log(f"  primary   {PRIMARY_RADIUS:.0f} um  (pre-specified for design 2)")
    log(f"  mass floor percentile {args.mass_floor_pct:.0f}")
    log()
    ok = validate_paths()
    if not ok:
        sys.exit("FATAL: inputs missing, see above.")
    if args.validate_only:
        read_floors(args.chain)
        log("\n--validate-only: paths and floors resolve. Stopping.")
        return

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)
    rows = []

    beat0_backup(not args.no_backup, rows)

    floors = read_floors(args.chain)

    log()
    log("loading annotated object (backed, obs only)")
    adata = ad.read_h5ad(ADATA, backed="r")
    log(f"  {adata.n_obs:,} beads x {adata.n_vars:,} genes")

    log()
    log("loading TCR tables")
    tcr = load_tcr(args.chain)

    log()
    log("loading Step08c design artifacts")
    pr = pd.read_csv(PAIRS, sep="\t")
    pr["puck"] = pr["puck"].astype(str)
    pairs = {p: list(zip(g["case"], g["control"]))
             for p, g in pr.groupby("puck")}
    tc = pd.read_csv(TESTABLE, sep="\t")
    tc["puck"] = tc["puck"].astype(str)
    cases_by_puck = {p: list(g["CB"]) for p, g in tc.groupby("puck")}
    controls_by_puck = {p: [b for _, b in v] for p, v in pairs.items()}
    for p in PUCKS:
        log(f"  puck {p}: {len(pairs.get(p, [])):,} pairs | "
            f"{len(cases_by_puck.get(p, [])):,} testable cases")

    beat1_inventory(adata, args.mass_floor_pct, rows)
    beat2_collinearity(adata, args.mass_floor_pct, tcr, rows)
    beat3_markers(adata, args.mass_floor_pct, tcr, rows, args.skip_markers)
    beat4_weight(adata, args.mass_floor_pct, rows)
    beat5_accounting(adata, tcr, floors, args.radii, args.mass_floor_pct,
                     cases_by_puck, controls_by_puck, rows)
    beat6_attrition(adata, tcr, floors, args.radii, args.mass_floor_pct,
                    pairs, rows)
    beat7_design2(adata, args.radii, args.mass_floor_pct, pairs, rows)
    beat8_design3(adata, args.radii, args.mass_floor_pct, cases_by_puck, rows)
    beat9_verdict(rows)

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "cd8_phase0_summary.tsv"),
                              sep="\t", index=False)
    log(f"\noutputs in {OUT}")
    log(f"  finished  {datetime.now():%Y-%m-%d %H:%M:%S}")


if __name__ == "__main__":
    main()
