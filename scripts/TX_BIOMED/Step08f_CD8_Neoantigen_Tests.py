#!/usr/bin/env python3
"""
Step08f_CD8_Neoantigen_Tests.py
================================================================================
The CD8-conditioned neoantigen tests. Every design decision was fixed in
Step08e BEFORE the outcomes in this script were computed, and the two that were
not are labelled as such in their own output.

WHAT STEP08e SETTLED
    DESIGN 1, CD8-weighted clonality: NOT RUNNABLE. Reported here as an
        infeasibility finding with its retention table, not as a negative.
        Two independent reasons, either sufficient on its own:
          (a) the CD8 share of T lineage rises ~3x across depth deciles
              (Spearman decile vs mean share +0.95 / +0.81 / +0.93), so
              weighting TCR UMIs by it reintroduces the exact depth confound
              the Step08a BEAT 2 rewrite existed to remove;
          (b) paired retention at R = 100 um is 0/30, 0/18, 2/55. The only
              radius clearing 50% is 200 um, where anchors 194-384 um apart
              stop being independent and 08c's 2R pool buffer would exclude
              most of each puck.
    DESIGN 2, CD8 abundance as outcome: RUNS. This script's PRIMARY.
    DESIGN 3, stratified clonality: RUNS, because it stratifies ANCHORS rather
        than reweighting UMIs, so the Step08c floors apply unweighted
        (retention 27/30, 16/18, 49/55 at R = 100).

PRE-SPECIFICATION, AND WHAT IS NOT NEW
    PART A at R = 100 um recomputes exactly what Step08e BEAT 7 reported. Its
    per-puck p-values WILL match that run by construction. This is not
    replication and must never be written up as such. What is new here is the
    cross-puck combination, the depth adjustment, and PARTS B and C.

    Pre-specified in Step08e, before any outcome was seen:
        primary outcome    mean c2l CD8 over beads within 100 um, ANCHOR
                           EXCLUDED (the anchor's own CD8 is a mediator: a T
                           cell touching a 10 um epithelial bead puts its
                           transcripts on that bead, which is why Step08c v3
                           removed cd8 from the matching covariates)
        primary test       paired Wilcoxon, one-sided, within puck
        combination        Stouffer across pucks. THREE PUCKS ARE THREE
                           PATIENTS with private repertoires: combine evidence,
                           never data
        negative control   local depth, same geometry

    Pre-specified in THIS script, before its outcomes were computed:
        PART B annulus     inner rings exceed outer rings if the effect is
                           local. Some spread is expected even for a true point
                           source because Step08b measured bead-to-bead
                           diffusion directly, so the prediction is "inner
                           exceeds outer", NOT "outer is zero"
        PART C             distance-matched controls at 200-600 um, widened if
                           pairs starve. ATTRITION IS REPORTED BEFORE THE TEST

    NOT pre-specified, and labelled exploratory in the output:
        PART E             the conditional Treg / B cell arm. It selects
                           anchors ON THE PRIMARY OUTCOME, so nothing in it is
                           an independent test of CD8. Hypothesis-generating
                           for the next experiment, and the script says so in
                           its own result table.

WHY THE ANNULUS TEST IS THE CONFOUND TEST AND NOT A ROTATION NULL
    A rotation null answers a related question but discards the matched
    controls, which are this design's best feature: cases and controls both
    passed identical SComatic calling filters and differ only in whether the
    mutation yields a predicted binder. Step08d's test C already showed
    rotation machinery produces nothing interpretable on these pucks.

    The annulus decomposition keeps the pairing and separates the two readings
    directly. On a SCALE-FREE basis the Step08e sweep already pointed at a
    local effect: the ratio of the paired difference to the control median
    peaks at or beside 100 um in all three pucks and falls by 200 um. The raw
    paired difference rises monotonically in puck 29 only because the control
    median itself nearly doubles over that range. That is why PART B's
    pre-specified contrast is on the LOG RATIO, not the raw difference.

MULTIPLICITY, SEPARATE FROM THE Step08d FAMILY
    The anchors overlap with 08d's but the hypothesis does not, so these are
    not folded into that BH family.
        primary      1 combined test. No correction.
        confound     PART B contrast + PART C. BH across 2.
        secondary    CD4, Treg, cd8_share. BH across 3 combined p-values.
        design 3     2 strata. BH across 2.
        conditional  none. Exploratory, selected on the primary.

STANDING LIMIT
    The anchors are still the 110 beads whose neoantigens were predicted
    against a fixed ten-allele panel rather than patient HLA, so a share of
    them are peptides the patient cannot present. Every result here is bounded
    by that. The fix is HLA typing from the uBAMs, not more statistics.

    c2l abundance is a MODEL OUTPUT, not a cell count, and Step08e BEAT 4
    showed it is depth-sensitive. "More CD8 abundance" and "more CD8 cells"
    are different claims and the writeup must keep them apart.

USAGE
    python Step08f_CD8_Neoantigen_Tests.py
    python Step08f_CD8_Neoantigen_Tests.py --chain TRB --radius 100
    python Step08f_CD8_Neoantigen_Tests.py --skip-design3
    python Step08f_CD8_Neoantigen_Tests.py --validate-only

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
================================================================================
"""

import argparse
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon, mannwhitneyu, norm, t as tdist, binomtest

try:
    import anndata as ad
except ImportError:
    sys.exit("anndata not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTPUTS = os.path.join(PROOT, "data/outputs")
FEAS = os.path.join(OUTPUTS, "11_neoantigen_tcr")
PHASE0 = os.path.join(OUTPUTS, "13_cd8_phase0")
OUT = os.path.join(OUTPUTS, "14_cd8_tests")
FIGDIR = os.path.join(OUT, "figures")

ADATA = os.path.join(OUTPUTS, "04_annotation/all_pucks_annotated_unified.h5ad")
NEO = os.path.join(OUTPUTS, "07_neoantigen/neoantigens_per_bead.tsv")
MUT = os.path.join(OUTPUTS, "05_mutations/SComatic/SingleCell/mutations_per_bead.tsv")
PAIRS = os.path.join(FEAS, "matched_pairs.tsv")
TCR_DIR = os.path.join(PROOT, "data/inputs/tcr/processed")

PUCKS = ["29", "37", "40"]
TCR_CSV = {p: f"B59_{p}_hTCR_tcr.csv" for p in PUCKS}

ANNOT_COL = "unified_annotation"
BASAL_T = "T_cell"
BASAL_B = "B_cell"
C_CD8 = "c2l_CD8-positive, alpha-beta T cell"
C_CD4 = "c2l_CD4-positive, alpha-beta T cell"
C_TREG = "c2l_regulatory T cell"
C_B_CANDIDATES = ["c2l_B cell", "c2l_B_cell", "c2l_plasma cell"]

PRIMARY_RADIUS = 100.0
DEFAULT_CHAIN = "TRB"

# PART B: four 50 um rings, per Jake. Finer than a two-way inner/outer contrast
# so the falloff SHAPE is visible rather than just its sign.
RING_EDGES = [0.0, 50.0, 100.0, 150.0, 200.0]
RING_INNER = [0, 1]      # 0-50, 50-100     pre-specified as INNER
RING_OUTER = [2, 3]      # 100-150, 150-200 pre-specified as OUTER
MIN_RING_BEADS = 3       # a ring with fewer beads is not a mean, it is a bead

# PART C: distance-matched controls. Far enough that the two neighbourhoods
# cannot overlap at R = 100, close enough to share regional context. Widened
# rather than dropped if pairs starve.
DIST_WINDOWS = [(200.0, 600.0), (200.0, 800.0), (200.0, 1200.0)]
DIST_MIN_PAIRS = 8
DEPTH_CALIPER = 0.25     # log-depth SD, the Step08c matching caliper
SMD_OK = 0.10            # Step08c balance threshold

# PART D: Step08d v4 settings, unchanged, because the stratified tests ARE the
# unweighted 08d tests run inside a stratum.
N_RAREFY = 400           # v4: 100 moved p by 0.24 -> 0.30 on Monte Carlo alone
DESIGN3_MIN_ANCHORS = 8  # Step08c's POOL_MIN_CASES

# PART E
TERTILE = 3

EPS = 1e-12
ALPHA = 0.05

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
C_CORAL = "#ed6a5a"     # cases, neoantigen
C_MUSTARD = "#F6D155"   # controls, mutation only
C_GRAY = "#d3d3d3"
C_DARK = "#4a4a4a"
C_BLUE = "#7eb0d5"      # CD8
C_GREEN = "#b2e061"     # CD4 / B
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
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    spool = np.sqrt((np.var(a, ddof=1) + np.var(b, ddof=1)) / 2.0)
    return float((a.mean() - b.mean()) / spool) if spool > 0 else np.nan


def stouffer(pvals, labels=None):
    """Combine one-sided p-values across pucks, equal weights.

    Three pucks are three patients with private TCR repertoires and private
    mutation sets. Evidence combines; data never does. Equal weights rather
    than n-weights, because the patients are the unit and weighting by bead
    count would let one puck's tissue area speak for another patient.
    """
    p = np.asarray([x for x in pvals if x is not None and np.isfinite(x)], float)
    if not len(p):
        return np.nan, np.nan
    p = np.clip(p, 1e-300, 1 - 1e-16)
    z = norm.isf(p)
    zc = z.sum() / np.sqrt(len(z))
    return float(zc), float(norm.sf(zc))


def bh(pvals):
    """Benjamini-Hochberg, returned in input order. NaNs pass through."""
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    out = np.full(len(p), np.nan)
    if not ok.sum():
        return out
    q = p[ok]
    n = len(q)
    order = np.argsort(q)
    ranked = q[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def paired_test(a, b, name="", quiet=False):
    """One-sided paired Wilcoxon plus the sign test, on complete pairs only."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    if len(a) < 5:
        return {"n": int(len(a)), "p_wilcoxon": np.nan, "p_sign": np.nan,
                "median_case": np.nan, "median_control": np.nan,
                "median_diff": np.nan, "n_case_gt": 0}
    try:
        _, pw = wilcoxon(a, b, alternative="greater")
    except ValueError:
        pw = np.nan
    npos = int((a > b).sum())
    try:
        ps = binomtest(npos, len(a), 0.5, alternative="greater").pvalue
    except Exception:
        ps = np.nan
    return {"n": int(len(a)), "p_wilcoxon": float(pw) if pw == pw else np.nan,
            "p_sign": float(ps) if ps == ps else np.nan,
            "median_case": float(np.median(a)),
            "median_control": float(np.median(b)),
            "median_diff": float(np.median(a - b)), "n_case_gt": npos}


def ols_intercept(y, x):
    """Fit y = b0 + b1*x, return b0 and a one-sided p for b0 > 0.

    This is the depth-adjusted version of the paired test: y is the per-pair
    difference in the outcome, x the per-pair difference in log local depth.
    The intercept is the outcome difference at zero depth difference, which is
    the quantity the unadjusted paired test assumes the matching already
    delivered. Reported ALONGSIDE the unadjusted primary, never instead of it.
    """
    y, x = np.asarray(y, float), np.asarray(x, float)
    ok = np.isfinite(y) & np.isfinite(x)
    y, x = y[ok], x[ok]
    n = len(y)
    if n < 6:
        return np.nan, np.nan, np.nan, n
    X = np.column_stack([np.ones(n), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - 2
    s2 = float(resid @ resid) / dof
    try:
        cov = s2 * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return float(beta[0]), np.nan, np.nan, n
    se0 = float(np.sqrt(cov[0, 0]))
    if not np.isfinite(se0) or se0 <= 0:
        return float(beta[0]), np.nan, np.nan, n
    tstat = beta[0] / se0
    return float(beta[0]), se0, float(tdist.sf(tstat, dof)), n


def log_ratio(case, control):
    """Per-pair log ratio, the scale-free contrast PART B is pre-specified on.

    The raw paired difference is not comparable across rings because the
    control level itself grows with radius: in puck 29 the control median
    nearly doubles from 100 to 200 um. A difference that looks like it is
    growing can simply be the denominator growing. log((case+e)/(control+e))
    removes that and is symmetric about zero.
    """
    c = np.asarray(case, float)
    k = np.asarray(control, float)
    e = 1e-6
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log((c + e) / (k + e))


# ==============================================================================
# loading
# ==============================================================================
def validate_paths():
    ok = True
    need = [("adata", ADATA), ("neoantigens", NEO), ("mutations", MUT),
            ("matched_pairs", PAIRS)]
    for p in PUCKS:
        need.append((f"tcr {p}", os.path.join(TCR_DIR, TCR_CSV[p])))
    for name, path in need:
        hit = os.path.exists(path)
        log(f"  {'OK ' if hit else '***'} {name:<18} {path}")
        ok &= hit
    for name, path in (("phase0 attrition", os.path.join(PHASE0, "pair_attrition.tsv")),):
        log(f"  {'OK ' if os.path.exists(path) else '(-)'} {name:<18} {path}")
    return ok


def read_floors(chain):
    path = os.path.join(FEAS, f"feasibility_summary_{chain}.tsv")
    if not os.path.exists(path):
        sys.exit(f"FATAL: {path} not found. Run Step08c --chain {chain} first. "
                 f"Floors are never hardcoded; that was the v2 bug.")
    d = pd.read_csv(path, sep="\t")
    rows = d[d.beat == "rarefaction_chosen"]
    if not len(rows):
        sys.exit(f"FATAL: no 'rarefaction_chosen' rows in {path}.")
    floors = {str(r.puck): int(r.floor) for r in rows.itertuples()}
    log(f"  floors from feasibility_summary_{chain}.tsv: {floors}")
    log("  NOTE: Step08e read these same values. If they differ from the "
        "20/75/50 in")
    log("  the Step08d docstring, 08c was rerun after 08d and the docstring is "
        "stale.")
    return floors


def load_tcr(chain):
    frames = []
    for p in PUCKS:
        d = pd.read_csv(os.path.join(TCR_DIR, TCR_CSV[p]))
        d["puck"] = p
        d["bead"] = d["bc"].astype(str) + f"-1_Puck_211214_{p}"
        frames.append(d)
    t = pd.concat(frames, ignore_index=True)
    if chain and chain != "all":
        keep = [c.strip() for c in chain.split(",")]
        t = t[t["chain"].isin(keep)].copy()
    log(f"  TCR rows after chain filter [{chain}]: {len(t):,}")
    return t.reset_index(drop=True)


def find_b_column(obs_cols):
    for c in C_B_CANDIDATES:
        if c in obs_cols:
            return c
    hits = [c for c in obs_cols if c.startswith("c2l_B")]
    return hits[0] if hits else None


def puck_frame(adata, p, b_col):
    m = adata.obs["puck_id"].astype(str).values == f"Puck_211214_{p}"
    idx = adata.obs_names[m]
    o = adata.obs.loc[idx]
    f = pd.DataFrame(index=idx)
    f["depth"] = o["total_counts"].values.astype(float)
    f["log_depth"] = np.log(np.clip(f["depth"].values, 1, None))
    f["annot"] = o[ANNOT_COL].astype(str).values
    f["is_B"] = (f["annot"] == BASAL_B).astype(float)
    f["is_T"] = (f["annot"] == BASAL_T).astype(float)
    f["cd8"] = o[C_CD8].values.astype(float)
    f["cd4"] = o[C_CD4].values.astype(float)
    f["treg"] = o[C_TREG].values.astype(float)
    f["t_mass"] = f[["cd8", "cd4", "treg"]].sum(axis=1)
    f["bcell"] = (o[b_col].values.astype(float) if b_col else np.nan)
    xy = o[["x_coord", "y_coord"]].values.astype(float)
    pos = {b: i for i, b in enumerate(idx)}
    return f, xy, pos


METRICS = ["cd8", "cd4", "treg", "bcell", "is_B", "t_mass", "depth", "log_depth"]


def disc_profile(tree, xy, vals, j, R, exclude_self=True):
    """Mean of each metric over beads within R of bead index j."""
    near = tree.query_ball_point(xy[j], R)
    if exclude_self:
        near = [k for k in near if k != j]
    if not near:
        return None
    near = np.asarray(near)
    o = {k: float(np.nanmean(v[near])) for k, v in vals.items()}
    o["n_beads"] = len(near)
    o["cd8_share"] = o["cd8"] / o["t_mass"] if o["t_mass"] > EPS else np.nan
    return o


def ring_profile(tree, xy, vals, j, edges):
    """Mean of each metric inside each annulus, anchor bead always excluded."""
    outer = tree.query_ball_point(xy[j], edges[-1])
    outer = np.asarray([k for k in outer if k != j])
    if not len(outer):
        return [None] * (len(edges) - 1)
    d = np.linalg.norm(xy[outer] - xy[j], axis=1)
    out = []
    for i in range(len(edges) - 1):
        m = (d >= edges[i]) & (d < edges[i + 1])
        if m.sum() < MIN_RING_BEADS:
            out.append(None)
            continue
        sel = outer[m]
        o = {k: float(np.nanmean(v[sel])) for k, v in vals.items()}
        o["n_beads"] = int(m.sum())
        o["cd8_share"] = o["cd8"] / o["t_mass"] if o["t_mass"] > EPS else np.nan
        out.append(o)
    return out


# ==============================================================================
# PART A
# ==============================================================================
def part_a_primary(adata, pairs, b_col, R, rows, store):
    banner("PART A  PRIMARY: local CD8 abundance, cases vs matched controls")
    log(f"  Outcome: mean c2l CD8 over beads within {R:.0f} um, ANCHOR EXCLUDED.")
    log("  Test: one-sided paired Wilcoxon within puck, Stouffer across pucks.")
    log()
    log("  THIS RECOMPUTES Step08e BEAT 7. The per-puck p-values below will")
    log("  MATCH that run by construction. It is not replication and must not")
    log("  be written up as such. New here: the cross-puck combination, the")
    log("  depth adjustment, and PARTS B and C.")

    prim_p, sec = [], {k: [] for k in ("cd4", "treg", "cd8_share")}
    neg_p = []
    for p in PUCKS:
        f, xy, pos = puck_frame(adata, p, b_col)
        tree = cKDTree(xy)
        vals = {k: f[k].values for k in METRICS}
        pp = pairs.get(p, [])
        ca, kb = [], []
        for a, b in pp:
            if a not in pos or b not in pos:
                continue
            va = disc_profile(tree, xy, vals, pos[a], R)
            vb = disc_profile(tree, xy, vals, pos[b], R)
            if va is None or vb is None:
                continue
            ca.append(va)
            kb.append(vb)
        if len(ca) < 5:
            log(f"\n  --- puck {p}: only {len(ca)} complete pairs, skipping ---")
            continue
        ca, kb = pd.DataFrame(ca), pd.DataFrame(kb)
        store[p] = (ca, kb)
        log(f"\n  --- puck {p}  ({len(ca)} complete pairs, median "
            f"{np.median(ca.n_beads):.0f} beads per neighbourhood) ---")
        log(f"    {'outcome':<12} {'role':<11} {'case med':>10} {'ctrl med':>10} "
            f"{'diff':>10} {'case>ctrl':>11} {'Wilcoxon':>11}")
        for metric, role in (("cd8", "PRIMARY"), ("cd4", "secondary"),
                             ("treg", "secondary"), ("cd8_share", "secondary"),
                             ("depth", "neg control")):
            r = paired_test(ca[metric].values, kb[metric].values)
            log(f"    {metric:<12} {role:<11} {r['median_case']:>10.4f} "
                f"{r['median_control']:>10.4f} {r['median_diff']:>+10.4f} "
                f"{r['n_case_gt']:>6}/{r['n']:<4} {r['p_wilcoxon']:>11.4g}")
            rows.append({"part": "A_primary", "puck": p, "radius": R,
                         "metric": metric, "role": role, **r})
            if metric == "cd8":
                prim_p.append(r["p_wilcoxon"])
            elif metric == "depth":
                neg_p.append(r["p_wilcoxon"])
            else:
                sec[metric].append(r["p_wilcoxon"])

        # depth-adjusted intercept, the required sensitivity
        dy = ca["cd8"].values - kb["cd8"].values
        dx = ca["log_depth"].values - kb["log_depth"].values
        b0, se0, pint, nfit = ols_intercept(dy, dx)
        log(f"    depth-adjusted: intercept {b0:+.5f} (se {se0:.5f}) "
            f"one-sided p = {pint:.4g}  [n={nfit}]")
        log("      intercept = CD8 difference at zero depth difference. "
            "Sensitivity,")
        log("      reported alongside the unadjusted primary, never instead "
            "of it.")
        rows.append({"part": "A_depth_adjusted", "puck": p, "radius": R,
                     "metric": "cd8", "intercept": b0, "se": se0,
                     "p_wilcoxon": pint, "n": nfit})

    banner("PART A  combined across pucks")
    z, pc = stouffer(prim_p)
    log(f"  PRIMARY  local CD8, per-puck p = "
        f"{', '.join(f'{x:.4g}' for x in prim_p)}")
    log(f"           Stouffer z = {z:+.3f}, combined one-sided p = {pc:.4g}")
    rows.append({"part": "A_combined", "metric": "cd8", "role": "PRIMARY",
                 "radius": R, "stouffer_z": z, "p_combined": pc,
                 "n_pucks": len(prim_p)})

    zn, pn = stouffer(neg_p)
    log(f"\n  NEGATIVE CONTROL  local depth, per-puck p = "
        f"{', '.join(f'{x:.4g}' for x in neg_p)}")
    log(f"           Stouffer z = {zn:+.3f}, combined p = {pn:.4g}")
    log("    A LARGE p HERE IS THE PASS. It means cases do not have more local")
    log("    coverage than controls. A small p would mean the CD8 result is")
    log("    partly a coverage readout.")
    rows.append({"part": "A_combined", "metric": "depth", "role": "neg control",
                 "radius": R, "stouffer_z": zn, "p_combined": pn,
                 "n_pucks": len(neg_p)})

    sp, sk = [], []
    for k, v in sec.items():
        zz, pp_ = stouffer(v)
        sk.append(k)
        sp.append(pp_)
        rows.append({"part": "A_combined", "metric": k, "role": "secondary",
                     "radius": R, "stouffer_z": zz, "p_combined": pp_,
                     "n_pucks": len(v)})
    adj = bh(sp)
    log(f"\n  SECONDARY family, BH across {len(sk)}:")
    for k, raw, a in zip(sk, sp, adj):
        log(f"    {k:<12} combined p = {raw:.4g} | BH = {a:.4g}")
        for r in rows:
            if r.get("part") == "A_combined" and r.get("metric") == k:
                r["p_bh"] = float(a)
    log()
    log("  Secondary outcomes are DESCRIPTORS, not competitors for the primary")
    log("  alpha. If CD4 and Treg move with CD8, part of what is being seen is")
    log("  general immune proximity rather than CD8-specific recruitment, and")
    log("  PART B is what distinguishes those two readings.")
    return store


# ==============================================================================
# PART B
# ==============================================================================
def part_b_annulus(adata, pairs, b_col, rows):
    banner("PART B  CONFOUND TEST: annulus decomposition")
    log(f"  Four 50 um rings: " +
        ", ".join(f"{RING_EDGES[i]:.0f}-{RING_EDGES[i+1]:.0f}"
                  for i in range(len(RING_EDGES) - 1)) + " um.")
    log("  PRE-SPECIFIED CONTRAST: mean per-pair log ratio over the INNER two")
    log("  rings minus the OUTER two, one-sided paired Wilcoxon, > 0 expected")
    log("  if the effect is local.")
    log()
    log("  LOCAL EFFECT      inner log ratio exceeds outer.")
    log("  REGIONAL CONFOUND inner and outer comparable: the whole region is")
    log("                    elevated and the anchor is incidental to it.")
    log()
    log("  Spread is EXPECTED even for a true point source, because Step08b")
    log("  measured bead-to-bead diffusion directly. The prediction is 'inner")
    log("  exceeds outer', NOT 'outer is zero'. A nonzero outer ring is not")
    log("  evidence against a local effect on its own.")
    log()
    log("  The contrast is on the LOG RATIO because the control level itself")
    log("  grows with radius: puck 29's control median nearly doubles from 100")
    log("  to 200 um, so a raw difference that looks like it is growing can")
    log("  just be its denominator growing.")

    contrast_p, detail = [], []
    for p in PUCKS:
        f, xy, pos = puck_frame(adata, p, b_col)
        tree = cKDTree(xy)
        vals = {k: f[k].values for k in METRICS}
        pp = pairs.get(p, [])
        lr = {i: [] for i in range(len(RING_EDGES) - 1)}
        raw = {i: [] for i in range(len(RING_EDGES) - 1)}
        kept = 0
        for a, b in pp:
            if a not in pos or b not in pos:
                continue
            ra = ring_profile(tree, xy, vals, pos[a], RING_EDGES)
            rb = ring_profile(tree, xy, vals, pos[b], RING_EDGES)
            if any(x is None for x in ra) or any(x is None for x in rb):
                # a pair enters the contrast only if EVERY ring is populated in
                # both members, so inner and outer are compared on identical
                # pairs rather than on two different subsets
                continue
            kept += 1
            for i in range(len(RING_EDGES) - 1):
                lr[i].append(float(log_ratio(ra[i]["cd8"], rb[i]["cd8"])))
                raw[i].append(ra[i]["cd8"] - rb[i]["cd8"])
        log(f"\n  --- puck {p}  ({kept} of {len(pp)} pairs with all four rings "
            f"populated in both members) ---")
        if kept < 5:
            log("    too few complete pairs for the contrast")
            continue
        log(f"    {'ring (um)':<14} {'median log ratio':>18} "
            f"{'median raw diff':>18} {'pairs case>ctrl':>17}")
        for i in range(len(RING_EDGES) - 1):
            v = np.array(lr[i], float)
            w = np.array(raw[i], float)
            lab = f"{RING_EDGES[i]:.0f}-{RING_EDGES[i+1]:.0f}"
            log(f"    {lab:<14} {np.median(v):>+18.4f} {np.median(w):>+18.5f} "
                f"{int((v > 0).sum()):>12}/{len(v):<4}")
            rows.append({"part": "B_ring", "puck": p, "ring": lab,
                         "median_log_ratio": float(np.median(v)),
                         "median_raw_diff": float(np.median(w)),
                         "n_pairs": len(v), "n_pos": int((v > 0).sum())})
            detail.append({"puck": p, "ring": lab, "ring_idx": i,
                           "median_log_ratio": float(np.median(v)),
                           "n_pairs": len(v)})
        inner = np.mean([np.array(lr[i], float) for i in RING_INNER], axis=0)
        outer = np.mean([np.array(lr[i], float) for i in RING_OUTER], axis=0)
        r = paired_test(inner, outer)
        log(f"\n    CONTRAST inner({RING_EDGES[0]:.0f}-{RING_EDGES[2]:.0f}) vs "
            f"outer({RING_EDGES[2]:.0f}-{RING_EDGES[4]:.0f}):")
        log(f"      inner median {r['median_case']:+.4f} | outer median "
            f"{r['median_control']:+.4f} | diff {r['median_diff']:+.4f}")
        log(f"      {r['n_case_gt']}/{r['n']} pairs inner > outer | "
            f"one-sided Wilcoxon p = {r['p_wilcoxon']:.4g}")
        contrast_p.append(r["p_wilcoxon"])
        rows.append({"part": "B_contrast", "puck": p, **r})

    z, pc = stouffer(contrast_p)
    log(f"\n  COMBINED contrast: per-puck p = "
        f"{', '.join(f'{x:.4g}' for x in contrast_p)}")
    log(f"  Stouffer z = {z:+.3f}, combined one-sided p = {pc:.4g}")
    rows.append({"part": "B_combined", "stouffer_z": z, "p_combined": pc,
                 "n_pucks": len(contrast_p)})

    if detail:
        d = pd.DataFrame(detail)
        fig, axes = plt.subplots(1, len(PUCKS), figsize=(13 * len(PUCKS), 11),
                                 sharey=True)
        for i, p in enumerate(PUCKS):
            ax = axes[i] if len(PUCKS) > 1 else axes
            s = d[d.puck == p].sort_values("ring_idx")
            if not len(s):
                style(ax, title=f"Puck {p} (n/a)")
                continue
            cols = [C_BLUE if k in RING_INNER else C_GRAY
                    for k in s.ring_idx.values]
            ax.bar(np.arange(len(s)), s.median_log_ratio.values, color=cols)
            ax.axhline(0, color=C_DARK, lw=2)
            ax.set_xticks(np.arange(len(s)))
            ax.set_xticklabels(s.ring.values, rotation=45, ha="right")
            style(ax, title=f"Puck {p}", xlabel="ring (µm)",
                  ylabel="median log ratio, case / control" if i == 0 else None)
        fig.tight_layout()
        savefig(fig, "Fig_08f_B_annulus_profile")
    return pc


# ==============================================================================
# PART C
# ==============================================================================
def build_control_pool(adata):
    """Mutation-carrying beads with NO predicted binder, the Step08c controls."""
    neo = pd.read_csv(NEO, sep="\t")
    mut = pd.read_csv(MUT, sep="\t")
    ncol = "CB" if "CB" in neo.columns else neo.columns[0]
    mcol = "CB" if "CB" in mut.columns else mut.columns[0]
    cases = set(neo[ncol].astype(str))
    allmut = set(mut[mcol].astype(str))
    controls = allmut - cases
    obs = set(adata.obs_names)
    controls = {b for b in controls if b in obs}
    cases = {b for b in cases if b in obs}
    log(f"  neoantigen beads {len(cases):,} | mutation-only controls "
        f"{len(controls):,} (both present in adata)")
    log("  mutations_per_bead.tsv is NOT epithelial-only: 189 epithelial, 82")
    log("  myeloid, 15 T_cell, 14 B_cell, 12 fibroblast. Controls inherit that")
    log("  composition, same as Step08c.")
    return cases, controls


def part_c_distance_matched(adata, cases, controls, b_col, R, rows):
    banner("PART C  CONFIRMATORY: distance-matched controls")
    log("  Each case is re-matched to a mutation-only control drawn from a")
    log("  DISTANCE ANNULUS: far enough that the two neighbourhoods cannot")
    log("  overlap at the primary radius, close enough to share regional")
    log("  context. Regional structure then cancels WITHIN the pair by")
    log("  construction, which is what makes this the confirmatory arm.")
    log(f"  Depth matched on log depth, caliper {DEPTH_CALIPER} SD, greedy")
    log("  nearest neighbour without replacement, exactly as Step08c BEAT 3b.")
    log()
    log("  ATTRITION IS REPORTED BEFORE ANY TEST. Step08c measured anchors at")
    log("  194-384 um apart, so the window is tight and this arm may starve.")
    log("  Windows are widened rather than dropped, and the window actually")
    log("  used is reported per puck.")

    part_p = []
    for p in PUCKS:
        f, xy, pos = puck_frame(adata, p, b_col)
        tree = cKDTree(xy)
        vals = {k: f[k].values for k in METRICS}
        suf = f"Puck_211214_{p}"
        cs = [b for b in cases if b.endswith(suf) and b in pos]
        ks = [b for b in controls if b.endswith(suf) and b in pos]
        log(f"\n  --- puck {p}  ({len(cs)} cases, {len(ks)} control candidates) ---")
        if len(cs) < 5 or len(ks) < 5:
            log("    too few beads to match")
            continue
        ld = f["log_depth"].values
        sd = float(np.std(np.r_[ld[[pos[b] for b in cs]],
                                ld[[pos[b] for b in ks]]], ddof=1))
        cal = DEPTH_CALIPER * sd

        chosen = None
        for (dmin, dmax) in DIST_WINDOWS:
            used, pairs_c = set(), []
            cand = {}
            for a in cs:
                ja = pos[a]
                d = np.linalg.norm(xy[[pos[b] for b in ks]] - xy[ja], axis=1)
                m = (d >= dmin) & (d <= dmax)
                cand[a] = [(ks[i], abs(ld[pos[ks[i]]] - ld[ja]), d[i])
                           for i in np.where(m)[0]
                           if abs(ld[pos[ks[i]]] - ld[ja]) <= cal]
            # scarcest cases first, so a case with one candidate is not starved
            for a in sorted(cs, key=lambda x: len(cand[x])):
                opts = [o for o in cand[a] if o[0] not in used]
                if not opts:
                    continue
                best = min(opts, key=lambda o: o[1])
                used.add(best[0])
                pairs_c.append((a, best[0], best[2]))
            log(f"    window {dmin:.0f}-{dmax:.0f} um -> {len(pairs_c)} pairs "
                f"of {len(cs)} cases ({100*len(pairs_c)/len(cs):.0f}%)")
            if len(pairs_c) >= DIST_MIN_PAIRS:
                chosen = ((dmin, dmax), pairs_c)
                break
            chosen = ((dmin, dmax), pairs_c)
        (dmin, dmax), pairs_c = chosen
        if len(pairs_c) < DIST_MIN_PAIRS:
            log(f"    *** ARM INFEASIBLE in this puck: {len(pairs_c)} pairs "
                f"after widening to {dmax:.0f} um, below the {DIST_MIN_PAIRS} "
                f"minimum. Reported, not tested.")
            rows.append({"part": "C_attrition", "puck": p, "window_min": dmin,
                         "window_max": dmax, "n_pairs": len(pairs_c),
                         "n_cases": len(cs), "feasible": False})
            continue
        dsep = np.array([d for _, _, d in pairs_c])
        a_ld = np.array([ld[pos[a]] for a, _, _ in pairs_c])
        b_ld = np.array([ld[pos[b]] for _, b, _ in pairs_c])
        s = smd(a_ld, b_ld)
        log(f"    USING window {dmin:.0f}-{dmax:.0f} um | {len(pairs_c)} pairs "
            f"| median separation {np.median(dsep):.0f} um")
        log(f"    log-depth SMD after matching: {s:+.3f} "
            f"({'balanced' if abs(s) < SMD_OK else '*** IMBALANCED ***'})")
        rows.append({"part": "C_attrition", "puck": p, "window_min": dmin,
                     "window_max": dmax, "n_pairs": len(pairs_c),
                     "n_cases": len(cs), "median_separation_um": float(np.median(dsep)),
                     "logdepth_smd": s, "feasible": True})

        ca, kb = [], []
        for a, b, _ in pairs_c:
            va = disc_profile(tree, xy, vals, pos[a], R)
            vb = disc_profile(tree, xy, vals, pos[b], R)
            if va is None or vb is None:
                continue
            ca.append(va)
            kb.append(vb)
        if len(ca) < 5:
            log("    too few complete neighbourhoods to test")
            continue
        ca, kb = pd.DataFrame(ca), pd.DataFrame(kb)
        log(f"    {'outcome':<12} {'case med':>10} {'ctrl med':>10} "
            f"{'diff':>10} {'case>ctrl':>11} {'Wilcoxon':>11}")
        for metric in ("cd8", "cd4", "treg", "depth"):
            r = paired_test(ca[metric].values, kb[metric].values)
            log(f"    {metric:<12} {r['median_case']:>10.4f} "
                f"{r['median_control']:>10.4f} {r['median_diff']:>+10.4f} "
                f"{r['n_case_gt']:>6}/{r['n']:<4} {r['p_wilcoxon']:>11.4g}")
            rows.append({"part": "C_test", "puck": p, "metric": metric,
                         "window_min": dmin, "window_max": dmax, **r})
            if metric == "cd8":
                part_p.append(r["p_wilcoxon"])

    z, pc = stouffer(part_p)
    log(f"\n  COMBINED distance-matched CD8: per-puck p = "
        f"{', '.join(f'{x:.4g}' for x in part_p) if part_p else 'none'}")
    log(f"  Stouffer z = {z:+.3f}, combined one-sided p = {pc:.4g}")
    log(f"  Based on {len(part_p)} of {len(PUCKS)} pucks. If fewer than three,")
    log("  this arm speaks for the pucks it covered and no others.")
    rows.append({"part": "C_combined", "metric": "cd8", "stouffer_z": z,
                 "p_combined": pc, "n_pucks": len(part_p)})
    return pc


# ==============================================================================
# PART D
# ==============================================================================
def rarefied_top_fraction(labels, floor, n_rep=N_RAREFY):
    """Rarefied top-clone fraction, the Step08d test A statistic.

    Rarefaction is not optional: a maximum computed from fewer UMIs is biased
    upward, and on an even repertoire a 20-UMI neighbourhood shows a top-clone
    fraction of 0.100 against 0.020 at 600 UMIs, a fivefold inflation from
    sample size alone.
    """
    n = len(labels)
    if n < floor:
        return None
    tops = np.empty(n_rep)
    for r in range(n_rep):
        s = RNG.choice(labels, floor, replace=False)
        _, cnts = np.unique(s, return_counts=True)
        tops[r] = (cnts / floor).max()
    return float(tops.mean())


def part_d_stratified(adata, tcr, pairs, floors, b_col, R, rows):
    banner("PART D  DESIGN 3: clonality stratified on depth-residualized CD8")
    log("  The tests inside each stratum are the UNWEIGHTED Step08d test A, so")
    log("  the Step08c floors apply unweighted and retention stays usable.")
    log("  That is the whole reason design 3 survived where design 1 did not:")
    log("  it stratifies ANCHORS, it does not reweight UMIs.")
    log()
    log("  STRATIFICATION VARIABLE: local CD8 share at the primary radius,")
    log("  RESIDUALIZED ON LOCAL LOG DEPTH within puck. Step08e BEAT 4 found")
    log("  the CD8 share climbs ~3x across depth deciles, so a raw-share split")
    log("  would be a depth split wearing a CD8 label.")
    log()
    log("  Pairs are stratified by the CASE member's residual. Stratifying the")
    log("  pairs rather than testable_cases.tsv avoids the Step08d membership")
    log("  mismatch: the two anchor sets total 103 each but are not the same")
    log("  anchors.")
    log()
    log("  Test A only. Tests B and D score each anchor against a background")
    log("  pool with BH across anchors within puck; at 8-9 anchors per stratum")
    log("  that correction has no resolution, so running them would be theatre.")

    strat_p = {"CD8_high": [], "CD8_low": []}
    for p in PUCKS:
        f, xy, pos = puck_frame(adata, p, b_col)
        tree = cKDTree(xy)
        vals = {k: f[k].values for k in METRICS}
        floor = floors.get(p)
        tp = tcr[tcr.puck == p].reset_index(drop=True)
        idx_set = set(f.index)
        bead_rows = {bead: sub.index.values
                     for bead, sub in tp[tp.bead.isin(idx_set)].groupby(
                         "bead", sort=False)}
        clone = tp["cloneId"].values
        idx = f.index.values

        def nb_labels(anchor):
            if anchor not in pos:
                return None
            near = tree.query_ball_point(xy[pos[anchor]], R)
            acc = [bead_rows[idx[j]] for j in near if idx[j] in bead_rows]
            return clone[np.concatenate(acc)] if acc else np.array([], int)

        recs = []
        for a, b in pairs.get(p, []):
            if a not in pos or b not in pos:
                continue
            va = disc_profile(tree, xy, vals, pos[a], R)
            if va is None or not np.isfinite(va.get("cd8_share", np.nan)):
                continue
            la, lb = nb_labels(a), nb_labels(b)
            if la is None or lb is None:
                continue
            ta = rarefied_top_fraction(la, floor)
            tb = rarefied_top_fraction(lb, floor)
            recs.append({"case": a, "control": b, "share": va["cd8_share"],
                         "log_depth": va["log_depth"], "top_case": ta,
                         "top_ctrl": tb})
        if len(recs) < 10:
            log(f"\n  --- puck {p}: only {len(recs)} usable pairs, skipping ---")
            continue
        d = pd.DataFrame(recs)
        X = np.column_stack([np.ones(len(d)), d["log_depth"].values])
        beta, *_ = np.linalg.lstsq(X, d["share"].values, rcond=None)
        d["resid"] = d["share"].values - X @ beta
        med = float(np.median(d["resid"].values))
        d["stratum"] = np.where(d["resid"].values > med, "CD8_high", "CD8_low")

        hi = d[d.stratum == "CD8_high"]
        lo = d[d.stratum == "CD8_low"]
        s_depth = smd(hi["log_depth"].values, lo["log_depth"].values)
        s_share = smd(hi["share"].values, lo["share"].values)
        flag = abs(s_depth) >= SMD_OK
        log(f"\n  --- puck {p}  ({len(d)} pairs, floor {floor} UMIs) ---")
        log(f"    residual split at {med:+.4f} -> {len(hi)} CD8-high, "
            f"{len(lo)} CD8-low")
        log(f"    BALANCE  log-depth SMD between strata {s_depth:+.3f} "
            f"{'*** IMBALANCED, FLAGGED, PROCEEDING ***' if flag else '(balanced)'}")
        log(f"             raw CD8 share SMD between strata {s_share:+.3f} "
            f"(this one SHOULD be large; it is the split)")
        rows.append({"part": "D_balance", "puck": p, "n_pairs": len(d),
                     "split_value": med, "n_high": len(hi), "n_low": len(lo),
                     "logdepth_smd": s_depth, "share_smd": s_share,
                     "imbalanced": bool(flag)})

        for name, g in (("CD8_high", hi), ("CD8_low", lo)):
            ok = g.dropna(subset=["top_case", "top_ctrl"])
            n_floor_loss = len(g) - len(ok)
            if len(ok) < DESIGN3_MIN_ANCHORS:
                log(f"    {name:<9} {len(ok)} pairs above floor "
                    f"({n_floor_loss} lost to the floor) -> below the "
                    f"{DESIGN3_MIN_ANCHORS} minimum, NOT TESTED")
                rows.append({"part": "D_test", "puck": p, "stratum": name,
                             "n_pairs": len(ok), "n_lost_to_floor": n_floor_loss,
                             "tested": False})
                continue
            r = paired_test(ok["top_case"].values, ok["top_ctrl"].values)
            log(f"    {name:<9} {r['n']} pairs ({n_floor_loss} lost to floor) | "
                f"case {r['median_case']:.4f} vs ctrl {r['median_control']:.4f} "
                f"| diff {r['median_diff']:+.4f} | {r['n_case_gt']}/{r['n']} | "
                f"p = {r['p_wilcoxon']:.4g}")
            strat_p[name].append(r["p_wilcoxon"])
            rows.append({"part": "D_test", "puck": p, "stratum": name,
                         "n_lost_to_floor": n_floor_loss, "tested": True, **r})

    log()
    ps, ks = [], []
    for name in ("CD8_high", "CD8_low"):
        z, pc = stouffer(strat_p[name])
        ks.append(name)
        ps.append(pc)
        log(f"  {name:<9} per-puck p = "
            f"{', '.join(f'{x:.4g}' for x in strat_p[name]) if strat_p[name] else 'none'}"
            f" | Stouffer z = {z:+.3f}, combined p = {pc:.4g}")
        rows.append({"part": "D_combined", "stratum": name, "stouffer_z": z,
                     "p_combined": pc, "n_pucks": len(strat_p[name])})
    adj = bh(ps)
    log(f"  BH across the 2 strata: " +
        ", ".join(f"{k} {a:.4g}" for k, a in zip(ks, adj)))
    for k, a in zip(ks, adj):
        for r in rows:
            if r.get("part") == "D_combined" and r.get("stratum") == k:
                r["p_bh"] = float(a)
    log()
    log("  READING THIS: a CD8-high stratum that moves while CD8-low does not")
    log("  is the result design 3 exists to find. BOTH strata moving says the")
    log("  stratification is not the operative variable. NEITHER moving is the")
    log("  Step08d result reproduced inside a smaller anchor set, which is")
    log("  what §7 of the handoff predicts: the median bead carries two or")
    log("  three TCR molecules, and at that resolution one cell and two cells")
    log("  are not distinguishable.")


# ==============================================================================
# PART E
# ==============================================================================
def part_e_conditional(store, rows):
    banner("PART E  CONDITIONAL ARM: Treg and B cell at CD8-high anchors")
    log("  *** EXPLORATORY. NOT A TEST OF CD8. ***")
    log("  Anchors are selected on the PRIMARY OUTCOME (top tertile of per-pair")
    log("  CD8 difference), so nothing here is independent of PART A. This is")
    log("  hypothesis-generating for the next experiment and is labelled that")
    log("  way in the result table. No multiplicity correction is applied,")
    log("  because correcting exploratory output implies it was a test.")
    log()
    log("  THE QUESTION: if a neoantigen is enriching CD8, is the local Treg")
    log("  and B compartment distinct there? Co-elevated CD8, Treg and B in one")
    log("  neighbourhood is a tertiary-lymphoid-like signature. The contrast")
    log("  that makes it informative is three-way:")
    log("    (i)  CD8-high cases vs their OWN matched controls (paired)")
    log("    (ii) CD8-high cases vs CD8-LOW cases (unpaired)")
    log("  (i) alone cannot separate 'Tregs follow CD8' from 'Tregs are")
    log("  elevated wherever mutations are'. (ii) is what does that.")

    for p in PUCKS:
        if p not in store:
            continue
        ca, kb = store[p]
        diff = ca["cd8"].values - kb["cd8"].values
        n = len(diff)
        if n < 9:
            log(f"\n  --- puck {p}: {n} pairs, too few for tertiles ---")
            continue
        cut = np.quantile(diff, 1 - 1.0 / TERTILE)
        hi = diff >= cut
        lo = diff < np.quantile(diff, 1.0 / TERTILE)
        log(f"\n  --- puck {p}  ({n} pairs -> {int(hi.sum())} CD8-high, "
            f"{int(lo.sum())} CD8-low) ---")
        log(f"    (i) CD8-high cases vs their own matched controls, paired")
        log(f"      {'outcome':<12} {'case med':>10} {'ctrl med':>10} "
            f"{'diff':>10} {'case>ctrl':>11} {'Wilcoxon':>11}")
        for metric in ("treg", "bcell", "is_B", "cd4"):
            if ca[metric].isna().all():
                continue
            r = paired_test(ca.loc[hi, metric].values, kb.loc[hi, metric].values)
            log(f"      {metric:<12} {r['median_case']:>10.4f} "
                f"{r['median_control']:>10.4f} {r['median_diff']:>+10.4f} "
                f"{r['n_case_gt']:>6}/{r['n']:<4} {r['p_wilcoxon']:>11.4g}")
            rows.append({"part": "E_paired", "puck": p, "metric": metric,
                         "exploratory": True, **r})

        log(f"    (ii) CD8-high cases vs CD8-low cases, unpaired")
        log(f"      {'outcome':<12} {'high med':>10} {'low med':>10} "
            f"{'diff':>10} {'MWU':>11}")
        for metric in ("treg", "bcell", "is_B", "cd4", "cd8"):
            if ca[metric].isna().all():
                continue
            a = ca.loc[hi, metric].values.astype(float)
            b = ca.loc[lo, metric].values.astype(float)
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            if len(a) < 3 or len(b) < 3:
                continue
            try:
                _, pv = mannwhitneyu(a, b, alternative="greater")
            except ValueError:
                pv = np.nan
            log(f"      {metric:<12} {np.median(a):>10.4f} {np.median(b):>10.4f} "
                f"{np.median(a)-np.median(b):>+10.4f} {pv:>11.4g}")
            rows.append({"part": "E_unpaired", "puck": p, "metric": metric,
                         "median_high": float(np.median(a)),
                         "median_low": float(np.median(b)),
                         "n_high": len(a), "n_low": len(b),
                         "p_mwu": float(pv) if pv == pv else None,
                         "exploratory": True})
        log("      NOTE: cd8 is listed in (ii) ONLY to show how large the")
        log("      selection effect is. It is guaranteed significant by")
        log("      construction and carries no information.")


# ==============================================================================
# PART F
# ==============================================================================
def part_f_design1(rows):
    banner("PART F  DESIGN 1: measured infeasibility, carried forward")
    log("  This is NOT a negative result. A negative would be 'CD8-weighted")
    log("  clonality shows no enrichment'. What was found is that the test")
    log("  CANNOT BE RUN at any radius where anchors remain independent.")
    path = os.path.join(PHASE0, "pair_attrition.tsv")
    if not os.path.exists(path):
        log(f"  {path} not found; run Step08e first to regenerate the table.")
        return
    d = pd.read_csv(path, sep="\t")
    s = d[d.weight.isin(["w_flat", "w_share"])].copy()
    s["frac"] = s.n_retained / s.n_pairs.clip(lower=1)
    log(f"\n  Paired retention, Step08e BEAT 6 (methods-ready):")
    log(f"    {'puck':>5} {'R um':>6} {'unweighted':>13} {'CD8-weighted':>15}")
    for p in sorted(s.puck.astype(str).unique()):
        for R in sorted(s[s.puck.astype(str) == p].radius.unique()):
            a = s[(s.puck.astype(str) == p) & (s.radius == R) &
                  (s.weight == "w_flat")]
            b = s[(s.puck.astype(str) == p) & (s.radius == R) &
                  (s.weight == "w_share")]
            if not len(a) or not len(b):
                continue
            log(f"    {p:>5} {R:>6.0f} "
                f"{int(a.n_retained.iloc[0]):>6}/{int(a.n_pairs.iloc[0]):<6} "
                f"{int(b.n_retained.iloc[0]):>8}/{int(b.n_pairs.iloc[0]):<6}")
            rows.append({"part": "F_design1_retention", "puck": p, "radius": R,
                         "n_pairs": int(a.n_pairs.iloc[0]),
                         "retained_unweighted": int(a.n_retained.iloc[0]),
                         "retained_cd8_weighted": int(b.n_retained.iloc[0])})
    log()
    log("  Two independent reasons, either sufficient alone:")
    log("    1. the CD8 share of T lineage rises ~3x across depth deciles")
    log("       (Spearman +0.95 / +0.81 / +0.93), so weighting UMIs by it")
    log("       reintroduces the depth confound Step08a BEAT 2 removed;")
    log("    2. retention above is 0/30, 0/18, 2/55 at 100 um. The only radius")
    log("       clearing 50% is 200 um, where anchors 194-384 um apart stop")
    log("       being independent.")
    log()
    log("  THE USEFUL PART: Step08e BEAT 5 gives the shortfall as a number.")
    log("  Weighted ESS medians at 100 um were 12.7 / 23.5 / 30.0 against")
    log("  floors of 30 / 100 / 75, so roughly 2-4x more usable CD8-attributable")
    log("  TCR depth is needed. That is a target for the realignment pass, not")
    log("  a dead end.")


# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="CD8-conditioned neoantigen tests (designs 2 and 3).")
    ap.add_argument("--chain", default=DEFAULT_CHAIN)
    ap.add_argument("--radius", type=float, default=PRIMARY_RADIUS)
    ap.add_argument("--skip-design3", action="store_true")
    ap.add_argument("--skip-partC", action="store_true")
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()
    R = args.radius

    banner("Step08f_CD8_Neoantigen_Tests")
    log(f"  started  {datetime.now():%Y-%m-%d %H:%M:%S}")
    log(f"  chain    {args.chain}")
    log(f"  radius   {R:.0f} um  (pre-specified in Step08e)")
    log()
    if not validate_paths():
        sys.exit("FATAL: inputs missing, see above.")
    if args.validate_only:
        read_floors(args.chain)
        log("\n--validate-only: paths and floors resolve. Stopping.")
        return

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)
    rows, store = [], {}

    floors = read_floors(args.chain)
    log()
    log("loading annotated object (backed, obs only)")
    adata = ad.read_h5ad(ADATA, backed="r")
    log(f"  {adata.n_obs:,} beads x {adata.n_vars:,} genes")
    b_col = find_b_column(adata.obs.columns)
    if b_col:
        log(f"  B cell abundance column: {b_col}")
    else:
        log("  B cell abundance column: NONE FOUND. PART E falls back to the "
            "local")
        log("  fraction of B_cell-annotated beads, which is a different and "
            "coarser")
        log("  measure; say so if it ends up in the writeup.")

    log()
    log("loading TCR tables")
    tcr = load_tcr(args.chain)

    log()
    log("loading matched pairs")
    pr = pd.read_csv(PAIRS, sep="\t")
    pr["puck"] = pr["puck"].astype(str)
    pairs = {p: list(zip(g["case"], g["control"])) for p, g in pr.groupby("puck")}
    for p in PUCKS:
        log(f"  puck {p}: {len(pairs.get(p, [])):,} pairs")

    store = part_a_primary(adata, pairs, b_col, R, rows, store)
    p_b = part_b_annulus(adata, pairs, b_col, rows)

    p_c = np.nan
    if not args.skip_partC:
        cases, controls = build_control_pool(adata)
        p_c = part_c_distance_matched(adata, cases, controls, b_col, R, rows)

    banner("CONFOUND FAMILY  BH across PART B and PART C")
    adj = bh([p_b, p_c])
    log(f"  PART B annulus contrast   combined p = {p_b:.4g} | BH = {adj[0]:.4g}")
    log(f"  PART C distance-matched   combined p = {p_c:.4g} | BH = {adj[1]:.4g}")
    rows.append({"part": "confound_family", "metric": "B_contrast",
                 "p_combined": p_b, "p_bh": float(adj[0])})
    rows.append({"part": "confound_family", "metric": "C_distance_matched",
                 "p_combined": p_c, "p_bh": float(adj[1])})

    if not args.skip_design3:
        part_d_stratified(adata, tcr, pairs, floors, b_col, R, rows)
    part_e_conditional(store, rows)
    part_f_design1(rows)

    # primary figure
    if store:
        fig, axes = plt.subplots(1, len(PUCKS), figsize=(13 * len(PUCKS), 12),
                                 sharey=False)
        for i, p in enumerate(PUCKS):
            ax = axes[i] if len(PUCKS) > 1 else axes
            if p not in store:
                style(ax, title=f"Puck {p} (n/a)")
                continue
            ca, kb = store[p]
            for x1, x2 in zip(kb["cd8"].values, ca["cd8"].values):
                ax.plot([0, 1], [x1, x2], color=C_GRAY, alpha=0.5, lw=1.5)
            ax.scatter(np.zeros(len(kb)), kb["cd8"].values, s=90,
                       color=C_MUSTARD, zorder=3, label="mutation only")
            ax.scatter(np.ones(len(ca)), ca["cd8"].values, s=90,
                       color=C_CORAL, zorder=3, label="neoantigen")
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["control", "case"])
            style(ax, title=f"Puck {p}",
                  ylabel=f"local c2l CD8, {R:.0f} µm" if i == 0 else None)
            if i == 0:
                ax.legend(fontsize=FS_ANNOT, frameon=False)
        fig.tight_layout()
        savefig(fig, "Fig_08f_A_primary_cd8_paired")

    pd.DataFrame(rows).to_csv(os.path.join(OUT, f"cd8_test_results_{args.chain}.tsv"),
                              sep="\t", index=False)
    log(f"\noutputs in {OUT}")
    log(f"  finished {datetime.now():%Y-%m-%d %H:%M:%S}")


if __name__ == "__main__":
    main()
