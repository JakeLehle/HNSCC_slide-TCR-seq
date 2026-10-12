#!/usr/bin/env python3
"""
Step08d_Neoantigen_TCR_Tests.py
================================================================================
The neoantigen-anchored TCR tests. Every setting was fixed in Step08c BEFORE
any outcome was computed. Nothing here is tuned to the result.

TEST ORDER: A, then B, then D. C is independent.
    The background pool excludes beads near CASE anchors but NOT beads near
    controls, because controls are null by hypothesis. Test A is what CHECKS
    that assumption. If mutation-only beads do attract T cells, the pool is
    contaminated and B is conservative. So A must be read first.

    A  local clonality contrast. Rarefied top-clone fraction, each case against
       its matched control. Paired Wilcoxon signed-rank is primary (the pairs
       are real; discarding them for Mann-Whitney costs power). MWU reported.
    B  PRIMARY. Per-anchor maximum clone enrichment against that clone's
       puck-wide UMI frequency, scored against the background pool, then
       Benjamini-Hochberg across anchors within puck.
    D  in-situ expansion. B restricted to clones spanning >= k distinct beads
       inside the neighbourhood, so one deeply covered bead cannot pass as a
       locally expanding population. k = 2 primary, k = 3 strict.
    C  neoantigen potential field. Kernel-weighted neoantigen density at every
       epithelial bead, regressed on local clonality with log depth as a
       covariate. Null by ROTATION, which preserves the spatial autocorrelation
       of both layers while destroying their alignment.

WHY NOT A MAX STATISTIC WITH A PERMUTATION NULL
    An earlier draft used exactly that and it is INVALID: whichever bead holds
    the largest value lands in the permuted case set half the time regardless
    of its true label, so the p-value floors near 0.5. One anchor shifted +1.00
    returned p = 0.477, type I 0.000, power 0.000. Bootstrapping the null from
    the ~30 control values is worse: a bootstrap maximum is capped at the
    largest control observed, giving type I 0.472. The per-anchor design below
    is calibrated at 0.030-0.046 (Step08c BEAT 6).

TWO ANCHOR SETS, DELIBERATELY NOT HARMONISED
    A uses matched_pairs.tsv  (30 / 18 / 55 = 103)
    B and D use testable_cases.tsv (31 / 18 / 54 = 103)
    The totals coincide by accident; the membership does not. Each test reports
    its own set explicitly so no one compares them assuming they are the same
    anchors.

FIXED SETTINGS (Step08c)
    TCR unit        ALL TCR-bearing beads, never source beads: using source
                    beads assumes one cell per clone, the hypothesis under test
    chain           TRB primary (TRB/TRA = 2.77, TRA is the weaker measurement)
    radii           75 um primary, 100 um reported
    floors          per puck, 29 -> 20, 37 -> 75, 40 -> 50 UMIs
    rarefaction     100 replicates; a maximum from fewer UMIs is biased upward
    enrichment      denominator is the clone's puck-wide UMI frequency, never
                    zero: the top clone sits on 3.5-4.1% of TCR beads and lands
                    on ~6-7 beads of any 75 um neighbourhood by chance
    multiplicity    BH across anchors within puck; Simes to a puck-level p;
                    Stouffer across pucks (three patients, private repertoires,
                    so evidence is combined, never data)

PRE-REGISTERED BOUND (Step08c BEAT 6)
    DETECTABLE  shift >= 0.20 in top-clone fraction (one clone holding ~28% of
                local UMIs against a typical 8%). Power 0.93-1.00 at ONE anchor.
    MARGINAL    0.10-0.20. Power 0.14-0.17 at one anchor.
    NOT VISIBLE < 0.10.
    A null result is bounded at 0.20, NOT evidence that neoantigens do not
    recruit T cells.

USAGE
    python Step08d_Neoantigen_TCR_Tests.py
    python Step08d_Neoantigen_TCR_Tests.py --chain all
    python Step08d_Neoantigen_TCR_Tests.py --radii 75 100 150
    python Step08d_Neoantigen_TCR_Tests.py --skip-C
    python Step08d_Neoantigen_TCR_Tests.py --validate-only

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
================================================================================
"""

import argparse
import os
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon, norm, spearmanr

try:
    import anndata as ad
except ImportError:
    sys.exit("anndata not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTPUTS = os.path.join(PROOT, "data/outputs")
FEAS = os.path.join(OUTPUTS, "11_neoantigen_tcr")
OUT = os.path.join(OUTPUTS, "12_neoantigen_tcr_tests")
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
EPITHELIAL = "epithelial"

# v2: floors are READ from the Step08c summary for the SAME CHAIN, never
# hardcoded. v1 hardcoded 20/75/50, which Step08c had chosen on chain=all
# (269,594 rows) while 08d ran chain=TRB (198,037 rows, 73.5%). The same
# neighbourhoods therefore held ~26% fewer UMIs than the floors were sized for,
# and pair attrition at R=75 hit 8 of 18 in puck 37 against the 94-95% BEAT 4
# predicted. Reading the floors makes that mismatch structurally impossible.
FLOORS_FALLBACK = {"29": 20, "37": 75, "40": 50}
RADII = [100.0, 75.0]      # 100 um is PRIMARY: 75 um is floor-limited in every puck
K_VALUES = [2, 3]
N_RAREFY = 100

# v2: a single pre-specified lambda. v1 combined three lambdas per puck with
# Stouffer at n=9, treating three views of the same data on the same puck as
# independent. They are strongly correlated, so that overstated the evidence.
KERNEL_LAMBDA = 100.0                       # pre-specified, middle of the range
KERNEL_LAMBDAS_REPORT = [50.0, 100.0, 200.0]   # reported, not combined
N_ROT = 1000
ALPHA = 0.05

# v2: test B's statistic is now the rarefied TOP-CLONE FRACTION.
# v1 used max-over-clones enrichment with no restriction, which is degenerate:
# a clone seen once puck-wide has frequency 1/N, so one copy in a sample of
# `floor` scores N/floor, while a genuinely expanded clone at 3.5% puck-wide
# holding 30% of the neighbourhood scores 8.6. The rarest clone always wins and
# its value is a fixed ceiling. Simulating 500 deliberately different
# neighbourhoods returned exactly ONE distinct value, producing mass ties and
# p = 1 in five of six cells. A minimum LOCAL abundance does not fix it: a clone
# with m copies locally and m puck-wide scores N/floor for any m.
# Top-clone fraction has no denominator, cannot saturate, and is the statistic
# the pre-registered bound is written in.
#
# The enrichment statistic survives in test D, where it is restricted to clones
# that are abundant enough puck-wide for the denominator to be estimated.
# v3: LOOSENED. At 10 / 3 the restriction over-killed the test: at R=100, 100%
# of puck 29 anchors and 86-94% of puck 40 anchors had NO qualifying clone, so
# the statistic was 0 nearly everywhere and BH still divided by all anchors.
# These were guesses, not derived. The script now reports the qualifying rate,
# so if it is still high at 5 / 2 the restriction, not the depth, is the problem.
MIN_PUCK_UMI = 5
MIN_LOCAL_UMI = 2
MIN_POOL = 1000        # below this the resolution floor is too coarse to trust

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
C_CORAL = "#ed6a5a"
C_MUSTARD = "#F6D155"
C_GRAY = "#d3d3d3"
C_DARK = "#4a4a4a"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

RNG = np.random.default_rng(20261009)


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


def bh_adjust(pv):
    pv = np.asarray(pv, float); n = len(pv)
    if n == 0:
        return pv
    o = np.argsort(pv); adj = np.empty(n)
    adj[o] = np.minimum.accumulate((pv[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(adj, 1.0)


def simes(pv):
    """Puck-level p for the global null 'no anchor is enriched'. Valid under
    positive dependence, which per-anchor tests against a shared pool have."""
    pv = np.sort(np.asarray(pv, float)); n = len(pv)
    return float(np.min(pv * n / np.arange(1, n + 1))) if n else np.nan


def stouffer(pvals):
    pv = [p for p in pvals if p is not None and p == p and 0 < p < 1]
    if not pv:
        return np.nan
    z = np.sum([norm.isf(p) for p in pv]) / np.sqrt(len(pv))
    return float(norm.sf(z))


def load_floors(chain):
    """Floors come from the Step08c run for THIS chain. Hardcoding them is what
    produced the v1 mismatch, so a missing file is a hard stop, not a fallback."""
    f = os.path.join(FEAS, f"feasibility_summary_{chain}.tsv")
    if not os.path.exists(f):
        log(f"  MISSING {f}")
        log(f"  Step08c has not been run with --chain {chain}. The rarefaction")
        log(f"  floors must come from the SAME chain set as the data, or")
        log(f"  neighbourhoods are compared against a floor sized for a")
        log(f"  different number of UMIs.")
        log(f"    sbatch Step08c_Neoantigen_TCR_Feasibility.sh --chain {chain} "
            f"--radius 100")
        return None
    d = pd.read_csv(f, sep="\t")
    sel = d[d.get("beat", pd.Series(dtype=str)) == "rarefaction_chosen"]
    if not len(sel):
        log(f"  {f} has no rarefaction_chosen rows"); return None
    out = {str(r.puck): int(r.floor) for r in sel.itertuples()}
    log(f"  floors read from {os.path.basename(f)}: {out}")
    return out


def load_tcr(p, chain):
    d = pd.read_csv(os.path.join(TCR_DIR, TCR_CSV[p]), sep=None, engine="python")
    d.columns = [c.strip() for c in d.columns]
    if chain and chain != "all":
        d = d[d.chain == chain]
    d["bead"] = d["bc"].astype(str) + f"-1_Puck_211214_{p}"
    return d


def puck_view(adata, p):
    m = adata.obs["puck_id"].astype(str).values == f"Puck_211214_{p}"
    idx = list(adata.obs_names[m])
    xy = adata.obs.loc[idx, ["x_coord", "y_coord"]].values.astype(float)
    f = pd.DataFrame(index=idx)
    f["depth"] = adata.obs.loc[idx, "total_counts"].values.astype(float)
    f["annot"] = adata.obs.loc[idx, ANNOT_COL].astype(str).values
    return idx, f, xy, {b: i for i, b in enumerate(idx)}


# ==============================================================================
# neighbourhood machinery
# ==============================================================================
def build_index(tcr_p, idx_set):
    """bead -> array of cloneIds, one entry per UMI row. Pre-building this makes
    per-anchor neighbourhood assembly a concatenate instead of a pandas filter,
    which matters because the background pool is ~3,000 anchors per puck."""
    sub = tcr_p[tcr_p.bead.isin(idx_set)]
    out = {}
    for b, g in sub.groupby("bead", sort=False):
        out[b] = g["cloneId"].values
    return out


def neighbourhood(tree, xy, idx, pos, bead_clones, anchor, radius):
    if anchor not in pos:
        return None, None
    near = tree.query_ball_point(xy[pos[anchor]], radius)
    chunks, cb = [], defaultdict(int)
    for j in near:
        b = idx[j]
        arr = bead_clones.get(b)
        if arr is None:
            continue
        chunks.append(arr)
        for c in np.unique(arr):          # distinct beads per clone
            cb[c] += 1
    if not chunks:
        return np.array([], dtype=object), {}
    return np.concatenate(chunks), dict(cb)


def rarefied(labels, floor, puck_freq, puck_umi, cb, k_list, n_rep=N_RAREFY):
    """Subsample to `floor` UMIs, then compute the statistics. Rarefaction is
    not optional: a maximum computed from fewer UMIs is biased upward, and on an
    even repertoire a 20-UMI neighbourhood shows a top-clone fraction of 0.100
    against 0.020 at 600 UMIs, a fivefold inflation from sample size alone."""
    n = len(labels)
    if n < floor:
        return None
    tops = np.empty(n_rep)
    top3 = np.empty(n_rep)
    simp = np.empty(n_rep)
    enrk = {k: np.empty(n_rep) for k in k_list}
    qual = {k: np.empty(n_rep) for k in k_list}
    for r in range(n_rep):
        s = RNG.choice(labels, floor, replace=False)
        vals, cnts = np.unique(s, return_counts=True)
        fr = cnts / floor
        tops[r] = fr.max()
        # secondary metrics: top-clone fraction is a MAXIMUM, so one odd clone
        # moves it. If two clones respond, only the larger is seen.
        top3[r] = np.sort(fr)[-3:].sum()
        simp[r] = float((fr ** 2).sum())     # concentration: higher = more clonal
        e = fr / np.array([puck_freq.get(v, 1e-9) for v in vals])
        # restrict to clones abundant enough LOCALLY and PUCK-WIDE for the
        # enrichment ratio to mean anything (see MIN_* notes above)
        ok = (cnts >= MIN_LOCAL_UMI) & np.array(
            [puck_umi.get(v, 0) >= MIN_PUCK_UMI for v in vals])
        for k in k_list:
            m = ok & np.array([cb.get(v, 0) >= k for v in vals])
            enrk[k][r] = e[m].max() if m.any() else 0.0
            qual[k][r] = m.sum()
    out = {"n_umi": n, "top_frac": float(tops.mean()),
           "top3_frac": float(top3.mean()), "simpson": float(simp.mean())}
    for k in k_list:
        out[f"max_enr_k{k}"] = float(enrk[k].mean())
        out[f"n_qualifying_k{k}"] = float(qual[k].mean())
    return out


def winning_clone(labels, puck_freq, puck_umi, cb, k=None):
    """The clone behind the statistic, on the FULL neighbourhood, for reporting."""
    if not len(labels):
        return None
    vals, cnts = np.unique(labels, return_counts=True)
    fr = cnts / cnts.sum()
    e = fr / np.array([puck_freq.get(v, 1e-9) for v in vals])
    if k is not None:
        m = (np.array([cb.get(v, 0) >= k for v in vals]) &
             (cnts >= MIN_LOCAL_UMI) &
             np.array([puck_umi.get(v, 0) >= MIN_PUCK_UMI for v in vals]))
        if not m.any():
            return None
        e = np.where(m, e, -1.0)
    i = int(np.argmax(e))
    return {"clone": int(vals[i]), "clone_umis": int(cnts[i]),
            "clone_beads": int(cb.get(vals[i], 0)), "enrichment": float(e[i])}


def score_anchors(anchors, tree, xy, idx, pos, bead_clones, floor, puck_freq,
                  puck_umi, k_list):
    rows = []
    for a in anchors:
        labels, cb = neighbourhood(tree, xy, idx, pos, bead_clones, a, R_CUR)
        if labels is None:
            continue
        st = rarefied(labels, floor, puck_freq, puck_umi, cb, k_list)
        if st is None:
            rows.append({"anchor": a, "below_floor": True,
                         "n_umi": 0 if labels is None else len(labels)})
            continue
        st.update({"anchor": a, "below_floor": False})
        rows.append(st)
    return pd.DataFrame(rows)


# ==============================================================================
def preflight(adata, neo, pairs, testable, rows):
    banner("PREFLIGHT  anchor accounting: what was dropped, and why")
    log("  BREADCRUMB. Three filters reduce 110 neoantigen beads to the sets the")
    log("  tests actually use. Each is OUTCOME-BLIND: none looks at TCR data, so")
    log("  none can bias toward a positive result. Dropping anchors always")
    log("  invites that question, so the full accounting is printed here.")
    log()
    log("    1. CELL TYPE   4 of 110 are not epithelial (2 B_cell, 1 myeloid,")
    log("       1 endothelial). Controls are restricted to epithelial, so cases")
    log("       are too, or the contrast is not like-for-like.")
    log("    2. DEPTH MATCHING (test A set). A case with no control inside the")
    log("       0.25 SD log-depth caliper has no partner. 3 dropped.")
    log("    3. BACKGROUND POOL (test B/D set). A case whose neighbourhood depth")
    log("       has no comparable background anywhere on the puck cannot be")
    log("       scored fairly; and one further case was dropped to bring puck")
    log("       37's pool into depth balance.")

    ann = adata.obs[ANNOT_COL].astype(str)
    dep = adata.obs["total_counts"]
    all_cb = set(neo["CB"].astype(str))
    epi = {b for b in all_cb if ann.get(b, "") == EPITHELIAL}
    pair_cases = set(pairs["case"].astype(str))
    test_cases = set(testable["CB"].astype(str))

    info = neo.drop_duplicates("CB").set_index("CB")
    log(f"\n  110 neoantigen beads -> {len(epi)} epithelial -> "
        f"A set {len(pair_cases)}, B/D set {len(test_cases)}")

    log(f"\n  {'bead':<34} {'gene':<10} {'hgvs_p':<16} {'depth':>8}  excluded from")
    for b in sorted(all_cb):
        in_a, in_b = b in pair_cases, b in test_cases
        if in_a and in_b:
            continue
        if b not in epi:
            why = f"both (cell type: {ann.get(b, '?')})"
        elif not in_a and not in_b:
            why = "both"
        elif not in_a:
            why = "A only (no depth-matched control)"
        else:
            why = "B/D only (no comparable background)"
        g = info.loc[b, "gene"] if b in info.index else "?"
        h = info.loc[b, "hgvs_p"] if b in info.index else "?"
        d = dep.get(b, float("nan"))
        log(f"  {b:<34} {str(g):<10} {str(h):<16} {d:>8,.0f}  {why}")
        rows.append({"beat": "preflight_excluded", "bead": b, "gene": g,
                     "hgvs_p": h, "depth": float(d), "excluded_from": why})

    log(f"\n  per-puck anchor counts")
    log(f"    {'puck':>5} {'epithelial':>11} {'A set':>7} {'B/D set':>9}")
    for p in PUCKS:
        sfx = f"Puck_211214_{p}"
        log(f"    {p:>5} {sum(1 for b in epi if b.endswith(sfx)):>11} "
            f"{sum(1 for b in pair_cases if b.endswith(sfx)):>7} "
            f"{sum(1 for b in test_cases if b.endswith(sfx)):>9}")
    log("\n  The two sets are NOT harmonised, by decision. Totals coincide at 103")
    log("  but membership differs. Each test reports its own set so the two are")
    log("  never compared as if they covered the same anchors.")


# ==============================================================================
def test_A(ct, kt, pairs_p, rows, p, R, chain):
    log(f"    TEST A  clonality contrast, rarefied top-clone fraction")
    c = ct[~ct.below_floor].set_index("anchor")
    k = kt[~kt.below_floor].set_index("anchor")
    pairs_ok = [(a, b) for a, b in pairs_p if a in c.index and b in k.index]
    log(f"      anchor set: matched pairs, {len(pairs_p)} total, "
        f"{len(pairs_ok)} with both members above the floor")
    if len(pairs_ok) < 5:
        log("      too few complete pairs to test")
        return None
    ca = np.array([c.loc[a, "top_frac"] for a, _ in pairs_ok])
    kb = np.array([k.loc[b, "top_frac"] for _, b in pairs_ok])
    try:
        _, pw = wilcoxon(ca, kb, alternative="greater")
    except ValueError:
        pw = np.nan
    try:
        _, pm = mannwhitneyu(c["top_frac"].values, k["top_frac"].values,
                             alternative="greater")
    except ValueError:
        pm = np.nan
    log(f"      case median {np.median(ca):.4f} | control median {np.median(kb):.4f} "
        f"| diff {np.median(ca - kb):+.4f}")
    log(f"      paired Wilcoxon (PRIMARY) p = {pw:.4g} | Mann-Whitney p = {pm:.4g}")
    rows.append({"test": "A_clonality", "puck": p, "radius": R, "chain": chain,
                 "n_pairs": len(pairs_ok), "median_case": float(np.median(ca)),
                 "median_control": float(np.median(kb)),
                 "median_diff": float(np.median(ca - kb)),
                 "p_wilcoxon": float(pw), "p_mwu": float(pm)})
    return {"p": float(pw), "ca": ca, "kb": kb}


def test_B(ct, pt, rows, p, R, chain, col, label, k=None):
    tag = (f"k>={k}, >= {MIN_LOCAL_UMI} local / >= {MIN_PUCK_UMI} puck-wide UMIs"
           if k else "rarefied top-clone fraction")
    log(f"    {label}  per-anchor enrichment vs background pool  ({tag})")
    c = ct[~ct.below_floor]
    pool = pt[~pt.below_floor][col].values
    log(f"      anchor set: testable cases, {len(ct)} total, "
        f"{len(c)} above the floor | pool {len(pool):,} above the floor")
    if len(c) < 3 or len(pool) < 200:
        log("      insufficient anchors or pool")
        return None
    if len(pool) < MIN_POOL:
        log(f"      *** pool is only {len(pool):,} (< {MIN_POOL}). Resolution floor")
        log(f"          {1/(len(pool)+1):.2e} is too coarse; treat this arm as")
        log(f"          underpowered and do not combine it as if comparable. ***")
    vals = c[col].values
    if k is not None:
        # An anchor with a zero statistic has NO qualifying clone: it carries no
        # information, and leaving it in only inflates the BH denominator.
        keep = vals > 0
        n_zero = int((~keep).sum())
        if n_zero:
            log(f"      {n_zero}/{len(vals)} anchors have a zero statistic (no")
            log(f"        qualifying clone) and are excluded from the BH family")
        c = c[keep]; vals = vals[keep]
        if len(vals) < 3:
            log("      too few informative anchors after exclusion")
            return None
    pv = np.array([(1 + np.sum(pool >= v)) / (len(pool) + 1) for v in vals])
    adj = bh_adjust(pv)
    sig = adj < ALPHA
    pk = simes(pv)
    log(f"      resolution floor 1/(pool+1) = {1/(len(pool)+1):.2e}")
    log(f"      min raw p {pv.min():.3g} | min BH p {adj.min():.3g} | "
        f"hits {int(sig.sum())}/{len(c)}")
    log(f"      puck-level p (Simes, global null) = {pk:.4g}")
    res = pd.DataFrame({"anchor": c["anchor"].values, "stat": vals,
                        "p": pv, "p_bh": adj, "hit": sig})
    rows.append({"test": label, "puck": p, "radius": R, "chain": chain,
                 "k": k, "n_anchors": len(c), "pool": len(pool),
                 "min_p": float(pv.min()), "min_p_bh": float(adj.min()),
                 "n_hits": int(sig.sum()), "p": pk})
    return {"p": pk, "res": res}



def anchor_independence(cases, xy, pos, R, p, rows):
    """Adjacent anchors give OVERLAPPING neighbourhoods, so they are not
    independent observations. BH and Simes stay valid under positive dependence,
    but the EFFECTIVE number of anchors is lower than the count, and both tests
    become conservative. Reported so the count is never read as an n."""
    pts = np.array([xy[pos[a]] for a in cases if a in pos])
    if len(pts) < 2:
        return
    t = cKDTree(pts)
    pairs = t.query_pairs(2 * R)          # neighbourhoods touch at all
    over = t.query_pairs(R)               # each contains the other's centre
    d = t.query(pts, k=2)[0][:, 1]
    # single-linkage clusters of anchors whose neighbourhoods touch
    seen, clusters = set(), 0
    adj = defaultdict(set)
    for i, j in pairs:
        adj[i].add(j); adj[j].add(i)
    for i in range(len(pts)):
        if i in seen:
            continue
        clusters += 1
        stack = [i]
        while stack:
            u = stack.pop()
            if u in seen:
                continue
            seen.add(u)
            stack.extend(adj[u] - seen)
    log(f"      anchor independence: {len(pts)} anchors | nearest-neighbour "
        f"median {np.median(d):,.0f} um")
    log(f"        pairs within {2*R:.0f} um (overlapping): {len(pairs)} | "
        f"within {R:.0f} um (centres inside): {len(over)}")
    log(f"        effective independent groups: {clusters} of {len(pts)}")
    rows.append({"test": "anchor_independence", "puck": p, "radius": R,
                 "n_anchors": len(pts), "median_nn_um": float(np.median(d)),
                 "overlapping_pairs": len(pairs), "effective_groups": clusters})


def test_E_dose_response(ct, neo_info, xy, pos, cases, R, p, chain, rows):
    """If the association is real it should SCALE with antigen strength, not
    merely with presence. A binary case/control difference of 0.007 is easy to
    dismiss; a monotone trend across predicted binding affinity is not.
    Artifacts rarely produce dose-response."""
    log(f"    TEST E  dose-response")
    c = ct[~ct.below_floor].set_index("anchor")
    if len(c) < 8:
        log(f"      only {len(c)} anchors above floor, too few"); return

    rec = []
    for a in c.index:
        if a not in neo_info.index:
            continue
        ic50 = float(neo_info.loc[a, "mut_ic50"])
        strong = bool(neo_info.loc[a, "is_strong"])
        # local neoantigen burden: how many OTHER anchors sit within R
        burden = sum(1 for b in cases
                     if b != a and b in pos and a in pos
                     and np.hypot(*(xy[pos[b]] - xy[pos[a]])) <= R)
        rec.append({"anchor": a, "top_frac": float(c.loc[a, "top_frac"]),
                    "simpson": float(c.loc[a, "simpson"]),
                    "neg_log_ic50": -np.log10(max(ic50, 1e-3)),
                    "is_strong": strong, "burden": burden})
    if len(rec) < 8:
        log("      too few anchors with affinity annotation"); return
    d = pd.DataFrame(rec)

    for xcol, lab in (("neg_log_ic50", "binding strength (-log10 IC50)"),
                      ("burden", f"neoantigen burden within {R:.0f} um")):
        if d[xcol].nunique() < 3:
            log(f"      {lab}: too few distinct values"); continue
        rho, pv = spearmanr(d[xcol].values, d["top_frac"].values)
        rho_s, pv_s = spearmanr(d[xcol].values, d["simpson"].values)
        log(f"      {lab}: n={len(d)}")
        log(f"        vs top-clone fraction  rho {rho:+.3f}  p {pv:.4g}")
        log(f"        vs Simpson             rho {rho_s:+.3f}  p {pv_s:.4g}")
        rows.append({"test": f"E_dose_{xcol}", "puck": p, "radius": R,
                     "chain": chain, "n": len(d), "rho": float(rho),
                     "p": float(pv), "rho_simpson": float(rho_s),
                     "p_simpson": float(pv_s)})

    ns, nw = int(d.is_strong.sum()), int((~d.is_strong).sum())
    if ns >= 3 and nw >= 3:
        try:
            _, pv = mannwhitneyu(d.loc[d.is_strong, "top_frac"],
                                 d.loc[~d.is_strong, "top_frac"],
                                 alternative="greater")
        except ValueError:
            pv = np.nan
        log(f"      strong binders (n={ns}) vs weak (n={nw}): "
            f"median {d.loc[d.is_strong,'top_frac'].median():.4f} vs "
            f"{d.loc[~d.is_strong,'top_frac'].median():.4f}, p {pv:.4g}")
        rows.append({"test": "E_strong_vs_weak", "puck": p, "radius": R,
                     "chain": chain, "n_strong": ns, "n_weak": nw, "p": float(pv)})
    else:
        log(f"      strong/weak split too uneven ({ns}/{nw}) to test")


def power_test_A(ct, kt, pairs_p, p, R, rows, n_sim=2000):
    """Power for the PAIRED design. The Step08c bound (shift >= 0.20 detectable)
    was computed for the PER-ANCHOR design and does NOT apply here: a paired test
    across all anchors is powered for a small shift present in most of them,
    which is the opposite alternative."""
    c = ct[~ct.below_floor].set_index("anchor")
    k = kt[~kt.below_floor].set_index("anchor")
    ok = [(a, b) for a, b in pairs_p if a in c.index and b in k.index]
    if len(ok) < 5:
        return
    base = np.array([k.loc[b, "top_frac"] for _, b in ok])
    sd = np.std(np.array([c.loc[a, "top_frac"] for a, _ in ok]) - base, ddof=1)
    n = len(ok)
    log(f"    POWER for test A (paired, n={n}, observed paired SD {sd:.4f})")
    log(f"      {'shift':>8} {'rel. to base':>13} {'power':>8}")
    detect = None
    for shift in (0.002, 0.005, 0.0075, 0.010, 0.015, 0.020, 0.030):
        hits = 0
        for _ in range(n_sim // 4):
            diffs = RNG.normal(shift, sd, n)
            try:
                _, pv = wilcoxon(diffs, alternative="greater")
                hits += pv < ALPHA
            except ValueError:
                pass
        pw = hits / (n_sim // 4)
        log(f"      {shift:>8.4f} {100*shift/max(1e-9, np.median(base)):>12.0f}% "
            f"{pw:>8.2f}")
        if detect is None and pw >= 0.80:
            detect = shift
        rows.append({"test": "A_power", "puck": p, "radius": R, "shift": shift,
                     "power": round(pw, 3), "n_pairs": n, "paired_sd": float(sd)})
    log(f"      smallest shift at 80% power: "
        + (f"{detect:.4f} ({100*detect/max(1e-9,np.median(base)):.0f}% of the "
           f"control median)" if detect else "> 0.030"))


def test_C(adata, tcr_p, cases, p, R, chain, rows):
    log(f"    TEST C  neoantigen potential field")
    idx, f, xy, pos = puck_view(adata, p)
    epi = np.array([i for i, b in enumerate(idx) if f.iloc[i]["annot"] == EPITHELIAL])
    cxy = np.array([xy[pos[a]] for a in cases if a in pos])
    if len(epi) < 200 or len(cxy) < 3:
        log("      too few beads or anchors"); return None

    tcr_beads = [b for b in tcr_p.bead.unique() if b in pos]
    if len(tcr_beads) < 50:
        log("      too few TCR beads"); return None
    umi = tcr_p.groupby("bead")["umi"].nunique()
    tb = np.array([pos[b] for b in tcr_beads])
    w = np.array([umi.get(b, 0) for b in tcr_beads], dtype=float)

    # local TCR UMI density at every epithelial bead, computed ONCE
    etree = cKDTree(xy[epi])
    dens = np.zeros(len(epi))
    for src, j in enumerate(tb):
        for d in etree.query_ball_point(xy[j], R):
            dens[d] += w[src]

    depth = f.iloc[epi]["depth"].values
    y = np.log1p(dens)
    centroid = xy[epi].mean(axis=0)
    E = xy[epi]

    log(f"      {len(epi):,} epithelial beads, {len(cxy)} anchors, "
        f"{len(tb):,} TCR beads, {N_ROT} rotations")
    for lam in KERNEL_LAMBDAS_REPORT:
        def potential(pts):
            d = np.linalg.norm(E[:, None, :] - pts[None, :, :], axis=2)
            return np.exp(-d / lam).sum(axis=1)

        pot = potential(cxy)
        X = np.vstack([np.ones(len(epi)), np.log1p(depth), pot]).T
        obs = float(np.linalg.lstsq(X, y, rcond=None)[0][2])

        # ROTATE THE ANCHORS, not the TCR layer: mathematically equivalent and
        # far cheaper, since anchors number tens and TCR beads number thousands.
        null = np.empty(N_ROT)
        rel = cxy - centroid
        for i in range(N_ROT):
            th = RNG.uniform(0, 2 * np.pi)
            cth, sth = np.cos(th), np.sin(th)
            rot = np.column_stack([rel[:, 0]*cth - rel[:, 1]*sth,
                                   rel[:, 0]*sth + rel[:, 1]*cth]) + centroid
            Xn = np.vstack([np.ones(len(epi)), np.log1p(depth), potential(rot)]).T
            null[i] = np.linalg.lstsq(Xn, y, rcond=None)[0][2]
        pv = (1 + int(np.sum(null >= obs))) / (N_ROT + 1)
        log(f"      lambda {lam:6.0f} um: beta {obs:+.5f} | null mean "
            f"{null.mean():+.5f} sd {null.std():.5f} | p = {pv:.4g}"
            + ("   <- PRE-SPECIFIED, combined" if lam == KERNEL_LAMBDA
               else "   (reported only)"))
        primary = (lam == KERNEL_LAMBDA)
        rows.append({"test": "C_potential_field" if primary
                             else "C_potential_field_reported",
                     "puck": p, "radius": R, "chain": chain, "lambda_um": lam,
                     "beta": obs, "null_mean": float(null.mean()),
                     "null_sd": float(null.std()), "p": pv,
                     "is_primary_lambda": primary})
    log(f"      Only lambda = {KERNEL_LAMBDA:.0f} um enters the combination. v1")
    log("      combined all three per puck with Stouffer at n=9, treating three")
    log("      views of the same data on the same puck as independent. They are")
    log("      strongly correlated, so that overstated the evidence.")
    log("      Rotation preserves the spatial autocorrelation of both layers")
    log("      while destroying their alignment. Rotated anchors may land off")
    log("      tissue, which makes the null slightly conservative.")


# ==============================================================================
def main():
    global R_CUR
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--chain", default="TRB", help="TRB (default) or all")
    ap.add_argument("--radii", type=float, nargs="+", default=RADII)
    ap.add_argument("--skip-C", action="store_true")
    ap.add_argument("--rarefy-reps", type=int, default=N_RAREFY)
    ap.add_argument("--force-fallback", action="store_true",
                    help="use hardcoded floors if Step08c output for this chain "
                         "is absent; results are NOT trustworthy")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    needed = [ADATA, NEO, PAIRS, TESTABLE, POOL] + \
             [os.path.join(TCR_DIR, TCR_CSV[p]) for p in PUCKS]
    missing = [f for f in needed if not os.path.exists(f)]
    log(f"path check: {len(missing)} missing")
    for f in missing:
        log(f"   MISSING {f}")
        if f in (PAIRS, TESTABLE, POOL):
            log("     run Step08c first; it writes this")
    if args.validate_only:
        sys.exit(1 if missing else 0)
    if missing:
        sys.exit("cannot proceed")

    log(f"\nloading {ADATA}")
    adata = ad.read_h5ad(ADATA)
    log(f"  {adata.n_obs:,} beads x {adata.n_vars:,} genes")

    neo = pd.read_csv(NEO, sep="\t")
    pairs = pd.read_csv(PAIRS, sep="\t"); pairs["puck"] = pairs["puck"].astype(str)
    testable = pd.read_csv(TESTABLE, sep="\t"); testable["puck"] = testable["puck"].astype(str)
    pool_df = pd.read_csv(POOL, sep="\t"); pool_df["puck"] = pool_df["puck"].astype(str)

    FLOORS = load_floors(args.chain)
    if FLOORS is None:
        if args.force_fallback:
            FLOORS = FLOORS_FALLBACK
            log(f"  *** --force-fallback: using {FLOORS}, which were chosen on a")
            log(f"      DIFFERENT chain set. Results are not trustworthy. ***")
        else:
            sys.exit("cannot proceed without matched floors (--force-fallback to override)")
    missing_f = [p for p in PUCKS if p not in FLOORS]
    if missing_f:
        sys.exit(f"no floor for puck(s) {missing_f}")

    log(f"\nchain {args.chain} | radii {args.radii} (primary {args.radii[0]:.0f} um) | "
        f"floors {FLOORS} | k {K_VALUES} | rarefaction {args.rarefy_reps} reps")
    log(f"enrichment restricted to clones with >= {MIN_LOCAL_UMI} local and "
        f">= {MIN_PUCK_UMI} puck-wide UMIs")
    log("TEST B statistic: rarefied TOP-CLONE FRACTION (the units the")
    log("  pre-registered bound is written in). The v1 enrichment statistic was")
    log("  degenerate: saturated by singleton clones, one distinct value, p = 1.")

    rows = []
    preflight(adata, neo, pairs, testable, rows)

    info = neo.drop_duplicates("CB").set_index("CB")
    hits, a_store = [], {}
    log("\n  REPORTING ROLES, pre-specified: B is PRIMARY. A was specified as the")
    log("  control check. If B is null and A is not, A is a SECONDARY finding and")
    log("  must be reported as such.")

    for R in args.radii:
        R_CUR = R
        banner(f"RADIUS {R:.0f} um   (chain {args.chain})")
        for p in PUCKS:
            idx, f, xy, pos = puck_view(adata, p)
            tree = cKDTree(xy)
            tp = load_tcr(p, args.chain)
            tp = tp[tp.bead.isin(set(idx))]
            if not len(tp):
                log(f"\n  puck {p}: no TCR rows"); continue
            tot = tp.groupby("cloneId")["umi"].nunique()
            puck_freq = (tot / tot.sum()).to_dict()
            puck_umi = tot.to_dict()
            bead_clones = build_index(tp, set(idx))
            floor = FLOORS[p]

            pairs_p = list(zip(pairs.loc[pairs.puck == p, "case"],
                               pairs.loc[pairs.puck == p, "control"]))
            tcases = list(testable.loc[testable.puck == p, "CB"])
            pbeads = list(pool_df.loc[pool_df.puck == p, "CB"])

            log(f"\n  --- puck {p}: floor {floor} UMIs, {len(pairs_p)} pairs, "
                f"{len(tcases)} testable, pool {len(pbeads):,} ---")

            ct_A = score_anchors([a for a, _ in pairs_p], tree, xy, idx, pos,
                                 bead_clones, floor, puck_freq, puck_umi, K_VALUES)
            kt_A = score_anchors([b for _, b in pairs_p], tree, xy, idx, pos,
                                 bead_clones, floor, puck_freq, puck_umi, K_VALUES)
            rA = test_A(ct_A, kt_A, pairs_p, rows, p, R, args.chain)
            if rA:
                a_store[(p, R)] = rA
            power_test_A(ct_A, kt_A, pairs_p, p, R, rows)

            ct_B = score_anchors(tcases, tree, xy, idx, pos, bead_clones,
                                 floor, puck_freq, puck_umi, K_VALUES)
            pt_B = score_anchors(pbeads, tree, xy, idx, pos, bead_clones,
                                 floor, puck_freq, puck_umi, K_VALUES)
            for k in K_VALUES:
                col = f"n_qualifying_k{k}"
                if col in ct_B.columns:
                    sub = ct_B[~ct_B.below_floor]
                    if len(sub):
                        log(f"      clones qualifying for enrichment at k>={k}: "
                            f"median {sub[col].median():.1f} per case anchor, "
                            f"{100*(sub[col] < 1).mean():.0f}% of anchors have none")

            rB = test_B(ct_B, pt_B, rows, p, R, args.chain, "top_frac",
                        "B_topclone")
            for k in K_VALUES:
                rD = test_B(ct_B, pt_B, rows, p, R, args.chain, f"max_enr_k{k}",
                            f"D_in_situ_k{k}", k=k)
                if rD is None:
                    continue
                for _, h in rD["res"][rD["res"].hit].iterrows():
                    labels, cb = neighbourhood(tree, xy, idx, pos, bead_clones,
                                               h.anchor, R)
                    wc = winning_clone(labels, puck_freq, puck_umi, cb, k=k) or {}
                    ii = info.loc[h.anchor] if h.anchor in info.index else None
                    hits.append({"puck": p, "radius": R, "k": k,
                                 "anchor": h.anchor, "p": h.p, "p_bh": h.p_bh,
                                 "statistic": h.stat,
                                 "gene": ii["gene"] if ii is not None else None,
                                 "hgvs_p": ii["hgvs_p"] if ii is not None else None,
                                 "peptide": ii["mut_peptide"] if ii is not None else None,
                                 "allele": ii["best_allele"] if ii is not None else None,
                                 **wc})

            anchor_independence(tcases, xy, pos, R, p, rows)
            test_E_dose_response(ct_B, info, xy, pos, tcases, R, p,
                                 args.chain, rows)
            if not args.skip_C:
                test_C(adata, tp, tcases, p, R, args.chain, rows)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, f"test_results_{args.chain}.tsv"), sep="\t", index=False)

    banner("COMBINED ACROSS PUCKS (Stouffer)")
    log("  Three patients, private repertoires: evidence is combined, never data.")
    SKIP_COMBINE = {"C_potential_field_reported"}
    for test in sorted(t for t in df.get("test", pd.Series(dtype=str)).dropna().unique()):
        if test in SKIP_COMBINE:
            log(f"  {test:<22} not combined (correlated lambdas on the same data)")
            continue
        for R in args.radii:
            s = df[(df.test == test) & (df.radius == R)]
            pv = []
            if "p" in s.columns:
                pv += s["p"].dropna().tolist()
            if test == "A_clonality" and "p_wilcoxon" in s.columns:
                pv = s["p_wilcoxon"].dropna().tolist()
            if len(pv) >= 2:
                comb = stouffer(pv)
                log(f"  {test:<22} R={R:>5.0f}  n={len(pv)}  combined p = {comb:.4g}"
                    + ("   <-- significant" if comb < ALPHA else ""))
                rows.append({"test": f"{test}_combined", "radius": R,
                             "chain": args.chain, "p": comb, "n_pucks": len(pv)})

    if hits:
        hdf = pd.DataFrame(hits).sort_values("p_bh")
        hdf.to_csv(os.path.join(OUT, f"hits_{args.chain}.tsv"), sep="\t", index=False)
        banner(f"HITS  ({len(hdf)} anchor-clone pairs at BH p < {ALPHA})")
        for _, h in hdf.iterrows():
            log(f"  puck {h.puck} R={h.radius:.0f} k={h.k}  BH p={h.p_bh:.4g}  "
                f"enrichment {h.get('enrichment', float('nan')):.2f}")
            log(f"    {h.gene} {h.hgvs_p} | {h.peptide} | {h.allele}")
            log(f"    clone {h.get('clone')} on {h.get('clone_beads')} beads, "
                f"{h.get('clone_umis')} UMIs")
            log(f"    anchor {h.anchor}")
        log("\n  CANDIDATES, not confirmed reactive pairs. Each needs its CDR3")
        log("  inspected and, ideally, orthogonal validation.")
    else:
        banner("NO HITS")
        log("  No anchor reached BH-adjusted p < 0.05.")
        log("  Report as a BOUNDED negative. Step08c power: a shift of 0.20 in")
        log("  top-clone fraction is detectable at 0.93-1.00 even if only ONE")
        log("  anchor is affected; a shift of 0.10 at only 0.14-0.17. Effects")
        log("  below 0.10 are outside what this design can see. This is NOT")
        log("  evidence that neoantigens do not recruit T cells.")

    # figures
    if a_store:
        fig, axes = plt.subplots(1, len(PUCKS), figsize=(30, 9), sharey=True)
        for i, p in enumerate(PUCKS):
            ax = axes[i]
            r = a_store.get((p, args.radii[0]))
            if r is None:
                style(ax, title=f"Puck {p} (n/a)"); continue
            for ca, kb in zip(r["ca"], r["kb"]):
                ax.plot([0, 1], [kb, ca], color=C_GRAY, alpha=0.5, linewidth=1.5)
            ax.scatter(np.zeros(len(r["kb"])), r["kb"], s=90, color=C_MUSTARD,
                       zorder=3, label="mutation only")
            ax.scatter(np.ones(len(r["ca"])), r["ca"], s=90, color=C_CORAL,
                       zorder=3, label="neoantigen")
            ax.set_xticks([0, 1]); ax.set_xticklabels(["control", "case"])
            style(ax, title=f"Puck {p}  p={r['p']:.3g}",
                  ylabel="Rarefied top-clone fraction" if i == 0 else None)
            if i == 0:
                ax.legend(fontsize=FS_ANNOT, frameon=False)
        fig.tight_layout(); savefig(fig, f"Fig_08d_A_clonality_{args.chain}")

    pd.DataFrame(rows).to_csv(
        os.path.join(OUT, f"test_results_{args.chain}.tsv"), sep="\t", index=False)
    log(f"\noutputs in {OUT}")


if __name__ == "__main__":
    main()
