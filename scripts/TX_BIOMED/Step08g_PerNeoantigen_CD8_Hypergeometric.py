#!/usr/bin/env python3
"""
Step08g_PerNeoantigen_CD8_Hypergeometric.py
================================================================================
PER-ANCHOR test, run on EVERY predicted neoantigen: is any TCR clone locally
over-represented on CD8-dominated beads around this specific neoantigen?

WHY THIS IS RUNNABLE WHEN DESIGN 1 WAS NOT
    Design 1 died on the rarefaction floor. The top-clone fraction is a MAXIMUM,
    a maximum computed from fewer draws is biased upward, so it must be rarefied
    to a common floor, and the CD8-restricted pool could not reach that floor
    (retention 0/30, 0/18, 2/55 at 100 um).

    The hypergeometric needs no rarefaction. It CONDITIONS ON n ANALYTICALLY:
    given n local draws from a known urn, is clone c present more often than
    chance? Valid at n = 6 as at n = 600, just less powerful. It is the one
    statistic in the family that survives the CD8 restriction.

    Three consequences:
      - no floor, so ALL 110 neoantigen anchors are tested, not the 101-103 in
        testable_cases.tsv, whose membership was defined by a floor that no
        longer applies
      - conditioning on n absorbs most of the Step08e BEAT 4 depth confound: a
        deeper neighbourhood yields a larger n, and n is conditioned on rather
        than compared
      - it is natively per-anchor, which is the question being asked

WHAT THIS CAN AND CANNOT CLAIM
    The MiXCR table carries NO subtype label, and TRB is expressed by CD8, CD4
    and regulatory T cells alike. There is no measurement anywhere in this
    dataset of "a CD8 cell's TCR". The closest honest construction is: TRB
    signal on beads whose transcriptome is CD8-dominated.

    SO THE CLAIM IS "TCR CLONALITY RESTRICTED TO CD8-DOMINATED BEADS",
    NEVER "CD8 CLONAL EXPANSION". Step08e BEAT 3 showed the c2l split only
    weakly separates CD8A/B-positive from negative T beads once the 2x depth
    difference between those groups is accounted for.

REGIONAL VS LOCAL, AND WHY PART C's FAILURE DOES NOT TRANSFER
    Step08f PART C failed because c2l CD8 ABUNDANCE has regional structure that
    spatially unconstrained matching does not control. TCR clonality has no such
    structure to control: Step08a measured radius of gyration at 1.03-1.05x the
    random-labelling null at p = 1.0000, and neighbours per clone at 30 um at
    0.018-0.019. A puck-wide background is therefore defensible here in a way it
    was not there. ARM B runs the regional version anyway as the sensitivity; if
    the two disagree, that contradicts Step08a and needs chasing before any hit
    is believed.

THE UNIT: SPATIAL FOCI, AND WHAT IT COSTS
    UMIs are not error-corrected upstream (3.9-5.0% of rows are within one
    substitution of another on the SAME bead and clone), and one cell's
    transcripts diffuse onto neighbouring beads. Counting UMIs therefore counts
    one molecule twice and one cell many times.

    A FOCUS is a connected component of one clone's beads, linked whenever two
    beads sit within the diffusion radius Step08b MEASURED. One cell is one
    focus. Two cells 80 um apart are two foci.

    THE COST, STATED PLAINLY: clonal expansion by local proliferation puts
    daughter cells ADJACENT to each other, and adjacent daughter cells merge
    into ONE focus. At 10 um beads with ~10 um cells, diffusion and adjacent
    proliferation are not geometrically separable in this data. The focus unit
    is therefore CONSERVATIVE against the single most biologically expected
    pattern. Beads and UMIs run as sensitivity arms so the cost is visible as a
    number rather than asserted.

    The molecular alternative (flagging foci whose UMI count exceeds a one-cell
    reference) is deliberately NOT built: the only available reference is
    single-bead singleton clones, and Step08b BEAT 10 found those are enriched
    for error split off real clones, so the reference is biased low and the flag
    would over-call. Dropped at Jake's instruction rather than shipped biased.

THE URN MUST MATCH THE UNIT
    A naive implementation scores LOCAL FOCI against a PUCK-WIDE UMI total. That
    is invalid. Both numerator and denominator are foci here, which requires the
    focus decomposition to run PUCK-WIDE, not just inside each neighbourhood.
    ARM B's urn is every focus within REGIONAL_RADIUS of the anchor, with the
    local disc NESTED INSIDE it, which is a proper hypergeometric setup; an urn
    that excluded the local region would not be.

WITHIN-ANCHOR MULTIPLICITY: BONFERRONI IS NOT ENOUGH, AND THIS WAS MEASURED
    Taking the minimum p over the clones present locally makes a MAXIMUM
    statistic. The first build corrected it analytically, Bonferroni over the
    number of distinct local clones. A type I check on a RANDOM FIELD with no
    clonal structure whatsoever killed that:

        3,000 clones / 4,000 foci   raw 0.182   Bonferroni 0.182
        9,000 clones / 18,000 foci  raw 0.723   Bonferroni 0.492
                                                 (nominal 0.05)

    Bonferroni fails because its reference set is DATA-DEPENDENT. Only clones
    already present locally with k >= 2 are tested, so the correction is applied
    over clones that were SELECTED FOR HAVING SHOWN THE SIGNAL. The valid
    reference is every clone that could have landed in that neighbourhood, which
    has no closed form here.

    So the operative per-anchor p is EMPIRICAL, calibrated against the identical
    statistic computed at ROTATED anchor positions. Rotation keeps the anchor
    constellation rigid, so inter-anchor overlap survives into the null. p_raw
    and p_bonf are still reported as columns; THE HIT LIST IS BUILT ON p_emp,
    and the spike-in uses the same calibrated threshold so power and type I sit
    on one scale.

    This is why the regional arm B, which cannot be rotated cheaply, reports
    uncalibrated counts that must NOT be read the way the primary arm's are.

STRINGENCY LADDER, BH LAST
    Per Jake: start broad, add stringency, keep every rung visible.
      rung 1  raw per-anchor p and hit count          (type I 0.18-0.72, see above)
      rung 2  within-anchor Bonferroni                  (type I 0.18-0.49, fails)
      rung 3  ROTATION-CALIBRATED per-anchor p_emp      <== the operative rung
              Anchors are not independent: 31/17/53 anchors form 28/13/36
              effective groups once overlapping neighbourhoods are counted, so
              a uniform null would be wrong
      rung 4  p-value distribution of p_emp vs uniform, the only rung with power
              against a DIFFUSE subset (five anchors at p ~ 0.05 that BH flags
              individually at none and test A dilutes against 48 nulls)
      rung 5  BH on p_emp, reported as a final column

    Rungs 1 and 2 are kept visible because Jake asked to start broad, but they
    are DIAGNOSTICS OF THE MECHANISM, not candidate lists. The gap between the
    raw count and the calibrated count IS the k=2 diffusion artefact, measured.

WHAT TO EXPECT, WRITTEN DOWN BEFORE RUNNING
    At R = 100 um restricted to the top CD8 tertile, local focus counts will be
    thin. The test WILL generate raw hits, most at k = 2, and the type I numbers
    above say the great majority of those are the mechanism rather than biology.
    Expect the calibrated count to collapse toward 5% of testable anchors. If it
    does not, that excess is the result.

    COUNTERINTUITIVE: restricting the urn to CD8 beads SHRINKS the denominator,
    which RAISES every clone's background frequency, which REDUCES the evidence
    any local count carries. PART 3's detectability calculation accounts for it.
    Intuition does not.

PRE-SPECIFICATION RECORD
    This is the THIRD parallel design on the same 110 anchors, and the most
    permissive yet. 08d was the original; 08f was post hoc relative to 08d;
    08g is post hoc relative to both. Any hit here is A CANDIDATE FOR FOLLOW-UP,
    NOT A FINDING, and the output says so in its own hit table.

STANDING LIMIT
    The anchors are still predicted against a fixed ten-allele panel rather than
    patient HLA, so a share are peptides the patient cannot present. HLA typing
    from the uBAMs is upstream of any immunogenicity or pMHC-TCR work.

USAGE
    python Step08g_PerNeoantigen_CD8_Hypergeometric.py
    python Step08g_PerNeoantigen_CD8_Hypergeometric.py --diffusion-radius 30
    python Step08g_PerNeoantigen_CD8_Hypergeometric.py --n-perm 200 --n-spike 200
    python Step08g_PerNeoantigen_CD8_Hypergeometric.py --validate-only

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
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import hypergeom, kstest

try:
    import anndata as ad
except ImportError:
    sys.exit("anndata not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTPUTS = os.path.join(PROOT, "data/outputs")
SALVAGE = os.path.join(OUTPUTS, "10_tcr_salvage")
FEAS = os.path.join(OUTPUTS, "11_neoantigen_tcr")
OUT = os.path.join(OUTPUTS, "15_per_neoantigen_cd8")
FIGDIR = os.path.join(OUT, "figures")

ADATA = os.path.join(OUTPUTS, "04_annotation/all_pucks_annotated_unified.h5ad")
NEO = os.path.join(OUTPUTS, "07_neoantigen/neoantigens_per_bead.tsv")
TCR_DIR = os.path.join(PROOT, "data/inputs/tcr/processed")

PUCKS = ["29", "37", "40"]
TCR_CSV = {p: f"B59_{p}_hTCR_tcr.csv" for p in PUCKS}

C_CD8 = "c2l_CD8-positive, alpha-beta T cell"
C_CD4 = "c2l_CD4-positive, alpha-beta T cell"
C_TREG = "c2l_regulatory T cell"

# Jake moved the primary from 50 to 100 um specifically to accommodate the
# measured diffusion: at 50 um the neighbourhood is barely wider than the
# linking kernel and the focus count has no dynamic range.
PRIMARY_RADIUS = 100.0
RADII = [50.0, 75.0, 100.0, 150.0, 200.0]
REGIONAL_RADIUS = 600.0   # ARM B urn. Local disc is NESTED inside it.
DEFAULT_CHAIN = "TRB"

CD8_TERTILE = 2.0 / 3.0   # top tertile of depth-residualized CD8 share
ALPHA = 0.05
N_PERM = 200
N_SPIKE = 200
SPIKE_K = [2, 3, 5, 8, 12]          # LOCAL foci of the synthetic clone
# PUCK-WIDE focus count of the synthetic clone. K = k ("local") is the clone
# that exists ONLY at this anchor, which is the most detectable case and the
# only one the first build measured. A real expanded tumour-infiltrating clone
# is present THROUGHOUT the tumour, so K is large and the same local count
# carries far less evidence. Sweeping K is what turns the power statement from
# optimistic to honest.
SPIKE_K_PUCKWIDE = [None, 10, 25, 50, 100]   # None means K = k
MIN_LOCAL_FOCI = 2                  # below this an anchor cannot produce a hit

EPS = 1e-12

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
C_CORAL = "#ed6a5a"
C_MUSTARD = "#F6D155"
C_GRAY = "#d3d3d3"
C_DARK = "#4a4a4a"
C_BLUE = "#7eb0d5"
C_GREEN = "#b2e061"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

RNG = np.random.default_rng(20261028)


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


def bh(pvals):
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


def resolve_diffusion_radius(override):
    """Read the diffusion radius from Step08b. NEVER hardcoded.

    Same discipline as the rarefaction floors, for the same reason: a kernel
    measured on one data definition and applied to another is the v2 bug in a
    new costume. If the schema cannot be resolved the script STOPS and asks for
    an explicit value rather than guessing one.
    """
    if override is not None:
        log(f"  diffusion radius: {override:.1f} um (EXPLICIT --diffusion-radius)")
        log("    Confirm this against Step08b before trusting any hit.")
        return float(override)
    path = os.path.join(SALVAGE, "salvage_summary.tsv")
    if not os.path.exists(path):
        sys.exit(f"FATAL: {path} not found and no --diffusion-radius given. "
                 f"Run Step08b, or pass the measured value explicitly.")
    d = pd.read_csv(path, sep="\t")
    key = "beat" if "beat" in d.columns else d.columns[0]
    pat = d[key].astype(str).str.lower()
    hits = d[pat.str.contains("diffus|kernel|spread|halo", regex=True, na=False)]
    numcols = [c for c in d.columns
               if any(t in c.lower() for t in ("radius", "um", "dist", "value"))]
    for _, r in hits.iterrows():
        for c in numcols:
            try:
                v = float(r[c])
            except (TypeError, ValueError):
                continue
            if np.isfinite(v) and 5.0 <= v <= 200.0:
                log(f"  diffusion radius: {v:.1f} um "
                    f"(from salvage_summary.tsv, beat='{r[key]}', col='{c}')")
                return v
    log("  Could not resolve a diffusion radius from salvage_summary.tsv.")
    log(f"  Rows whose key matched diffusion/kernel/spread: {len(hits)}")
    if len(hits):
        log(hits.to_string(max_rows=20))
    sys.exit("FATAL: pass the measured Step08b value with --diffusion-radius. "
             "Not guessed, not defaulted.")


# ==============================================================================
# loading
# ==============================================================================
def validate_paths():
    ok = True
    need = [("adata", ADATA), ("neoantigens", NEO)]
    for p in PUCKS:
        need.append((f"tcr {p}", os.path.join(TCR_DIR, TCR_CSV[p])))
    for name, path in need:
        hit = os.path.exists(path)
        log(f"  {'OK ' if hit else '***'} {name:<16} {path}")
        ok &= hit
    for name, path in (("salvage summary", os.path.join(SALVAGE, "salvage_summary.tsv")),):
        log(f"  {'OK ' if os.path.exists(path) else '(-)'} {name:<16} {path}")
    return ok


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
    log("  Clones grouped on cloneId ONLY. vGene, cGene and aaSeqCDR3 all vary")
    log("  within one cloneId from low-depth error and are never grouping keys.")
    return t.reset_index(drop=True)


def puck_frame(adata, p):
    m = adata.obs["puck_id"].astype(str).values == f"Puck_211214_{p}"
    idx = adata.obs_names[m]
    o = adata.obs.loc[idx]
    f = pd.DataFrame(index=idx)
    f["depth"] = o["total_counts"].values.astype(float)
    f["log_depth"] = np.log(np.clip(f["depth"].values, 1, None))
    f["cd8"] = o[C_CD8].values.astype(float)
    f["cd4"] = o[C_CD4].values.astype(float)
    f["treg"] = o[C_TREG].values.astype(float)
    f["t_mass"] = f[["cd8", "cd4", "treg"]].sum(axis=1)
    f["share"] = np.where(f["t_mass"].values > EPS,
                          f["cd8"].values / np.maximum(f["t_mass"].values, EPS),
                          0.0)
    xy = o[["x_coord", "y_coord"]].values.astype(float)
    pos = {b: i for i, b in enumerate(idx)}
    return f, xy, pos


def cd8_tertile_mask(f):
    """Top tertile of CD8 share, RESIDUALIZED ON LOG DEPTH within puck.

    Step08e BEAT 4: the raw share climbs ~3x across depth deciles (Spearman
    +0.95 / +0.81 / +0.93), so a raw-share tertile would be a depth tertile
    wearing a CD8 label. Step08f PART D's puck 40 split still came back at
    SMD +0.426 after this same residualization, so the correction is partial,
    not complete. The balance is reported per puck below, not assumed.
    """
    X = np.column_stack([np.ones(len(f)), f["log_depth"].values])
    beta, *_ = np.linalg.lstsq(X, f["share"].values, rcond=None)
    resid = f["share"].values - X @ beta
    cut = float(np.quantile(resid, CD8_TERTILE))
    return resid > cut, resid, cut


# ==============================================================================
# focus decomposition
# ==============================================================================
def decompose_clone(xy_sub, link_radius):
    """Connected components of one clone's beads, linked within link_radius."""
    n = len(xy_sub)
    if n == 1:
        return np.zeros(1, int), 1
    tree = cKDTree(xy_sub)
    pairs = np.array(sorted(tree.query_pairs(link_radius)), dtype=int)
    if not len(pairs):
        return np.arange(n), n
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                   shape=(n, n))
    k, lab = connected_components(g, directed=False)
    return lab, k


def build_foci(tcr_puck, bead_xy, link_radius):
    """Puck-wide focus table for one arm.

    PUCK-WIDE is the point: the urn and the local count must be the same unit,
    so the decomposition cannot be done only inside neighbourhoods.
    """
    recs = []
    for clone, sub in tcr_puck.groupby("cloneId", sort=False):
        beads = sub["bead"].unique()
        pts = np.array([bead_xy[b] for b in beads])
        lab, k = decompose_clone(pts, link_radius)
        umi_by_bead = sub.groupby("bead")["umi"].nunique()
        for fid in range(k):
            m = lab == fid
            p = pts[m]
            bs = beads[m]
            if len(p) > 1:
                d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1)
                dmax = float(d.max())
            else:
                dmax = 0.0
            recs.append({"cloneId": clone, "focus": fid,
                         "x": float(p[:, 0].mean()), "y": float(p[:, 1].mean()),
                         "n_beads": int(m.sum()),
                         "n_umis": int(umi_by_bead.reindex(bs).fillna(0).sum()),
                         "max_pair_um": dmax})
    return pd.DataFrame(recs)


# ==============================================================================
# the per-anchor test
# ==============================================================================
def anchor_test(tree, fx, clone_arr, counts_total, N_total, axy, R,
                extra_clone=None, extra_k=0, extra_K=0):
    """Hypergeometric over FOCI for every clone present within R of an anchor.

    Returns the minimum raw p, its clone, the within-anchor Bonferroni-corrected
    p, n_local, and the number of distinct local clones tested.

    The minimum over clones is a MAXIMUM statistic, which the correction log
    says cannot be tested by label permutation. It is corrected ANALYTICALLY
    here (Bonferroni on the number of distinct local clones). The ACROSS-anchor
    structure is handled by rotation in PART 4, not here.
    """
    near = tree.query_ball_point(axy, R)
    n = len(near) + extra_k
    if n < MIN_LOCAL_FOCI:
        return None
    best_p, best_clone, best_k, best_K = 1.0, None, 0, 0
    tested = 0
    if near:
        loc = clone_arr[np.asarray(near)]
        uc, cc = np.unique(loc, return_counts=True)
        tested = len(uc)
        for c, k in zip(uc, cc):
            if k < 2:
                continue          # k = 1 is uninformative at any urn size
            K = counts_total.get(c, 0)
            if K < k:
                K = k
            pv = float(hypergeom.sf(k - 1, N_total, K, n))
            if pv < best_p:
                best_p, best_clone, best_k, best_K = pv, c, int(k), int(K)
    if extra_clone is not None and extra_k >= 2:
        tested += 1
        pv = float(hypergeom.sf(extra_k - 1, N_total + extra_K,
                                extra_K, n))
        if pv < best_p:
            best_p, best_clone, best_k, best_K = pv, extra_clone, extra_k, extra_K
    if best_clone is None:
        return {"n_local": int(n), "n_clones": int(tested), "p_raw": np.nan,
                "p_bonf": np.nan, "clone": None, "k": 0, "K": 0}
    return {"n_local": int(n), "n_clones": int(tested), "p_raw": best_p,
            "p_bonf": float(min(1.0, best_p * max(1, tested))),
            "clone": best_clone, "k": best_k, "K": best_K}


def min_detectable_k(N_total, counts_total, n, n_clones, alpha=ALPHA):
    """Smallest k that could ever clear the threshold at this anchor.

    The most favourable REALIZABLE clone is one whose entire puck-wide presence
    is local, K = k, which is exactly the expansion hypothesis and exactly what
    the spike-in injects. Scanning K = 1 instead would be a bug: a clone with
    K = 1 cannot appear twice locally, so hypergeom returns p = 0 for an
    IMPOSSIBLE event and every anchor would report a min detectable k of 2.

    If the answer exceeds n, this anchor is structurally incapable of producing
    a hit and is reported as UNTESTED rather than null.
    """
    if n < MIN_LOCAL_FOCI:
        return None
    thresh = alpha / max(1, n_clones)
    for k in range(2, int(n) + 1):
        if float(hypergeom.sf(k - 1, N_total, k, n)) < thresh:
            return k
    return None


# ==============================================================================
def main():
    ap = argparse.ArgumentParser(
        description="Per-neoantigen CD8-restricted hypergeometric on foci.")
    ap.add_argument("--chain", default=DEFAULT_CHAIN)
    ap.add_argument("--radius", type=float, default=PRIMARY_RADIUS)
    ap.add_argument("--radii", type=float, nargs="+", default=RADII)
    ap.add_argument("--diffusion-radius", type=float, default=None,
                    help="Step08b measured value. Resolved from "
                         "salvage_summary.tsv when omitted; never defaulted.")
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    ap.add_argument("--n-spike", type=int, default=N_SPIKE)
    ap.add_argument("--skip-spike", action="store_true")
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()
    R0 = args.radius

    banner("Step08g_PerNeoantigen_CD8_Hypergeometric")
    log(f"  started  {datetime.now():%Y-%m-%d %H:%M:%S}")
    log(f"  chain    {args.chain}")
    log(f"  primary  {R0:.0f} um   sweep {args.radii}")
    log()
    if not validate_paths():
        sys.exit("FATAL: inputs missing.")
    link = resolve_diffusion_radius(args.diffusion_radius)
    if args.validate_only:
        log("\n--validate-only: paths and diffusion radius resolve. Stopping.")
        return

    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)
    rows, hits, spike_rows = [], [], []

    log()
    log("loading annotated object (backed, obs only)")
    adata = ad.read_h5ad(ADATA, backed="r")
    log(f"  {adata.n_obs:,} beads x {adata.n_vars:,} genes")
    tcr = load_tcr(args.chain)
    neo = pd.read_csv(NEO, sep="\t")
    ncol = "CB" if "CB" in neo.columns else neo.columns[0]
    log(f"  neoantigen anchors in file: {neo[ncol].nunique()} unique beads")
    log("  ALL of them are tested. testable_cases.tsv is NOT used: its")
    log("  membership was defined by a rarefaction floor this test does not have.")

    ARMS = [("A_cd8_puckwide", "CD8 tertile, puck-wide urn", True, False),
            ("B_cd8_regional", "CD8 tertile, regional urn", True, True),
            ("C_all_puckwide", "unrestricted TRB, puck-wide urn", False, False)]

    for p in PUCKS:
        f, xy, pos = puck_frame(adata, p)
        bead_xy = {b: xy[i] for b, i in pos.items()}
        hi, resid, cut = cd8_tertile_mask(f)
        d_hi = float(np.mean(f.loc[hi, "log_depth"]) - np.mean(f.loc[~hi, "log_depth"]))
        sd = float(np.sqrt((np.var(f.loc[hi, "log_depth"], ddof=1) +
                            np.var(f.loc[~hi, "log_depth"], ddof=1)) / 2))
        smd_depth = d_hi / sd if sd > 0 else np.nan
        cd8_beads = set(f.index[hi])
        tp = tcr[tcr.puck == p]
        tp = tp[tp.bead.isin(set(f.index))].copy()
        anchors = [b for b in neo[ncol].astype(str).unique()
                   if b.endswith(f"Puck_211214_{p}") and b in pos]

        banner(f"PUCK {p}")
        log(f"  beads {len(f):,} | CD8 top tertile {int(hi.sum()):,} "
            f"(residual cut {cut:+.4f})")
        log(f"  log-depth SMD, CD8-high vs rest: {smd_depth:+.3f} "
            f"{'(balanced)' if abs(smd_depth) < 0.10 else '*** IMBALANCED, FLAGGED ***'}")
        log(f"    Residualization is PARTIAL, not complete. Step08f PART D saw")
        log(f"    +0.426 in puck 40 after the same correction.")
        log(f"  neoantigen anchors on this puck: {len(anchors)}")
        rows.append({"part": "setup", "puck": p, "n_beads": len(f),
                     "n_cd8_tertile": int(hi.sum()), "resid_cut": cut,
                     "logdepth_smd": smd_depth, "n_anchors": len(anchors),
                     "link_radius_um": link})

        bead_tree = cKDTree(xy)   # for the local depth / CD8 diagnostic

        for arm, arm_lab, restrict, regional in ARMS:
            sub = tp[tp.bead.isin(cd8_beads)] if restrict else tp
            if sub.empty:
                log(f"\n  --- {arm}: no rows, skipping ---")
                continue
            foci = build_foci(sub, bead_xy, link)
            fx = foci[["x", "y"]].values
            clone_arr = foci["cloneId"].values
            tree = cKDTree(fx)
            N_all = len(foci)
            counts_all = foci.groupby("cloneId").size().to_dict()
            log(f"\n  --- {arm}  ({arm_lab}) ---")
            log(f"    UMIs {len(sub):,} -> beads {sub.bead.nunique():,} -> "
                f"FOCI {N_all:,} (link {link:.0f} um)")
            log(f"    urn shrink UMI->foci: {len(sub)/max(1,N_all):.2f}x. A "
                f"smaller urn RAISES every clone's background frequency and")
            log(f"    REDUCES the evidence any local count carries.")

            for R in args.radii:
                per = []
                for a in anchors:
                    axy = xy[pos[a]]
                    if regional:
                        reg = tree.query_ball_point(axy, REGIONAL_RADIUS)
                        if len(reg) < MIN_LOCAL_FOCI:
                            continue
                        reg = np.asarray(reg)
                        sub_fx = fx[reg]
                        sub_cl = clone_arr[reg]
                        st = cKDTree(sub_fx)
                        cnt = pd.Series(sub_cl).value_counts().to_dict()
                        r = anchor_test(st, sub_fx, sub_cl, cnt, len(reg), axy, R)
                    else:
                        r = anchor_test(tree, fx, clone_arr, counts_all,
                                        N_all, axy, R)
                    # local bead context, for the depth diagnostic. Computed
                    # over BEADS, not foci: the question is whether the hits
                    # are simply the deepest neighbourhoods, which would mean
                    # conditioning on n did not do its job.
                    nb = bead_tree.query_ball_point(axy, R)
                    loc_depth = float(np.mean(f["depth"].values[nb])) if nb else np.nan
                    loc_cd8 = float(np.mean(f["cd8"].values[nb])) if nb else np.nan
                    if r is None:
                        per.append({"anchor": a, "n_local": 0, "n_clones": 0,
                                    "p_raw": np.nan, "p_bonf": np.nan,
                                    "clone": None, "k": 0, "K": 0,
                                    "testable": False, "min_detectable_k": None,
                                    "local_depth": loc_depth,
                                    "local_cd8": loc_cd8,
                                    "n_local_beads": len(nb)})
                        continue
                    mk = min_detectable_k(N_all, counts_all, r["n_local"],
                                          max(1, r["n_clones"]))
                    r["testable"] = mk is not None and mk <= r["n_local"]
                    r["min_detectable_k"] = mk
                    r["anchor"] = a
                    r["local_depth"] = loc_depth
                    r["local_cd8"] = loc_cd8
                    r["n_local_beads"] = len(nb)
                    per.append(r)
                pf = pd.DataFrame(per)
                n_test = int(pf["testable"].sum())
                # CAPABILITY vs REALIZATION. "testable" means k = 2 is
                # arithmetically achievable given n_local. "realized" means a
                # clone ACTUALLY appeared twice. The first run reported only
                # capability, which overstated: puck 29 arm A showed 32
                # testable but only 9 anchors had any local clone recurrence
                # at all. Both are reported now.
                n_real = int(pf["p_raw"].notna().sum())
                n_hit_raw = int((pf["p_raw"] < ALPHA).sum())
                n_hit_bonf = int((pf["p_bonf"] < ALPHA).sum())
                tag = "  <== PRIMARY" if (R == R0 and arm == "A_cd8_puckwide") else ""
                log(f"    R={R:>5.0f}  anchors {len(pf):>3} | testable {n_test:>3} "
                    f"| realized {n_real:>3} "
                    f"| median n_local {np.nanmedian(pf.n_local):>4.0f} "
                    f"| hits raw {n_hit_raw:>3} | hits bonf {n_hit_bonf:>3}{tag}")
                rows.append({"part": "arm_summary", "puck": p, "arm": arm,
                             "radius": R, "n_anchors": len(pf),
                             "n_testable": n_test, "n_realized": n_real,
                             "median_n_local": float(np.nanmedian(pf.n_local)),
                             "n_hits_raw": n_hit_raw, "n_hits_bonf": n_hit_bonf})

                if R != R0:
                    continue

                # ==========================================================
                # ROTATION CALIBRATION. Not an extra rung: on this statistic
                # it is the ONLY VALID ONE.
                #
                # Measured type I on a random field with no clonal structure
                # whatsoever: raw 0.18-0.72, Bonferroni 0.18-0.49, against a
                # nominal 0.05. Bonferroni fails because its reference set is
                # DATA-DEPENDENT: only clones already present locally with
                # k >= 2 are tested, so the correction is applied over clones
                # that were selected for having shown the signal. The right
                # reference is every clone that COULD have landed there, which
                # has no closed form here.
                #
                # So the per-anchor p is calibrated EMPIRICALLY against the
                # same statistic computed at rotated anchor positions. Rotation
                # keeps the anchor constellation rigid, so inter-anchor overlap
                # (31/17/53 anchors forming 28/13/36 effective groups) survives
                # into the null. p_raw and p_bonf are still reported, but
                # p_emp is the operative one and the hit list is built on it.
                # ==========================================================
                thr_emp = np.nan
                if regional or args.n_perm <= 0:
                    log("      rotation skipped for this arm; p_emp unavailable, "
                        "so its hit counts are NOT calibrated and must not be "
                        "read as the primary arm's are")
                    pf["p_emp"] = np.nan
                else:
                    cxy = xy.mean(axis=0)
                    A = np.array([xy[pos[a]] for a in anchors]) - cxy
                    null = []
                    for s in range(args.n_perm):
                        th = RNG.uniform(0, 2 * np.pi)
                        ref = -1.0 if RNG.random() < 0.5 else 1.0
                        Rm = np.array([[np.cos(th), -np.sin(th)],
                                       [np.sin(th), np.cos(th)]])
                        B = (A * np.array([ref, 1.0])) @ Rm.T + cxy
                        for bxy in B:
                            rr = anchor_test(tree, fx, clone_arr, counts_all,
                                             N_all, bxy, R)
                            if rr and np.isfinite(rr["p_raw"]):
                                null.append(rr["p_raw"])
                    null = np.asarray(null, float)
                    if len(null) < 50:
                        log(f"      rotation produced only {len(null)} null "
                            f"values; p_emp unreliable, reported but not used")
                        pf["p_emp"] = np.nan
                    else:
                        ns = np.sort(null)
                        pf["p_emp"] = [
                            (1 + np.searchsorted(ns, v, side="right")) / (1 + len(ns))
                            if np.isfinite(v) else np.nan
                            for v in pf["p_raw"].values]
                        thr_emp = float(np.quantile(ns, ALPHA))
                        n_hit_emp = int((pf["p_emp"] < ALPHA).sum())
                        exp_hits = ALPHA * int(pf["p_raw"].notna().sum())
                        log(f"      ROTATION NULL  {len(ns):,} null statistics "
                            f"from {args.n_perm} rotations")
                        log(f"        empirical 5% threshold on p_raw: "
                            f"{thr_emp:.3e}")
                        log(f"        raw hits {n_hit_raw} | Bonferroni "
                            f"{n_hit_bonf} | CALIBRATED {n_hit_emp} "
                            f"(null expects {exp_hits:.1f})")
                        log("        If calibrated hits collapse toward the")
                        log("        expectation, the raw and Bonferroni counts")
                        log("        were the k=2 diffusion mechanism, not")
                        log("        biology.")
                        rows.append({"part": "rotation", "puck": p, "arm": arm,
                                     "radius": R, "n_null": len(ns),
                                     "thr_emp": thr_emp,
                                     "n_hits_raw": n_hit_raw,
                                     "n_hits_bonf": n_hit_bonf,
                                     "n_hits_emp": n_hit_emp,
                                     "expected_hits": exp_hits,
                                     "n_perm": args.n_perm})

                        pv = pf["p_emp"].dropna().values
                        if len(pv) >= 8:
                            ks = kstest(pv, "uniform", alternative="greater")
                            log(f"      P-DISTRIBUTION on CALIBRATED p  n={len(pv)}"
                                f" | below 0.1: {(pv<0.1).sum()} "
                                f"(expects {0.1*len(pv):.1f}) | KS p = "
                                f"{ks.pvalue:.4g}")
                            log("        The only rung with power against a")
                            log("        DIFFUSE subset: five anchors at p~0.05")
                            log("        that BH flags at none and test A")
                            log("        dilutes against 48 nulls. Run on p_emp")
                            log("        because p_raw is not uniform under the")
                            log("        null here.")
                            rows.append({"part": "p_distribution", "puck": p,
                                         "arm": arm, "radius": R, "n": len(pv),
                                         "n_below_0.1": int((pv < 0.1).sum()),
                                         "expected_uniform": 0.1 * len(pv),
                                         "ks_p": float(ks.pvalue)})

                pf["p_bh"] = bh(pf["p_emp"].values if pf["p_emp"].notna().any()
                                else pf["p_bonf"].values)
                pf["puck"] = p
                pf["arm"] = arm
                pf["radius"] = R
                hits.append(pf)

                # ---- spike-in power, primary arm only ----------------------
                if arm == "A_cd8_puckwide" and not args.skip_spike \
                        and np.isfinite(thr_emp):
                    log(f"      SPIKE-IN POWER (synthetic clone, k FOCI, "
                        f"{args.n_spike} reps each)")
                    log("        Detection uses the SAME calibrated threshold")
                    log("        as the hit list, so power and type I are on")
                    log("        one scale.")
                    log("        The synthetic clone exists ONLY locally")
                    log("        (K = k), the expansion hypothesis in its most")
                    log("        detectable form. Power against a clone that")
                    log("        also appears elsewhere is LOWER.")
                    log(f"        {'k local':>8} |" + "".join(
                        f"{('K=k' if K is None else f'K={K}'):>9}"
                        for K in SPIKE_K_PUCKWIDE))
                    for k in SPIKE_K:
                        cells = []
                        for K in SPIKE_K_PUCKWIDE:
                            Kuse = k if K is None else max(K, k)
                            det = 0
                            for _ in range(args.n_spike):
                                a = anchors[RNG.integers(len(anchors))]
                                axy = xy[pos[a]]
                                rr = anchor_test(tree, fx, clone_arr,
                                                 counts_all, N_all, axy, R,
                                                 extra_clone=-1, extra_k=k,
                                                 extra_K=Kuse)
                                if rr and rr["clone"] == -1 \
                                        and np.isfinite(rr["p_raw"]) \
                                        and rr["p_raw"] < thr_emp:
                                    det += 1
                            rate = det / max(1, args.n_spike)
                            cells.append(rate)
                            spike_rows.append({"puck": p, "radius": R, "k": k,
                                               "K_puckwide": Kuse,
                                               "K_label": ("local"
                                                           if K is None
                                                           else str(K)),
                                               "detect_rate": rate,
                                               "n_rep": args.n_spike})
                            rows.append({"part": "spike_in", "puck": p,
                                         "arm": arm, "radius": R, "k_foci": k,
                                         "K_puckwide": Kuse,
                                         "detect_rate": rate})
                        log(f"        {k:>8} |" + "".join(
                            f"{100*c:>8.1f}%" for c in cells))
                    log("        ROWS are the clone's LOCAL foci, COLUMNS its")
                    log("        PUCK-WIDE foci. K=k is a clone that exists")
                    log("        only here, the most detectable and least")
                    log("        realistic case. A real expanded clone is")
                    log("        present throughout the tumour: read the")
                    log("        right-hand columns for that.")
    # ----------------------------------------------------------------------
    banner("HIT TABLE  (primary radius, all arms)")
    if hits:
        H = pd.concat(hits, ignore_index=True)
        ann = neo.drop_duplicates(subset=[ncol]).set_index(ncol)
        keep = [c for c in ("gene", "hgvs_p", "mut_peptide", "best_allele",
                            "mut_ic50", "delta_ic50", "is_strong", "is_tcw_ct",
                            "unified_annotation") if c in ann.columns]
        H = H.join(ann[keep], on="anchor")

        # clone geometry and CDR3 consensus for anything that cleared Bonferroni
        extra = []
        hcol = "p_emp" if ("p_emp" in H.columns and H["p_emp"].notna().any()) \
            else "p_bonf"
        for _, r in H[H[hcol] < ALPHA].iterrows():
            sub = tcr[(tcr.puck == r.puck) & (tcr.cloneId == r.clone)]
            if not len(sub):
                continue
            vc = sub["aaSeqCDR3"].astype(str).value_counts()
            extra.append({"puck": r.puck, "arm": r.arm, "anchor": r.anchor,
                          "consensus_cdr3": vc.index[0],
                          "cdr3_stability": float(vc.iloc[0] / vc.sum()),
                          "clone_umis_puck": int(len(sub)),
                          "clone_beads_puck": int(sub.bead.nunique())})
        if extra:
            H = H.merge(pd.DataFrame(extra),
                        on=["puck", "arm", "anchor"], how="left")
        H["candidate_not_finding"] = True
        H.to_csv(os.path.join(OUT, f"per_anchor_results_{args.chain}.tsv"),
                 sep="\t", index=False)

        pcol = "p_emp" if ("p_emp" in H.columns and H["p_emp"].notna().any()) \
            else "p_bonf"
        prim = H[(H.arm == "A_cd8_puckwide") & (H[pcol] < ALPHA)].sort_values(pcol)
        log(f"  Anchors clearing the CALIBRATED threshold ({pcol}) in the "
            f"PRIMARY arm: {len(prim)}")
        if pcol != "p_emp":
            log("  *** p_emp unavailable, falling back to p_bonf, whose measured")
            log("  *** type I on a random field is 0.18-0.49. Treat this list as")
            log("  *** uncalibrated.")
        if len(prim):
            cols = [c for c in ("puck", "anchor", "gene", "hgvs_p",
                                "best_allele", "mut_ic50", "is_tcw_ct",
                                "n_local", "k", "K", "min_detectable_k",
                                "local_depth", "local_cd8",
                                "p_raw", "p_bonf", "p_emp", "p_bh",
                                "cdr3_stability") if c in prim.columns]
            log(prim[cols].to_string(index=False, max_rows=40))
            log()
            exp_tot = sum(r.get("expected_hits", 0) for r in rows
                          if r.get("part") == "rotation"
                          and r.get("arm") == "A_cd8_puckwide")
            log(f"  CONTEXT: {len(prim)} anchors cleared, chance alone expects "
                f"{exp_tot:.1f} across the three pucks.")
            log("  EVERY ROW IS A CANDIDATE FOR FOLLOW-UP, NOT A FINDING.")
            log("  This is the third parallel design on the same 110 anchors and")
            log("  the most permissive of the three. Read it against the")
            log("  rotation null above before believing any single row.")
            log()
            log("  cdr3_stability is the fraction of that clone's rows carrying")
            log("  its most common aaSeqCDR3. It is the GATE on whether any")
            log("  pMHC-TCR pair prediction (pMTnet, NetTCR, ERGO-II) is")
            log("  possible: a low value means the sequence fed to such a tool")
            log("  was not the sequence in the patient.")
            log()
            if "local_depth" in H.columns:
                armA = H[(H.arm == "A_cd8_puckwide") & H.local_depth.notna()]
                hit = armA[armA[pcol] < ALPHA]["local_depth"].values
                non = armA[armA[pcol] >= ALPHA]["local_depth"].values
                if len(hit) >= 2 and len(non) >= 5:
                    sp = np.sqrt((np.var(hit, ddof=1) if len(hit) > 1 else 0.0
                                  + np.var(non, ddof=1)) / 2)
                    d = ((np.mean(hit) - np.mean(non)) / sp
                         if sp > 0 else np.nan)
                    log(f"  DEPTH DIAGNOSTIC  median local depth: hits "
                        f"{np.median(hit):,.0f} vs non-hits "
                        f"{np.median(non):,.0f} | SMD {d:+.3f}")
                    log("    A large positive SMD means the hits are simply the")
                    log("    deepest neighbourhoods and conditioning on n did")
                    log("    not do its job.")
                    rows.append({"part": "depth_diagnostic", "arm": "A_cd8_puckwide",
                                 "median_depth_hits": float(np.median(hit)),
                                 "median_depth_nonhits": float(np.median(non)),
                                 "smd": float(d) if d == d else None,
                                 "n_hits": len(hit), "n_nonhits": len(non)})
        else:
            log("  None. With the rotation null above, that is the answer to the")
            log("  per-anchor question at this depth.")

    if spike_rows:
        S = pd.DataFrame(spike_rows)
        S.to_csv(os.path.join(OUT, "spike_in_power.tsv"), sep="\t", index=False)
        fig, axes = plt.subplots(1, len(PUCKS), figsize=(13 * len(PUCKS), 11),
                                 sharey=True)
        cols = [C_BLUE, C_GREEN, C_MUSTARD, C_CORAL, C_DARK]
        labs = sorted(S.K_label.unique(),
                      key=lambda v: (v != "local", int(v) if v.isdigit() else 0))
        for i, pk in enumerate(PUCKS):
            ax = axes[i] if len(PUCKS) > 1 else axes
            for j, lab in enumerate(labs):
                sl = S[(S.puck == pk) & (S.radius == R0) &
                       (S.K_label == lab)].sort_values("k")
                if len(sl):
                    ax.plot(sl.k, 100 * sl.detect_rate, "-o", lw=3, ms=12,
                            color=cols[j % len(cols)],
                            label=("K = k (local only)" if lab == "local"
                                   else f"K = {lab} puck-wide"))
            ax.axhline(80, color=C_DARK, ls="--", lw=2)
            ax.set_ylim(0, 102)
            style(ax, title=f"Puck {pk}",
                  xlabel="local foci of the clone (k)",
                  ylabel="detection rate (%)" if i == 0 else None)
            if i == 0:
                ax.legend(fontsize=FS_ANNOT, frameon=False, loc="lower right")
        fig.tight_layout()
        savefig(fig, "Fig_08g_A_spikein_power")
        log()
        log("  THIS FIGURE ANSWERS THE CONDITIONAL QUESTION. It converts the")
        log("  result into a sentence: this design detects a locally expanded")
        log("  clone of k or more SPATIALLY SEPARATED foci at X% power. If the")
        log("  smallest k reaching 80% is large, a null result is a statement")
        log("  about the assay, not about the biology, and TCR-binding")
        log("  prediction is the right next investment. If k is small and")
        log("  nothing was found, the biology is the answer.")

    pd.DataFrame(rows).to_csv(
        os.path.join(OUT, f"per_anchor_summary_{args.chain}.tsv"),
        sep="\t", index=False)
    log(f"\noutputs in {OUT}")
    log(f"  finished {datetime.now():%Y-%m-%d %H:%M:%S}")


if __name__ == "__main__":
    main()
