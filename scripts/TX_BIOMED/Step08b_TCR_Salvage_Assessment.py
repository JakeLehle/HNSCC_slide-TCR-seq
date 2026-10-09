#!/usr/bin/env python3
"""
Step08b_TCR_Salvage_Assessment.py   v2 (corrected)
================================================================================
READ-ONLY. Decides whether the processed TCR tables can support spatial
inference, and at what cost.

CORRECTIONS IN v2  (all six were misrepresentations in v1, not cosmetic)

  BEAT 1  v1 printed "mixcr nested within cloneId" whenever no mixcr clone
          spanned two cloneIds. That test passes TRIVIALLY when the mapping is
          one-to-one, which it is. v2 classifies the mapping explicitly as
          one-to-one / nested / crossing. It also notes that cloneId and
          mixcr_cloneId are numbered PER PUCK, so a nunique() over the pooled
          frame collides IDs across pucks and undercounts. The 43,174 figure
          quoted from Step08a was that pooling artifact.

  BEAT 2  v1 fit all four predictors in ONE model. They measure the same
          underlying quantity, so marker_T absorbed the signal and cd8, cd4 and
          is_T took small NEGATIVE residual coefficients. Those negatives were
          read as "CD8 is anti-correlated with TCR", which is wrong: they mean
          "adds nothing once marker_T is known". v2 reports MARGINAL
          (single-predictor) models alongside the joint model, plus the
          predictor correlation matrix, so collinearity is visible.

  BEAT 3  v1 gave a Wilson interval on the RATE but printed enrichment as a
          bare point estimate, so there was no way to see which bands were
          evidence. v2 propagates the interval onto the enrichment and prints
          the lower bound, which is the number that decides it.

  BEAT 4  v1 counted "local pairs" without separating the clone's own SOURCE
          bead from genuine NEIGHBOUR beads. With one source per clone, 1,532
          local pairs across 1,523 clones means only NINE real neighbours. v1
          made a source-bead selection look like a spatial neighbourhood.
          v2 reports source and neighbour counts and UMIs separately.

  BEAT 5  Following from BEAT 4, v1's "LOCAL" was ~99% source beads, so its
          2.1-2.4x enrichment could not be attributed. v2 splits the outcome
          into SOURCE, NEIGHBOUR, DIFFUSE and ALL so the two mechanisms are
          separable. If NEIGHBOUR has too few UMIs to fit, that is itself the
          result and is reported as such.

  BEAT 6  v1's depth matching FAILED in pucks 37 and 40 (median depth 873 vs
          600, 1354 vs 876) because min(n_paired, n_unpaired) per decile left
          the two sides unbalanced. v2 balances both sides per stratum, uses 20
          strata, and prints a QC line that FAILS LOUDLY if matching did not
          work. v1 also reported a RATIO of marker_T means; that variable is
          centred near zero and slightly negative, so the ratio was unstable
          and meaningless. v2 reports differences with Cohen's d and a
          Mann-Whitney p.

  BEAT 8  v2 adds the error quantification that was implicit in v1's output:
          distinct 14-mers recovered from the ONT reads vs real beads on the
          puck, and what fraction of recovered barcodes are in the whitelist.

USAGE
    python Step08b_TCR_Salvage_Assessment.py
    python Step08b_TCR_Salvage_Assessment.py --validate-only
    python Step08b_TCR_Salvage_Assessment.py --skip-ont
    python Step08b_TCR_Salvage_Assessment.py --radius 50

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
================================================================================
"""

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

try:
    import anndata as ad
except ImportError:
    sys.exit("anndata not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTPUTS = os.path.join(PROOT, "data/outputs")
OUT = os.path.join(OUTPUTS, "10_tcr_salvage")
FIGDIR = os.path.join(OUT, "figures")

ADATA = os.path.join(OUTPUTS, "04_annotation/all_pucks_annotated_unified.h5ad")
TCR_DIR = os.path.join(PROOT, "data/inputs/tcr/processed")
ONT_DIR = os.path.join(PROOT, "data/inputs/tcr/ont")
FASTQ_DIR = os.path.join(PROOT, "data/inputs/fastq")

PUCKS = ["29", "37", "40"]
TCR_CSV = {p: f"B59_{p}_hTCR_tcr.csv" for p in PUCKS}
TCR_ONT = {"29": "TCR_20220224_Puck_211214_29.gz",
           "37": "TCR_20220228_Puck_211214_37.gz",
           "40": "TCR_20220127_Puck_211214_40.gz"}

ANNOT_COL = "unified_annotation"
BASAL_T = "T_cell"
MARK = "marker_score_T_cell"
C_CD8 = "c2l_CD8-positive, alpha-beta T cell"
C_CD4 = "c2l_CD4-positive, alpha-beta T cell"

UP_LINKER = "TCTTCAGCGTTCCCGAGA"
UP_LINKER_RC = "TCTCGGGAACGCTGAAGA"
UMI_SPACE = 4 ** 9

LOCAL_RADII = [15.0, 30.0, 50.0, 100.0]
PRIMARY_RADIUS = 30.0
DIST_BINS = np.array([0, 15, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 2500, 5000], float)
MIN_BEADS_FOR_KERNEL = 3
CHAIN_MODES = [("all", None), ("TRB", ["TRB"])]
BANDS = [("3-5", 3, 5), ("6-20", 6, 20), ("21-100", 21, 100), ("101+", 101, 10**9)]
N_STRATA = 20                 # BEAT 6 depth strata; v1 used 10 and under-matched
MIN_UMI_FOR_GLM = 50

FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
C_CORAL = "#ed6a5a"
C_MUSTARD = "#F6D155"
C_GRAY = "#d3d3d3"
C_DARK = "#4a4a4a"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

RNG = np.random.default_rng(20261006)
COMP = str.maketrans("ACGTN", "TGCAN")


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


def rc(s):
    return s.translate(COMP)[::-1]


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z*z/n
    c = p + z*z/(2*n)
    s = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return (max(0.0, (c - s)/d), (c + s)/d)


def cohens_d(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    sp = np.sqrt(((na-1)*np.var(a, ddof=1) + (nb-1)*np.var(b, ddof=1)) / (na+nb-2))
    return (np.mean(a) - np.mean(b)) / sp if sp else np.nan


def load_tcr(p):
    df = pd.read_csv(os.path.join(TCR_DIR, TCR_CSV[p]), sep=None, engine="python")
    df.columns = [c.strip() for c in df.columns]
    df["puck"] = p
    df["bead"] = df["bc"].astype(str) + f"-1_Puck_211214_{p}"
    df["umi_purity"] = np.where(df.n_reads > 0, df.n_reads_clone / df.n_reads, np.nan)
    return df


def subset_chain(tcr, chains):
    return tcr if chains is None else tcr[tcr.chain.isin(chains)]


def puck_frame(adata, p):
    m = adata.obs["puck_id"].astype(str).values == f"Puck_211214_{p}"
    idx = adata.obs_names[m]
    f = pd.DataFrame(index=idx)
    f["depth"] = adata.obs.loc[idx, "total_counts"].values.astype(float)
    f["annot"] = adata.obs.loc[idx, ANNOT_COL].astype(str).values
    f["is_T"] = (f["annot"] == BASAL_T).astype(int)
    for col, nm in ((MARK, "marker_T"), (C_CD8, "cd8"), (C_CD4, "cd4")):
        if col in adata.obs.columns:
            f[nm] = adata.obs.loc[idx, col].values.astype(float)
    xy = adata.obs.loc[idx, ["x_coord", "y_coord"]].values.astype(float)
    return idx, f, xy, {b: i for i, b in enumerate(idx)}


def fit_glm(y, cols, frame, offset):
    """Poisson GLM, log-depth offset. Standardized beta = beta * SD(predictor),
    which is what makes a near-z-score, an abundance and a binary comparable."""
    try:
        import statsmodels.api as sm
    except ImportError:
        return None
    X = frame[cols].values.astype(float)
    sds = X.std(axis=0, ddof=1)
    Xc = sm.add_constant(X, has_constant="add")
    try:
        res = sm.GLM(y, Xc, family=sm.families.Poisson(), offset=offset).fit()
    except Exception:
        return None
    out = {}
    for j, c in enumerate(cols, start=1):
        b = float(res.params[j])
        out[c] = {"beta": b, "sd": float(sds[j-1]), "beta_std": b*float(sds[j-1]),
                  "exp_beta_std": float(np.exp(b*sds[j-1])), "p": float(res.pvalues[j])}
    return out


# ==============================================================================
def beat1_cloneid(tcr, rows):
    banner("BEAT 1  cloneId vs mixcr_cloneId")
    log("  v2 CORRECTION. v1's nesting test passed trivially when the mapping is")
    log("  one-to-one, and printed 'mixcr nested within cloneId' regardless.")
    log("  Both IDs are numbered PER PUCK, so a nunique() over the pooled frame")
    log("  collides them across pucks. The 43,174 figure quoted from Step08a was")
    log("  that artifact, not evidence of merging.")

    tot_c = sum(tcr[tcr.puck == p].cloneId.nunique() for p in PUCKS)
    tot_m = sum(tcr[tcr.puck == p].mixcr_cloneId.nunique() for p in PUCKS)
    log(f"\n  pooled nunique (WRONG, IDs collide): cloneId "
        f"{tcr.cloneId.nunique():,}, mixcr {tcr.mixcr_cloneId.nunique():,}")
    log(f"  sum of per-puck nunique (correct):   cloneId {tot_c:,}, mixcr {tot_m:,}")

    for p in PUCKS:
        d = tcr[tcr.puck == p]
        n_c, n_m = d.cloneId.nunique(), d.mixcr_cloneId.nunique()
        per_clone = d.groupby("cloneId")["mixcr_cloneId"].nunique()
        per_mixcr = d.groupby("mixcr_cloneId")["cloneId"].nunique()
        c_multi, m_multi = int((per_clone > 1).sum()), int((per_mixcr > 1).sum())

        if c_multi == 0 and m_multi == 0 and n_c == n_m:
            verdict = "ONE-TO-ONE (the two columns are the same partition)"
        elif m_multi == 0 and c_multi > 0:
            verdict = "mixcr NESTED within cloneId (cloneId is coarser, a real merge)"
        elif c_multi == 0 and m_multi > 0:
            verdict = "cloneId NESTED within mixcr (mixcr is coarser)"
        else:
            verdict = "CROSSING (neither nests in the other)"

        log(f"\n  puck {p}: {n_c:,} cloneId, {n_m:,} mixcr_cloneId")
        log(f"    cloneId spanning >1 mixcr: {c_multi:,} (max {int(per_clone.max())})")
        log(f"    mixcr spanning >1 cloneId: {m_multi:,} (max {int(per_mixcr.max())})")
        log(f"    VERDICT: {verdict}")
        if c_multi:
            big = d.groupby("cloneId").agg(nb=("bead", "nunique")).nlargest(5, "nb")
            log("    top 5 cloneId by beads, and their mixcr split:")
            for cid, r in big.iterrows():
                sub = d[d.cloneId == cid].groupby("mixcr_cloneId")["bead"].nunique() \
                        .sort_values(ascending=False)
                log(f"      cloneId {cid}: {int(r.nb):,} beads -> {len(sub)} mixcr, "
                    f"largest {int(sub.iloc[0]):,} ({100*sub.iloc[0]/r.nb:.0f}%)")
        else:
            log("    no merging: giant clones are genuine single MiXCR clones,")
            log("    so clone size is not an aggregation artifact")
        rows.append({"beat": "cloneid", "puck": p, "n_cloneId": n_c, "n_mixcr": n_m,
                     "cloneId_spanning_multi_mixcr": c_multi,
                     "mixcr_spanning_multi_cloneId": m_multi, "verdict": verdict})


# ==============================================================================
def beat2_standardized(tcr, adata, rows):
    banner("BEAT 2  effect sizes: MARGINAL and JOINT")
    log("  v2 CORRECTION. v1 fit all predictors in one model. They measure the")
    log("  same underlying quantity, so marker_T absorbed the signal and the")
    log("  others took small NEGATIVE residual coefficients. Those negatives do")
    log("  NOT mean the predictor is anti-correlated with TCR; they mean it adds")
    log("  nothing once marker_T is known. Marginal models are shown first.")

    for mode, chains in CHAIN_MODES:
        log(f"\n  ===== chain mode: {mode} =====")
        t = subset_chain(tcr, chains)
        for p in PUCKS:
            idx, f, xy, pos = puck_frame(adata, p)
            u = t[t.puck == p].groupby("bead")["umi"].nunique()
            y = u.reindex(f.index).fillna(0).values.astype(float)
            off = np.log(np.clip(f["depth"].values, 1, None))
            cols = [c for c in ("marker_T", "cd8", "cd4") if c in f.columns] + ["is_T"]

            log(f"\n    puck {p}  ({mode}): {int((y > 0).sum()):,} TCR-positive beads, "
                f"{y.sum():,.0f} UMIs")

            cm = f[cols].corr()
            log("      predictor correlations (why the joint model behaves as it does):")
            log("        " + "".join(f"{c:>10}" for c in cols))
            for c in cols:
                log(f"        {c:<10}" + "".join(f"{cm.loc[c, c2]:>10.3f}" for c2 in cols))

            log(f"      {'predictor':<11} {'MARGINAL':>22} | {'JOINT':>22}")
            log(f"      {'':<11} {'b*SD':>9} {'exp':>6} {'p':>6} | "
                f"{'b*SD':>9} {'exp':>6} {'p':>6}")
            joint = fit_glm(y, cols, f, off)
            for c in cols:
                marg = fit_glm(y, [c], f, off)
                if marg is None or joint is None:
                    log(f"      {c:<11}  fit failed"); continue
                m_, j_ = marg[c], joint[c]
                log(f"      {c:<11} {m_['beta_std']:>+9.4f} {m_['exp_beta_std']:>6.3f} "
                    f"{m_['p']:>6.0e} | {j_['beta_std']:>+9.4f} "
                    f"{j_['exp_beta_std']:>6.3f} {j_['p']:>6.0e}")
                rows.append({"beat": "effects", "puck": p, "chain_mode": mode,
                             "predictor": c,
                             "marginal_beta_std": round(m_["beta_std"], 5),
                             "marginal_exp": round(m_["exp_beta_std"], 4),
                             "marginal_p": m_["p"],
                             "joint_beta_std": round(j_["beta_std"], 5),
                             "joint_exp": round(j_["exp_beta_std"], 4),
                             "joint_p": j_["p"], "sd": round(j_["sd"], 4)})
            if joint:
                best = max(cols, key=lambda c: abs(joint[c]["beta_std"]))
                log(f"      largest joint standardized effect: {best}")
    log("\n  Read MARGINAL to ask 'does this measure track TCR signal at all'.")
    log("  Read JOINT to ask 'does it add anything beyond the others'.")


# ==============================================================================
def clone_geometry(tcr, adata, p, chains, min_beads=MIN_BEADS_FOR_KERNEL):
    idx, f, xy, pos = puck_frame(adata, p)
    d = subset_chain(tcr[tcr.puck == p], chains)
    d = d[d.bead.isin(set(idx))]
    recs = []
    for cid, g in d.groupby("cloneId"):
        pb = g.groupby("bead").agg(reads=("n_reads", "sum"), umi=("umi", "nunique"))
        if len(pb) < min_beads:
            continue
        src = pb["reads"].idxmax()
        si = pos[src]
        ii = [pos[b] for b in pb.index]
        dist = np.hypot(xy[ii, 0] - xy[si, 0], xy[ii, 1] - xy[si, 1])
        for (b, r), dd in zip(pb.iterrows(), dist):
            recs.append((cid, b, float(dd), int(r.umi), int(r.reads),
                         len(pb), b == src))
    return (pd.DataFrame(recs, columns=["cloneId", "bead", "dist_um", "umi",
                                        "reads", "clone_beads", "is_source"]),
            idx, f, xy, pos)


def beat3_nearfield(tcr, adata, rows):
    banner("BEAT 3  near-field enrichment with counts and enrichment intervals")
    log("  v2 CORRECTION. v1 gave a Wilson interval on the RATE but printed")
    log("  enrichment as a bare point estimate, so there was no way to see which")
    log("  bands were evidence. The interval is now propagated onto enrichment.")
    log("  The LOWER BOUND is the number that decides it: above 1.0 means the")
    log("  band is evidence of near-field structure, spanning 1.0 means it is not.")

    for mode, chains in CHAIN_MODES:
        log(f"\n  ===== chain mode: {mode} =====")
        for p in PUCKS:
            geo, idx, f, xy, pos = clone_geometry(tcr, adata, p, chains)
            if not len(geo):
                log(f"    puck {p}: no qualifying clones"); continue
            sat = geo[~geo.is_source]
            log(f"\n    puck {p}  ({mode})   total non-source beads within 15 um of "
                f"their source: {int((sat.dist_um <= 15).sum()):,}")
            log(f"      {'band':>7} {'clones':>7} {'obs':>5} {'avail':>8} "
                f"{'rate':>10} {'enrich':>8} {'enrich 95% CI':>22} {'verdict':>10}")
            for bn, lo, hi in BANDS:
                s = sat[(sat.clone_beads >= lo) & (sat.clone_beads <= hi)]
                srcs = geo[(geo.is_source) & (geo.clone_beads >= lo) &
                           (geo.clone_beads <= hi)]
                if not len(srcs):
                    continue
                avail = np.zeros(len(DIST_BINS) - 1)
                for b in srcs.bead:
                    si = pos[b]
                    dall = np.hypot(xy[:, 0] - xy[si, 0], xy[:, 1] - xy[si, 1])
                    avail += np.histogram(np.delete(dall, si), bins=DIST_BINS)[0]
                obs = np.histogram(s.dist_um.values, bins=DIST_BINS)[0]
                base = obs.sum() / avail.sum() if avail.sum() else np.nan
                k, n = int(obs[0]), int(avail[0])
                rate = k/n if n else np.nan
                lo_r, hi_r = wilson(k, n)
                if base and base == base and n:
                    e, e_lo, e_hi = rate/base, lo_r/base, hi_r/base
                    verdict = "EVIDENCE" if e_lo > 1.0 else "noise"
                else:
                    e = e_lo = e_hi = np.nan
                    verdict = "n/a"
                log(f"      {bn:>7} {len(srcs):>7,} {k:>5,} {n:>8,} {rate:>10.6f} "
                    f"{e:>8.2f} [{e_lo:>8.2f},{e_hi:>9.2f}] {verdict:>10}")
                rows.append({"beat": "nearfield", "puck": p, "chain_mode": mode,
                             "band": bn, "clones": len(srcs), "obs_0_15um": k,
                             "avail_0_15um": n,
                             "enrichment": float(e) if e == e else None,
                             "enrichment_lo": float(e_lo) if e_lo == e_lo else None,
                             "enrichment_hi": float(e_hi) if e_hi == e_hi else None,
                             "verdict": verdict})
    log("\n  The enrichment CI treats the baseline as fixed, so it is slightly")
    log("  optimistic. A band that fails here fails comfortably.")


# ==============================================================================
def beat4_decompose(tcr, adata, rows):
    banner("BEAT 4  decomposition: SOURCE beads vs genuine NEIGHBOUR beads")
    log("  v2 CORRECTION. v1 counted 'local pairs' without separating a clone's")
    log("  own SOURCE bead from genuine NEIGHBOURS. With one source per clone,")
    log("  1,532 local pairs across 1,523 clones means only NINE real neighbours.")
    log("  v1 made a source-bead SELECTION look like a spatial NEIGHBOURHOOD.")
    store = {}
    for mode, chains in CHAIN_MODES:
        log(f"\n  ===== chain mode: {mode} =====")
        for p in PUCKS:
            geo, idx, f, xy, pos = clone_geometry(tcr, adata, p, chains)
            if not len(geo):
                continue
            store[(mode, p)] = (geo, idx, f, xy, pos)
            nclone = geo.cloneId.nunique()
            src = geo[geo.is_source]
            sat = geo[~geo.is_source]
            tot = geo.umi.sum()
            log(f"\n    puck {p}  ({mode}): {nclone:,} clones, {len(geo):,} "
                f"clone-bead pairs, {tot:,} UMIs")
            log(f"      source beads: {len(src):,} carrying {int(src.umi.sum()):,} "
                f"UMIs ({100*src.umi.sum()/tot:.2f}%)")
            log(f"      {'R (um)':>7} {'neighbours':>11} {'nbr UMIs':>10} "
                f"{'% UMIs':>8} {'nbr per clone':>14}")
            for R in LOCAL_RADII:
                nb = sat[sat.dist_um <= R]
                log(f"      {R:>7.0f} {len(nb):>11,} {int(nb.umi.sum()):>10,} "
                    f"{100*nb.umi.sum()/tot:>7.2f}% {len(nb)/nclone:>14.4f}")
                rows.append({"beat": "decompose", "puck": p, "chain_mode": mode,
                             "radius_um": R, "n_clones": nclone,
                             "n_source_beads": len(src),
                             "source_umis": int(src.umi.sum()),
                             "n_neighbour_beads": len(nb),
                             "neighbour_umis": int(nb.umi.sum()),
                             "neighbours_per_clone": round(len(nb)/nclone, 5)})
            nb30 = sat[sat.dist_um <= PRIMARY_RADIUS]
            log(f"      -> at R={PRIMARY_RADIUS:.0f} um, 'local' is "
                f"{100*len(src)/max(1, len(src)+len(nb30)):.1f}% source beads")
    log("\n  If neighbours per clone is far below 1, there is no neighbourhood to")
    log("  filter on and any 'locality' effect is source-bead selection.")
    return store


# ==============================================================================
def beat5_salvage(store, rows):
    banner("BEAT 5  THE SALVAGE TEST, four outcomes")
    log("  v2 CORRECTION. v1's LOCAL was ~99% source beads, so its 2.1-2.4x")
    log("  enrichment could not be attributed to locality rather than to picking")
    log("  each clone's best-supported bead. The outcome is now split four ways:")
    log("    ALL        every TCR UMI on the bead")
    log("    SOURCE     UMIs from clones where this bead IS the source")
    log("    NEIGHBOUR  UMIs from clones where this bead is within R, not source")
    log("    DIFFUSE    the remainder")
    log("  SOURCE high with NEIGHBOUR unfittable means the effect is selection,")
    log("  not spatial proximity. That is a usable result but a different claim.")

    results = []
    for mode, _ in CHAIN_MODES:
        log(f"\n  ===== chain mode: {mode} =====")
        for p in PUCKS:
            key = (mode, p)
            if key not in store:
                continue
            geo, idx, f, xy, pos = store[key]
            off = np.log(np.clip(f["depth"].values, 1, None))
            cols = [c for c in ("marker_T", "cd8") if c in f.columns] + ["is_T"]

            src = geo[geo.is_source]
            nbr = geo[(~geo.is_source) & (geo.dist_um <= PRIMARY_RADIUS)]
            dif = geo[(~geo.is_source) & (geo.dist_um > PRIMARY_RADIUS)]

            def agg(fr):
                return fr.groupby("bead")["umi"].sum().reindex(f.index) \
                         .fillna(0).values.astype(float)

            outs = [("ALL", agg(geo)), ("SOURCE", agg(src)),
                    ("NEIGHBOUR", agg(nbr)), ("DIFFUSE", agg(dif))]
            log(f"\n    puck {p}  ({mode})  R = {PRIMARY_RADIUS:.0f} um")
            log("      UMIs: " + " | ".join(f"{nm} {y.sum():,.0f}" for nm, y in outs))
            log(f"      {'outcome':>10} {'predictor':<10} {'b*SD':>9} "
                f"{'exp(b*SD)':>10} {'p':>11}")
            for nm, y in outs:
                if y.sum() < MIN_UMI_FOR_GLM:
                    log(f"      {nm:>10} only {y.sum():,.0f} UMIs, below the "
                        f"{MIN_UMI_FOR_GLM} floor: NOT FITTABLE")
                    log(f"                 ^ this is the result, not a failure: there")
                    log(f"                   is essentially no neighbour signal to test")
                    rows.append({"beat": "salvage", "puck": p, "chain_mode": mode,
                                 "outcome": nm, "umis": float(y.sum()),
                                 "fittable": False})
                    continue
                res = fit_glm(y, cols, f, off)
                if res is None:
                    continue
                for c in cols:
                    r = res[c]
                    log(f"      {nm:>10} {c:<10} {r['beta_std']:>+9.4f} "
                        f"{r['exp_beta_std']:>10.3f} {r['p']:>11.2e}")
                    results.append({"puck": p, "chain_mode": mode, "outcome": nm,
                                    "predictor": c,
                                    "exp_beta_std": r["exp_beta_std"]})
                    rows.append({"beat": "salvage", "puck": p, "chain_mode": mode,
                                 "outcome": nm, "predictor": c, "fittable": True,
                                 "exp_beta_std": round(r["exp_beta_std"], 4),
                                 "p": r["p"], "umis": float(y.sum())})

            for nm, y in outs:
                if y.sum() < MIN_UMI_FOR_GLM:
                    continue
                dec = pd.qcut(f["depth"], 10, labels=False, duplicates="drop")
                rr = []
                for dd in np.unique(dec):
                    m = dec.values == dd
                    a, b = y[m & (f.is_T.values == 1)], y[m & (f.is_T.values == 0)]
                    if len(a) and len(b) and b.mean():
                        rr.append(a.mean()/b.mean())
                if rr:
                    log(f"      model-free {nm}: depth-stratified T/other ratio "
                        f"{np.nanmean(rr):.3f}")
                    rows.append({"beat": "salvage_modelfree", "puck": p,
                                 "chain_mode": mode, "outcome": nm,
                                 "stratified_ratio": round(float(np.nanmean(rr)), 4)})

    if results:
        rdf = pd.DataFrame(results)
        rdf.to_csv(os.path.join(OUT, "salvage_effects.tsv"), sep="\t", index=False)
        fig, axes = plt.subplots(1, len(CHAIN_MODES),
                                 figsize=(16*len(CHAIN_MODES), 10), sharey=True)
        axes = np.atleast_1d(axes)
        order = ["DIFFUSE", "ALL", "SOURCE", "NEIGHBOUR"]
        cols_ = [C_GRAY, C_MUSTARD, C_CORAL, C_DARK]
        for ax, (mode, _) in zip(axes, CHAIN_MODES):
            s = rdf[(rdf.chain_mode == mode) & (rdf.predictor == "marker_T")]
            w = 0.2
            for i, pk in enumerate(PUCKS):
                for j, o in enumerate(order):
                    v = s[(s.puck == pk) & (s.outcome == o)]["exp_beta_std"]
                    if len(v):
                        ax.bar(i + (j-1.5)*w, float(v.iloc[0]), width=w,
                               color=cols_[j], label=o if i == 0 else None)
            ax.axhline(1.0, color=C_DARK, linestyle="--", linewidth=2)
            ax.set_xticks(range(len(PUCKS)))
            ax.set_xticklabels([f"Puck {p}" for p in PUCKS])
            style(ax, title=f"{mode} chains",
                  ylabel="exp(beta × SD), marker_score_T_cell" if mode == "all" else None)
            ax.legend(fontsize=FS_ANNOT, frameon=False)
        fig.tight_layout(); savefig(fig, "Fig_Salvage_A_outcome_decomposition")

    log("\n  VERDICT GUIDE")
    log("    SOURCE >> ALL, NEIGHBOUR unfittable -> selection effect. Usable, but")
    log("      Methods must say 'best-supported bead per clone', not 'locality'.")
    log("    NEIGHBOUR fittable and >= SOURCE     -> genuine spatial locality.")
    log("    SOURCE ~ ALL                         -> nothing is concentrated.")


# ==============================================================================
def beat6_pairing(tcr, adata, rows):
    banner("BEAT 6  chain pairing, balanced depth matching")
    log("  v2 CORRECTIONS. (1) v1's matching FAILED in pucks 37 and 40, leaving")
    log("  median depth 873 vs 600 and 1354 vs 876, because min(n_paired,")
    log("  n_unpaired) per decile kept all paired beads while sampling fewer")
    log("  unpaired. v2 balances BOTH sides per stratum, uses 20 strata, and")
    log("  prints a QC line that fails loudly if matching did not work.")
    log("  (2) v1 reported a RATIO of marker_T means. That variable is centred")
    log("  near zero and slightly negative, so the ratio was unstable and")
    log("  meaningless. v2 reports differences, Cohen's d and Mann-Whitney.")

    for p in PUCKS:
        idx, f, xy, pos = puck_frame(adata, p)
        d = tcr[(tcr.puck == p) & (tcr.bead.isin(set(idx)))]
        per = d.groupby("bead")["chain"].apply(set)
        has_a = set(b for b, s in per.items() if "TRA" in s)
        has_b = set(b for b, s in per.items() if "TRB" in s)
        both = has_a & has_b
        N = len(idx)
        pa, pb = len(has_a)/N, len(has_b)/N
        obs, exp = len(both)/N, (len(has_a)/N)*(len(has_b)/N)
        ratio = obs/exp if exp else np.nan

        log(f"\n  puck {p}: {N:,} annotated beads")
        log(f"    P(TRA) {pa:.4f} | P(TRB) {pb:.4f} | P(both) {obs:.4f} vs "
            f"independent {exp:.4f} -> ratio {ratio:.3f}")
        if ratio > 1.15:
            log("    Above independence: the two chains co-occur more than chance,")
            log("    consistent with some beads holding a real T cell.")
        else:
            log("    At or near independence: chains arrive as unrelated events.")

        f["paired"] = [1 if b in both else 0 for b in f.index]
        f["tcr_pos"] = [1 if b in set(per.index) else 0 for b in f.index]
        pool = f[f.tcr_pos == 1].copy()
        if pool.paired.sum() < 50 or (pool.paired == 0).sum() < 50:
            log("    too few beads on one side to match"); continue
        pool["stratum"] = pd.qcut(pool["depth"], N_STRATA, labels=False,
                                  duplicates="drop")

        keep_p, keep_u = [], []
        dropped = 0
        for st, g in pool.groupby("stratum"):
            gp, gu = g[g.paired == 1], g[g.paired == 0]
            k = min(len(gp), len(gu))
            if k == 0:
                dropped += len(gp)
                continue
            keep_p.append(gp.sample(k, random_state=int(RNG.integers(1e9))))
            keep_u.append(gu.sample(k, random_state=int(RNG.integers(1e9))))
        if not keep_p:
            log("    matching produced no strata"); continue
        pr, mm = pd.concat(keep_p), pd.concat(keep_u)

        mdp, mdu = pr.depth.median(), mm.depth.median()
        qc = abs(mdp - mdu) / max(1.0, (mdp + mdu)/2)
        ok = (len(pr) == len(mm)) and qc < 0.05
        log(f"    depth-matched: {len(pr):,} paired vs {len(mm):,} unpaired "
            f"({dropped:,} paired beads dropped, no unpaired at that depth)")
        log(f"    MATCHING QC: median depth {mdp:,.0f} vs {mdu:,.0f}, "
            f"relative gap {100*qc:.2f}%  -> {'PASS' if ok else '*** FAIL, do not interpret below ***'}")

        log(f"      {'metric':<18} {'paired':>10} {'unpaired':>10} "
            f"{'diff':>9} {'d':>7} {'p':>10}")
        for col, nm in (("marker_T", "marker_T"), ("cd8", "c2l CD8"),
                        ("cd4", "c2l CD4"), ("is_T", "frac T_cell")):
            if col not in f.columns:
                continue
            a, b = pr[col].values, mm[col].values
            diff = float(np.mean(a) - np.mean(b))
            dd = cohens_d(a, b)
            try:
                _, pv = mannwhitneyu(a, b, alternative="two-sided")
            except ValueError:
                pv = np.nan
            log(f"      {nm:<18} {np.mean(a):>10.4f} {np.mean(b):>10.4f} "
                f"{diff:>+9.4f} {dd:>7.3f} {pv:>10.2e}")
            rows.append({"beat": "pairing", "puck": p, "metric": nm,
                         "paired_mean": float(np.mean(a)),
                         "unpaired_mean": float(np.mean(b)),
                         "difference": round(diff, 5),
                         "cohens_d": round(float(dd), 4) if dd == dd else None,
                         "p": float(pv) if pv == pv else None,
                         "matching_ok": bool(ok)})
        log("      Cohen's d below ~0.1 is negligible even at tiny p, because n")
        log("      is in the thousands. Judge by d, not by p.")

        rows.append({"beat": "pairing_cooccurrence", "puck": p,
                     "P_TRA": round(pa, 5), "P_TRB": round(pb, 5),
                     "P_both_obs": round(obs, 5), "P_both_indep": round(exp, 5),
                     "ratio_obs_over_indep": round(float(ratio), 4),
                     "n_paired": len(both)})

    log("\n  minor chains (retained, reported separately, never folded into TRA):")
    for p in PUCKS:
        ch = tcr[tcr.puck == p]["chain"].value_counts()
        minor = {k: int(v) for k, v in ch.items() if k not in ("TRA", "TRB")}
        log(f"    puck {p}: " + (", ".join(f"{k}:{v}" for k, v in minor.items()) or "none"))
    log("    TRD sits inside the TRA locus between TRAV and TRAJ, so TRA primers")
    log("    can amplify TRD by proximity.")


# ==============================================================================
def beat7_umi_sharing(tcr, rows):
    banner("BEAT 7  UMI sharing across beads within a clone")
    log("  The same UMI on several beads of one clone is an index-hopping or")
    log(f"  template-switching signature. A 9 nt UMI gives {UMI_SPACE:,} options,")
    log("  so chance collisions within a clone of n rows are about n^2/(2*space).")
    for p in PUCKS:
        d = tcr[tcr.puck == p]
        obs_s = exp_s = 0
        clones = 0
        for _, g in d.groupby("cloneId"):
            n = len(g)
            if n < 2:
                continue
            clones += 1
            obs_s += int((g.groupby("umi")["bead"].nunique() > 1).sum())
            exp_s += n*n / (2.0*UMI_SPACE)
        ratio = obs_s/exp_s if exp_s else np.nan
        log(f"  puck {p}: {clones:,} multi-row clones | shared observed {obs_s:,} "
            f"vs expected {exp_s:,.1f} -> {ratio:,.1f}x")
        rows.append({"beat": "umi_sharing", "puck": p, "clones": clones,
                     "shared_observed": obs_s, "shared_expected": round(exp_s, 2),
                     "ratio": round(float(ratio), 3) if ratio == ratio else None})
    log("\n  UPPER BOUND only: real UMI usage is not uniform, so chance sharing")
    log("  runs above this model even with no hopping. A few-fold excess is weak")
    log("  evidence; the absolute counts here are small either way.")


# ==============================================================================
def load_whitelist(p):
    """barcode_matching column 2 (corrected), suffix stripped: the real beads."""
    path = os.path.join(FASTQ_DIR, f"2022-01-28_Puck_211214_{p}",
                        "barcode_matching", f"Puck_211214_{p}_barcode_matching.txt.gz")
    if not os.path.exists(path):
        return set()
    s = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                s.add(f[1].split("-")[0])
    return s


def beat8_ont(tcr, rows, n_reads):
    banner("BEAT 8  ONT provenance and barcode error rate")
    log("  v2 adds the error quantification implicit in v1's output: distinct")
    log("  14-mers recovered from the ONT reads against the number of REAL beads")
    log("  on the puck. Recovering more distinct barcodes than beads exist means")
    log("  the excess is error-derived, which is the misassignment mechanism.")
    for p in PUCKS:
        path = os.path.join(ONT_DIR, TCR_ONT[p])
        if not os.path.exists(path):
            log(f"  puck {p}: {path} not found, skipping"); continue
        csv_bc = set(tcr.loc[tcr.puck == p, "bc"].astype(str))
        wl = load_whitelist(p)
        found, n, fwd, rev = set(), 0, 0, 0
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:
                    continue
                s = line.strip()
                n += 1
                pos = s.find(UP_LINKER)
                if pos >= 8 and pos + 24 <= len(s):
                    fwd += 1
                    bc = s[pos-8:pos] + s[pos+18:pos+24]
                    if len(bc) == 14:
                        found.add(bc)
                else:
                    pos = s.find(UP_LINKER_RC)
                    if pos >= 0:
                        rev += 1
                        seg = rc(s)
                        p2 = seg.find(UP_LINKER)
                        if p2 >= 8 and p2 + 24 <= len(seg):
                            bc = seg[p2-8:p2] + seg[p2+18:p2+24]
                            if len(bc) == 14:
                                found.add(bc)
                if n >= n_reads:
                    break
        inter = csv_bc & found
        log(f"\n  puck {p}: {n:,} reads sampled")
        log(f"    linker: forward {fwd:,} ({100*fwd/n:.2f}%), "
            f"reverse {rev:,} ({100*rev/n:.2f}%)")
        log(f"    distinct 14-mers recovered: {len(found):,}")
        log(f"    CSV barcodes for this puck:  {len(csv_bc):,}")
        log(f"    overlap with CSV: {len(inter):,} "
            f"({100*len(inter)/max(1, len(csv_bc)):.2f}% of CSV barcodes)")
        if wl:
            real = found & wl
            log(f"    real beads on this puck (whitelist): {len(wl):,}")
            log(f"    recovered 14-mers in the whitelist: {len(real):,} "
                f"({100*len(real)/max(1, len(found)):.2f}%)")
            log(f"    recovered 14-mers NOT in whitelist: {len(found)-len(real):,} "
                f"({100*(len(found)-len(real))/max(1, len(found)):.2f}%) <- error derived")
            log(f"    ratio recovered / real beads: {len(found)/max(1, len(wl)):.3f}")
        rows.append({"beat": "ont_provenance", "puck": p, "reads_sampled": n,
                     "linker_fwd": fwd, "linker_rev": rev,
                     "barcodes_recovered": len(found), "csv_barcodes": len(csv_bc),
                     "overlap": len(inter), "whitelist_size": len(wl),
                     "recovered_in_whitelist": len(found & wl) if wl else None,
                     "pct_csv_found": round(100*len(inter)/max(1, len(csv_bc)), 3)})
    log("\n  A substantial overlap with the CSV confirms these reads are the")
    log("  source. The off-whitelist fraction is the per-read barcode error our")
    log("  own calling would have to handle.")


# ==============================================================================
def main():
    global PRIMARY_RADIUS
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--skip-ont", action="store_true")
    ap.add_argument("--ont-reads", type=int, default=500000)
    ap.add_argument("--radius", type=float, default=PRIMARY_RADIUS)
    args = ap.parse_args()
    PRIMARY_RADIUS = args.radius

    os.makedirs(OUT, exist_ok=True)
    needed = [ADATA] + [os.path.join(TCR_DIR, TCR_CSV[p]) for p in PUCKS]
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

    log("loading TCR tables")
    tcr = pd.concat([load_tcr(p) for p in PUCKS], ignore_index=True)
    ch = tcr.chain.value_counts()
    log(f"  {len(tcr):,} rows")
    log(f"  chain census: " + "; ".join(f"{k}:{v:,}" for k, v in ch.items()))
    if "TRA" in ch and "TRB" in ch:
        log(f"  TRB/TRA ratio {ch['TRB']/ch['TRA']:.2f} "
            "(TRA capture is the weaker measurement)")
    log(f"  primary locality radius: {PRIMARY_RADIUS:.0f} um")

    rows = []
    beat1_cloneid(tcr, rows)
    beat2_standardized(tcr, adata, rows)
    beat3_nearfield(tcr, adata, rows)
    store = beat4_decompose(tcr, adata, rows)
    beat5_salvage(store, rows)
    beat6_pairing(tcr, adata, rows)
    beat7_umi_sharing(tcr, rows)
    if not args.skip_ont:
        beat8_ont(tcr, rows, args.ont_reads)
    else:
        log("\n  BEAT 8 skipped by flag")

    out = os.path.join(OUT, "salvage_summary.tsv")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    log(f"\nwrote {out}")

    banner("WHAT THIS DECIDES")
    log("  BEAT 4 first: how many genuine NEIGHBOUR beads exist per clone. If")
    log("  that is far below 1, there is no neighbourhood to filter on.")
    log("  BEAT 5 then: SOURCE vs NEIGHBOUR vs ALL. SOURCE high with NEIGHBOUR")
    log("  unfittable means the effect is best-bead selection, not locality.")
    log("  BEAT 2 marginal column: which T-cell measure tracks TCR signal.")
    log("  BEAT 3 lower bound: which bands are evidence.")
    log("  BEAT 6 QC line: whether the pairing comparison can be read at all.")


if __name__ == "__main__":
    main()
