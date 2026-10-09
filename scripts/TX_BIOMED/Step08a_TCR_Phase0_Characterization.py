#!/usr/bin/env python3
"""
Step08a_TCR_Phase0_Characterization.py
================================================================================
READ-ONLY characterization of the processed TCR tables. No hypothesis testing,
no filtering applied. Its output is what the Phase 1 and Phase 2 thresholds get
decided from.

PREFLIGHT (hard gate)
    The 110 neoantigen bead IDs are reconstructions from restore_cb() in
    Step05c. Membership in adata.obs_names is necessary but not sufficient: a
    bad restore most often yields a VALID barcode belonging to a DIFFERENT real
    bead, which membership testing passes. The gate therefore also compares the
    x_coord and y_coord carried in neoantigens_per_bead.tsv against adata.

BEATS
     1  join validation, bc -> adata.obs_names
     2  cell-type composition of TCR-bearing beads vs puck background
     3  UMI error inflation, adjacency within (bead, clone)
     4  per-bead UMI and read distributions
     5  clone size distribution
     6  DIFFUSION KERNEL, stratified by clone size band, with per-clone geometry
     7  clone compactness vs random labeling
     8  chain composition and per-bead pairing
     9  spatial layout of representative clones across all size bands
    10  singleton CDR3 adjacency: are one-bead clones sequencing error?

WHY BEAT 6 IS STRATIFIED
    Largest clones run 663 / 1,113 / 974 beads, up to 3.9% of a puck surface.
    Pooling those with 3-bead clones lets one giant clone dominate, since it
    contributes far more distance pairs than a thousand small ones. Bands are
    reported separately so different mechanisms stay visible.

WHY PER-CLONE GEOMETRY RATHER THAN PER-CLONE FITS
    A 3-bead clone yields two distance observations. A decay constant fitted to
    two points carries no information. Per-clone median/max distance, radius of
    gyration, and fraction of UMIs within 50 um are stable at n=2 and reveal
    whether a band hides two populations.

CONVENTIONS ESTABLISHED FROM THE DATA
    - Group clones on `cloneId` ONLY. vGene, cGene and aaSeqCDR3 all vary within
      one cloneId from low-depth assignment and sequencing error.
    - `n_reads` is reads for the (bc, umi) pair; `n_reads_clone` is the subset
      supporting this clone. The ratio is per-UMI purity, not a clone total.
    - UMIs are NOT error corrected upstream.
    - `bc` is the CORRECTED barcode, no suffix. Join key: f"{bc}-1_Puck_211214_{p}".
    - obs x_coord/y_coord are byte-identical to obsm['spatial'].

USAGE
    python Step08a_TCR_Phase0_Characterization.py
    python Step08a_TCR_Phase0_Characterization.py --validate-only
    python Step08a_TCR_Phase0_Characterization.py --min-beads-kernel 3

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
================================================================================
"""

import argparse
import os
import sys
from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

try:
    import anndata as ad
except ImportError:
    sys.exit("anndata not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTPUTS = os.path.join(PROOT, "data/outputs")
OUT = os.path.join(OUTPUTS, "09_tcr_phase0")
FIGDIR = os.path.join(OUT, "figures")

ADATA = os.path.join(OUTPUTS, "04_annotation/all_pucks_annotated_unified.h5ad")
NEOANTIGEN_BEADS = os.path.join(OUTPUTS, "07_neoantigen/neoantigens_per_bead.tsv")
MUTATION_BEADS = os.path.join(OUTPUTS, "05_mutations/SComatic/SingleCell/mutations_per_bead.tsv")
TCR_DIR = os.path.join(PROOT, "data/inputs/tcr/processed")

PUCKS = ["29", "37", "40"]
TCR_CSV = {p: f"B59_{p}_hTCR_tcr.csv" for p in PUCKS}

ANNOT_COL = "unified_annotation"
BASAL_T = "T_cell"
EXPECTED_NEOANTIGEN_BEADS = 110
COORD_TOL_UM = 0.5

# figure house style
FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
C_CORAL = "#ed6a5a"
C_MUSTARD = "#F6D155"
C_GRAY = "#d3d3d3"
C_DARK = "#4a4a4a"
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# kernel binning, microns. Bead pitch ~10 um, puck ~5 mm.
# BUG FIX (v2): the first bin was 0-10 um, but bead pitch is ~10 um, so a
# nearest neighbour sits at ~10.0 and np.histogram assigns it to the NEXT bin.
# Shell one was empty by construction and every near-field rate printed 0.00000.
# The first shell now spans 0-15 um so the first ring of neighbours lands in it.
DIST_BINS = np.array([0, 15, 30, 50, 75, 100, 150, 200, 300, 500, 1000, 2500, 5000],
                     dtype=float)
# clone size bands for BEAT 6 and BEAT 9
BANDS = [("3-5", 3, 5), ("6-20", 6, 20), ("21-100", 21, 100), ("101+", 101, 10**9)]
NEAR_FIELD_UM = 50.0          # "within 50 um of source" metric
ABUNDANT_MIN_BEADS = 5        # BEAT 10: what counts as a parent clone
N_PERM_COMPACT = 200
RNG = np.random.default_rng(20260927)


def log(msg=""):
    print(msg, flush=True)


def banner(msg):
    log()
    log("=" * 72)
    log(msg)
    log("=" * 72)


def savefig(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"), dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    log(f"    figure: {name}.pdf / .png")


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


def bead_id(bc, puck):
    return f"{bc}-1_Puck_211214_{puck}"


def hamming_le1(a, b):
    """True if a and b are equal length and differ by exactly one substitution."""
    if len(a) != len(b):
        return False
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > 1:
                return False
    return d == 1


def band_of(n):
    for name, lo, hi in BANDS:
        if lo <= n <= hi:
            return name
    return None


# ==============================================================================
# PREFLIGHT
# ==============================================================================
def preflight(adata):
    banner("PREFLIGHT  Step05c CB restore validation  (HARD GATE)")

    obs = set(adata.obs_names)
    log(f"  annotated beads: {len(obs):,}")

    if not os.path.exists(NEOANTIGEN_BEADS):
        log(f"  FAIL  missing {NEOANTIGEN_BEADS}")
        return False

    neo = pd.read_csv(NEOANTIGEN_BEADS, sep="\t")
    log(f"  neoantigens_per_bead.tsv: {len(neo):,} rows")

    bcol = next((c for c in ("CB", "bead", "bead_id", "cell_barcode", "barcode")
                 if c in neo.columns), None)
    if bcol is None:
        log(f"  FAIL  no recognizable bead column in {list(neo.columns)}")
        return False
    log(f"  bead column: {bcol}")

    beads = neo[bcol].astype(str).unique()
    log(f"  unique beads: {len(beads):,}  (documented {EXPECTED_NEOANTIGEN_BEADS})")
    if len(beads) != EXPECTED_NEOANTIGEN_BEADS:
        log(f"  WARN  count differs from the documented {EXPECTED_NEOANTIGEN_BEADS}")

    ok = True

    # 1. membership
    missing = [b for b in beads if b not in obs]
    log(f"  membership: {len(beads)-len(missing):,} / {len(beads):,} in adata.obs_names")
    if missing:
        ok = False
        log(f"  FAIL  {len(missing)} bead IDs are not in the annotated set")
        for b in missing[:10]:
            log(f"          {b}")

    # 2. form
    bare = sum(1 for b in beads if "-" not in b)
    withn = sum(1 for b in beads if "N" in b.split("-")[0])
    log(f"  form: {bare} without '-1' suffix (expect 0), {withn} containing N "
        "(expected, corrected barcodes carry N)")
    if bare:
        ok = False
        log("  FAIL  truncated barcodes present, restore_cb() did not run")

    # 3. puck tag
    pk = defaultdict(int)
    for b in beads:
        pk[b.split("_", 1)[-1] if "_" in b else "NO_PUCK_TAG"] += 1
    log("  puck distribution: " + ", ".join(f"{k}:{v}" for k, v in sorted(pk.items())))
    if "NO_PUCK_TAG" in pk:
        ok = False
        log("  FAIL  bead IDs without a puck tag, restore is incomplete")

    # 4. coordinates. A wrong-but-valid restore lands on a different real bead,
    #    which membership passes and this catches.
    if {"x_coord", "y_coord"}.issubset(neo.columns) and "x_coord" in adata.obs.columns:
        present = [b for b in beads if b in obs]
        sub = neo[neo[bcol].isin(present)].drop_duplicates(bcol).set_index(bcol)
        ref = adata.obs.loc[sub.index, ["x_coord", "y_coord"]]
        dx = (sub["x_coord"].values - ref["x_coord"].values)
        dy = (sub["y_coord"].values - ref["y_coord"].values)
        bad = int(((np.abs(dx) > COORD_TOL_UM) | (np.abs(dy) > COORD_TOL_UM)).sum())
        log(f"  coordinates: max |dx| {np.abs(dx).max():.4f}, "
            f"max |dy| {np.abs(dy).max():.4f}, {bad} beads off by > {COORD_TOL_UM} um")
        if bad:
            ok = False
            log("  FAIL  restored barcodes map to the wrong physical bead")
    else:
        log("  WARN  coordinate columns absent, skipping the wrong-bead check")

    # 5. cross-check upstream
    if os.path.exists(MUTATION_BEADS):
        mut = pd.read_csv(MUTATION_BEADS, sep="\t")
        mcol = next((c for c in ("CB", "bead", "bead_id", "cell_barcode", "barcode")
                     if c in mut.columns), None)
        if mcol:
            mb = set(mut[mcol].astype(str))
            orphan = len(set(beads) - mb)
            log(f"  mutations_per_bead.tsv: {len(mb):,} beads, "
                f"{orphan} neoantigen beads absent from it")

    # 6. pooled-calling artifact, reported not gated
    if {"gene", "hgvs_p"}.issubset(neo.columns):
        per = neo.assign(_p=neo[bcol].map(lambda b: b.split("_", 1)[-1])) \
                 .groupby(["gene", "hgvs_p"])["_p"].nunique().value_counts().sort_index()
        log("  mutations by puck count: " + ", ".join(f"{k} puck(s): {v}" for k, v in per.items()))
        log("    Variants were called on the POOLED BAM, so cross-puck appearance")
        log("    is not independent recurrence across patients.")

    log(f"\n  PREFLIGHT {'PASS' if ok else 'FAIL'}")
    return ok


# ==============================================================================
# LOAD
# ==============================================================================
def load_tcr(puck):
    df = pd.read_csv(os.path.join(TCR_DIR, TCR_CSV[puck]), sep=None, engine="python")
    df.columns = [c.strip() for c in df.columns]
    df["puck"] = puck
    df["bead"] = df["bc"].astype(str).map(lambda b: bead_id(b, puck))
    df["umi_purity"] = np.where(df["n_reads"] > 0, df["n_reads_clone"] / df["n_reads"], np.nan)
    return df


def puck_coords(adata, puck):
    """Coordinates for every annotated bead in one puck."""
    m = adata.obs["puck_id"].astype(str).values == f"Puck_211214_{puck}"
    idx = adata.obs_names[m]
    xy = adata.obs.loc[idx, ["x_coord", "y_coord"]].values.astype(float)
    return idx, xy, {b: i for i, b in enumerate(idx)}


# ==============================================================================
# BEATS 1-5, 7, 8
# ==============================================================================
def beat1_join(tcr, adata, rows):
    banner("BEAT 1  join validation")
    obs = set(adata.obs_names)
    for p in PUCKS:
        d = tcr[tcr.puck == p]
        beads = d["bead"].unique()
        hit = sum(1 for b in beads if b in obs)
        log(f"  puck {p}: {len(d):,} rows, {len(beads):,} beads, "
            f"{hit:,} annotated ({100*hit/len(beads):.1f}%)")
        rows.append({"beat": "join", "puck": p, "rows": len(d), "beads": len(beads),
                     "beads_in_adata": hit, "pct_in_adata": round(100*hit/len(beads), 2)})
    log("\n  Beads absent from the annotated set were filtered at QC, not lost in")
    log("  the join. A high rate here is expected.")


def beat2_celltype(tcr, adata, rows):
    """
    v2 REWRITE. The v1 test compared the fraction of T_cell beads carrying TCR
    against the puck background and returned ~1.00x in all three pucks. That
    test had no depth control, which is a real flaw: a bead with more total RNA
    captures more of everything, and T cells are small and RNA-poor relative to
    epithelial cells. A genuine signal could cancel against a depth deficit and
    land at exactly 1.0.

    This version models TCR UMIs per bead across ALL annotated beads in the
    puck (zeros included), against three measures of T-cell content, with
    sequencing depth as a covariate:

      PRIMARY    marker_score_T_cell   continuous, the unthresholded precursor
                 of unified_annotation. Same evidence, same gene sets, strictly
                 more sensitive than the argmax.
      SECONDARY  c2l CD8-positive      absolute abundance from cell2location.
                 CD8 specifically, because the MHCflurry panel is class I and
                 CD8 is the compartment that would recognize these neoantigens.
      ANCHOR     unified_annotation    the categorical call, retained so the
                 continuous result can be reported as a more sensitive version
                 of a test already run rather than a different test.

    Reported two ways: a model-free depth-stratified table, and a GLM. The
    stratified table is the one to trust if they disagree, since it assumes
    nothing about functional form.
    """
    banner("BEAT 2  T-cell association of TCR signal, depth controlled")

    C_CD8 = "c2l_CD8-positive, alpha-beta T cell"
    C_CD4 = "c2l_CD4-positive, alpha-beta T cell"
    MARK = "marker_score_T_cell"

    have = {"marker": MARK in adata.obs.columns,
            "cd8": C_CD8 in adata.obs.columns,
            "cd4": C_CD4 in adata.obs.columns,
            "depth": "total_counts" in adata.obs.columns}
    log("  available measures: " + ", ".join(f"{k}={v}" for k, v in have.items()))
    if not have["depth"]:
        log("  FATAL for this beat: total_counts absent, cannot control depth")
        return

    try:
        import statsmodels.api as sm
        HAVE_SM = True
    except ImportError:
        HAVE_SM = False
        log("  statsmodels not installed, GLM skipped; stratified table still runs")
        log("    conda install -c conda-forge statsmodels")

    ann = adata.obs[ANNOT_COL].astype(str)
    pid = adata.obs["puck_id"].astype(str)

    for p in PUCKS:
        log(f"\n  --- puck {p} ---")
        m = pid.values == f"Puck_211214_{p}"
        idx = adata.obs_names[m]
        df = pd.DataFrame(index=idx)
        df["depth"] = adata.obs.loc[idx, "total_counts"].values.astype(float)
        df["annot"] = ann.loc[idx].values
        df["is_T"] = (df["annot"] == BASAL_T).astype(int)
        if have["marker"]:
            df["marker_T"] = adata.obs.loc[idx, MARK].values.astype(float)
        if have["cd8"]:
            df["cd8"] = adata.obs.loc[idx, C_CD8].values.astype(float)
        if have["cd4"]:
            df["cd4"] = adata.obs.loc[idx, C_CD4].values.astype(float)

        # outcome: TCR UMIs per bead, zeros included
        t = tcr[tcr.puck == p].groupby("bead")["umi"].nunique()
        df["tcr_umi"] = t.reindex(df.index).fillna(0).values.astype(float)
        df["tcr_pos"] = (df["tcr_umi"] > 0).astype(int)

        # the diagnostic that motivated this rewrite
        dep_T = df.loc[df.is_T == 1, "depth"]
        dep_O = df.loc[df.is_T == 0, "depth"]
        log(f"    median total_counts: T_cell {dep_T.median():,.0f} (n={len(dep_T):,}) | "
            f"other {dep_O.median():,.0f} (n={len(dep_O):,}) | "
            f"ratio {dep_T.median()/max(1, dep_O.median()):.3f}")
        log("    If T_cell beads are shallower, equal TCR capture is a POSITIVE")
        log("    signal that the uncontrolled v1 test would have hidden.")

        # ---- model-free: stratify by depth decile ----
        df["dep_dec"] = pd.qcut(df["depth"], 10, labels=False, duplicates="drop")
        log(f"\n    depth-stratified TCR UMIs per bead (model free)")
        log(f"      {'decile':>7} {'median depth':>13} {'n':>7} "
            f"{'T_cell':>9} {'other':>9} {'ratio':>7}")
        ratios = []
        for dd, g in df.groupby("dep_dec"):
            a = g.loc[g.is_T == 1, "tcr_umi"]
            b = g.loc[g.is_T == 0, "tcr_umi"]
            if not len(a) or not len(b):
                continue
            r = a.mean() / b.mean() if b.mean() else np.nan
            ratios.append(r)
            log(f"      {int(dd):>7} {g.depth.median():>13,.0f} {len(g):>7,} "
                f"{a.mean():>9.3f} {b.mean():>9.3f} {r:>7.3f}")
        if ratios:
            log(f"      pooled across deciles: mean ratio {np.nanmean(ratios):.3f}")
            log("      (1.0 = TCR signal is independent of T-cell identity at")
            log("       matched depth; >1.0 = tracks T cells)")

        # same, stratified on the continuous primary instead of the argmax
        if have["marker"]:
            df["mark_q"] = pd.qcut(df["marker_T"], 4, labels=False, duplicates="drop")
            log(f"\n    marker_score_T_cell quartile x depth decile, "
                f"mean TCR UMIs per bead")
            log(f"      {'decile':>7} " + " ".join(f"{'Q'+str(q):>8}" for q in range(4))
                + f" {'Q4/Q1':>7}")
            q41 = []
            for dd, g in df.groupby("dep_dec"):
                means = [g.loc[g.mark_q == q, "tcr_umi"].mean() for q in range(4)]
                if any(mm != mm for mm in means):
                    continue
                r = means[3] / means[0] if means[0] else np.nan
                q41.append(r)
                log(f"      {int(dd):>7} " + " ".join(f"{mm:>8.3f}" for mm in means)
                    + f" {r:>7.3f}")
            if q41:
                log(f"      pooled Q4/Q1 across deciles: {np.nanmean(q41):.3f}")

        # ---- GLM ----
        glm_out = {}
        if HAVE_SM:
            log(f"\n    Poisson GLM, outcome TCR UMIs, offset log(total_counts)")
            off = np.log(np.clip(df["depth"].values, 1, None))
            specs = []
            if have["marker"]:
                specs.append(("PRIMARY   marker_score_T_cell", ["marker_T"]))
            if have["cd8"]:
                specs.append(("SECONDARY c2l CD8", ["cd8"]))
            if have["cd8"] and have["cd4"]:
                specs.append(("          c2l CD8 + CD4", ["cd8", "cd4"]))
            specs.append(("ANCHOR    unified_annotation == T_cell", ["is_T"]))
            for label, cols in specs:
                X = df[cols].values.astype(float)
                X = sm.add_constant(X, has_constant="add")
                try:
                    res = sm.GLM(df["tcr_umi"].values, X,
                                 family=sm.families.Poisson(), offset=off).fit()
                    parts = []
                    for j, c in enumerate(cols, start=1):
                        parts.append(f"{c} beta={res.params[j]:+.4f} "
                                     f"(exp {np.exp(res.params[j]):.3f}) "
                                     f"p={res.pvalues[j]:.3g}")
                        glm_out[f"{c}_beta"] = round(float(res.params[j]), 5)
                        glm_out[f"{c}_p"] = float(res.pvalues[j])
                    log(f"      {label}: " + "; ".join(parts))
                except Exception as e:
                    log(f"      {label}: FAILED ({e})")
            log("      exp(beta) is the multiplicative change in TCR UMIs per")
            log("      unit of the measure, at matched sequencing depth.")

        row = {"beat": "celltype_depth_controlled", "puck": p,
               "median_depth_Tcell": float(dep_T.median()),
               "median_depth_other": float(dep_O.median()),
               "depth_ratio": round(float(dep_T.median()/max(1, dep_O.median())), 4),
               "stratified_ratio_Tcell_vs_other": round(float(np.nanmean(ratios)), 4) if ratios else None}
        row.update(glm_out)
        rows.append(row)

        # keep the v1 statistic so the two are directly comparable
        sub = ann.loc[[b for b in tcr.loc[tcr.puck == p, "bead"].unique() if b in ann.index]]
        bg = ann[m]
        v1 = ((sub == BASAL_T).mean()) / max(1e-9, (bg == BASAL_T).mean())
        log(f"\n    v1 uncontrolled enrichment, for comparison: {v1:.3f}x")

    log("\n  Read the depth-stratified table first. If the within-decile ratio")
    log("  is ~1.0 and the Q4/Q1 marker-score ratio is ~1.0, TCR signal has no")
    log("  relationship to T-cell content even after controlling depth, and the")
    log("  processed table cannot support spatial inference.")


def beat3_umi_error(tcr, rows):
    banner("BEAT 3  UMI error inflation")
    log("  UMIs within one substitution of another on the same bead and clone.")
    for p in PUCKS:
        d = tcr[tcr.puck == p]
        adjacent = groups = 0
        for _, g in d.groupby(["bead", "cloneId"], sort=False):
            u = g["umi"].astype(str).unique()
            if len(u) < 2:
                continue
            groups += 1
            if len(u) > 60:
                continue
            seen = set()
            for a, b in combinations(u, 2):
                if hamming_le1(a, b):
                    seen.add(a); seen.add(b)
            adjacent += max(0, len(seen) - 1) if seen else 0
        pct = 100*adjacent/len(d) if len(d) else 0
        log(f"  puck {p}: {len(d):,} rows, {groups:,} multi-UMI groups, "
            f"{adjacent:,} redundant ({pct:.2f}%)")
        rows.append({"beat": "umi_error", "puck": p, "rows": len(d),
                     "multi_umi_groups": groups, "redundant_umis": adjacent,
                     "pct_inflation": round(pct, 3)})
    log("\n  This is the inflation directional adjacency collapse would remove.")


def beat4_per_bead(tcr, adata, rows):
    banner("BEAT 4  per-bead UMI distributions")
    ann = adata.obs[ANNOT_COL].astype(str)
    fig, axes = plt.subplots(1, 3, figsize=(30, 9))
    for i, p in enumerate(PUCKS):
        d = tcr[tcr.puck == p]
        per = d.groupby("bead").agg(n_umi=("umi", "nunique"), n_reads=("n_reads", "sum"),
                                    n_clones=("cloneId", "nunique"))
        per["annot"] = [ann.get(b, "not_in_adata") for b in per.index]
        q = per["n_umi"].quantile([0.5, 0.9, 0.99, 1.0])
        t = per.loc[per.annot == BASAL_T, "n_umi"]
        o = per.loc[(per.annot != BASAL_T) & (per.annot != "not_in_adata"), "n_umi"]
        log(f"  puck {p}: {len(per):,} beads, median {q[0.5]:.0f}, p90 {q[0.9]:.0f}, "
            f"p99 {q[0.99]:.0f}, max {q[1.0]:.0f}")
        log(f"    T_cell beads n={len(t):,} median {t.median() if len(t) else 'n/a'}; "
            f"other n={len(o):,} median {o.median() if len(o) else 'n/a'}")
        log(f"    single-UMI beads: {(per.n_umi==1).sum():,} ({100*(per.n_umi==1).mean():.1f}%)")
        rows.append({"beat": "per_bead", "puck": p, "beads": len(per),
                     "median_umi": float(q[0.5]), "p99_umi": float(q[0.99]),
                     "max_umi": float(q[1.0]),
                     "pct_single_umi": round(100*(per.n_umi==1).mean(), 2),
                     "median_umi_Tcell": float(t.median()) if len(t) else None,
                     "median_umi_other": float(o.median()) if len(o) else None})
        ax = axes[i]
        bins = np.logspace(0, np.log10(max(2, per.n_umi.max())), 40)
        if len(o):
            ax.hist(o, bins=bins, color=C_GRAY, alpha=0.9, label="other / unannotated")
        if len(t):
            ax.hist(t, bins=bins, color=C_CORAL, alpha=0.9, label="T_cell")
        ax.set_xscale("log"); ax.set_yscale("log")
        style(ax, title=f"Puck {p}", xlabel="UMIs per bead",
              ylabel="Beads" if i == 0 else None)
        if i == 0:
            ax.legend(fontsize=FS_ANNOT, frameon=False)
    fig.tight_layout(); savefig(fig, "Fig_Phase0_A_umi_per_bead")


def beat5_clone_size(tcr, rows):
    banner("BEAT 5  clone size distribution")
    for p in PUCKS:
        d = tcr[tcr.puck == p]
        cl = d.groupby("cloneId").agg(n_beads=("bead", "nunique"), n_umi=("umi", "nunique"))
        tot_u = int(cl.n_umi.sum())
        log(f"  puck {p}: {len(cl):,} clones, {tot_u:,} UMIs, max {int(cl.n_beads.max()):,} beads")
        for name, lo, hi in BANDS:
            s = cl[(cl.n_beads >= lo) & (cl.n_beads <= hi)]
            log(f"    band {name:>6}: {len(s):>6,} clones, {int(s.n_umi.sum()):>7,} UMIs "
                f"({100*s.n_umi.sum()/tot_u:5.2f}%)")
        sing = cl[cl.n_beads == 1]
        log(f"    singletons:  {len(sing):>6,} clones, {int(sing.n_umi.sum()):>7,} UMIs "
            f"({100*sing.n_umi.sum()/tot_u:5.2f}%), {sing.n_umi.mean():.2f} UMI each")
        rows.append({"beat": "clone_size", "puck": p, "clones": len(cl),
                     "singletons": len(sing), "singleton_umis": int(sing.n_umi.sum()),
                     "max_beads": int(cl.n_beads.max())})
    log("\n  Beads per clone is NOT clone size in cells until BEAT 6 settles the")
    log("  diffusion question. One T cell bleeding onto neighbours produces a")
    log("  multi-bead clone with no expansion at all.")


# ==============================================================================
# BEAT 6  stratified diffusion kernel + per-clone geometry
# ==============================================================================
def beat6_diffusion(tcr, adata, rows, min_beads):
    banner("BEAT 6  diffusion kernel, stratified by clone size band")
    log("  For each clone the highest-read bead is the putative source. Clone")
    log("  signal is measured against distance from it, normalized by the beads")
    log("  actually available at that distance. Bands are kept separate so one")
    log("  1,000-bead clone cannot dominate the pooled estimate.")

    centers = 0.5 * (DIST_BINS[:-1] + DIST_BINS[1:])
    per_clone_all = []
    band_curves = {}

    for p in PUCKS:
        idx, xy, pos = puck_coords(adata, p)
        d = tcr[(tcr.puck == p) & (tcr.bead.isin(set(idx)))]

        sig = {b[0]: np.zeros(len(centers)) for b in BANDS}
        avail = {b[0]: np.zeros(len(centers)) for b in BANDS}
        nclone = {b[0]: 0 for b in BANDS}

        for cid, g in d.groupby("cloneId"):
            per_bead = g.groupby("bead").agg(reads=("n_reads", "sum"), umi=("umi", "nunique"))
            nb = len(per_bead)
            if nb < min_beads:
                continue
            bn = band_of(nb)
            if bn is None:
                continue

            src = per_bead["reads"].idxmax()
            si = pos[src]
            sx, sy = xy[si]

            # distances from source to every bead in the puck, vectorized
            dall = np.hypot(xy[:, 0] - sx, xy[:, 1] - sy)
            avail[bn] += np.histogram(np.delete(dall, si), bins=DIST_BINS)[0]

            others = per_bead.drop(index=src)
            if not len(others):
                continue
            oi = [pos[b] for b in others.index]
            dist = dall[oi]
            sig[bn] += np.histogram(dist, bins=DIST_BINS)[0]
            nclone[bn] += 1

            u = others["umi"].values
            near = float(u[dist <= NEAR_FIELD_UM].sum() + per_bead.loc[src, "umi"])
            cxy = xy[[si] + oi]
            cen = cxy.mean(axis=0)
            rg = float(np.sqrt(((cxy - cen) ** 2).sum(axis=1).mean()))
            per_clone_all.append({
                "puck": p, "cloneId": cid, "band": bn, "n_beads": nb,
                "source_bead": src, "source_reads": int(per_bead.loc[src, "reads"]),
                "n_umi": int(per_bead["umi"].sum()),
                "median_dist_um": float(np.median(dist)),
                "max_dist_um": float(dist.max()),
                "Rg_um": round(rg, 2),
                "frac_umi_within_50um": round(near / float(per_bead["umi"].sum()), 4)})

        log(f"\n  --- puck {p} ---")
        for name, _, _ in BANDS:
            if not nclone[name]:
                log(f"    band {name:>6}: no clones")
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                rate = np.where(avail[name] > 0, sig[name] / avail[name], np.nan)
            band_curves[(p, name)] = rate

            # v2: report the FIRST POPULATED shell rather than a hardcoded index.
            # v1 hardcoded rate[0] on a 0-10 um bin that was empty by construction,
            # so every near-field value printed 0.00000 and every contrast 0.0x.
            pop = np.where(avail[name] > 0)[0]
            if not len(pop):
                log(f"    band {name:>6}: {nclone[name]:>5,} clones, no populated shells")
                continue
            k = int(pop[0])
            near = float(rate[k])
            far = float(np.nanmean(rate[-4:]))
            contrast = near / far if far and far == far else float("nan")

            # self-normalized: what the rate would be if clone beads were placed
            # uniformly across the puck. 1.0 means indistinguishable from random.
            base = sig[name].sum() / avail[name].sum() if avail[name].sum() else np.nan
            enr_near = near / base if base and base == base else float("nan")

            log(f"    band {name:>6}: {nclone[name]:>5,} clones")
            log(f"      first shell {DIST_BINS[k]:.0f}-{DIST_BINS[k+1]:.0f} um : "
                f"rate {near:.6f}   enrichment vs uniform {enr_near:6.2f}x")
            log(f"      far field  >{DIST_BINS[-4]:.0f} um      : rate {far:.6f}   "
                f"near/far contrast {contrast:,.2f}x")
            rows.append({"beat": "diffusion", "puck": p, "band": name,
                         "clones": nclone[name],
                         "first_shell_lo_um": float(DIST_BINS[k]),
                         "first_shell_hi_um": float(DIST_BINS[k+1]),
                         "rate_near": near,
                         "rate_far": far if far == far else None,
                         "enrichment_vs_uniform": float(enr_near) if enr_near == enr_near else None,
                         "contrast": float(contrast) if contrast == contrast else None})

    pc = pd.DataFrame(per_clone_all)
    if len(pc):
        pc.to_csv(os.path.join(OUT, "per_clone_geometry.tsv"), sep="\t", index=False)
        log(f"\n  wrote per_clone_geometry.tsv ({len(pc):,} clones)")
        log("\n  per-clone spread by band (fraction of clone UMIs within 50 um of source):")
        log(f"    {'puck':>5} {'band':>7} {'n':>6} {'p10':>7} {'p25':>7} {'p50':>7} "
            f"{'p75':>7} {'p90':>7} {'med Rg':>8}")
        for p in PUCKS:
            for name, _, _ in BANDS:
                s = pc[(pc.puck == p) & (pc.band == name)]
                if not len(s):
                    continue
                f = s.frac_umi_within_50um
                log(f"    {p:>5} {name:>7} {len(s):>6,} {f.quantile(.10):>7.3f} "
                    f"{f.quantile(.25):>7.3f} {f.quantile(.50):>7.3f} "
                    f"{f.quantile(.75):>7.3f} {f.quantile(.90):>7.3f} "
                    f"{s.Rg_um.median():>8.1f}")

    # kernel curves
    fig, axes = plt.subplots(1, 3, figsize=(30, 9), sharey=True)
    colors = [C_CORAL, C_MUSTARD, C_DARK, C_GRAY]
    for i, p in enumerate(PUCKS):
        ax = axes[i]
        for (name, _, _), col in zip(BANDS, colors):
            r = band_curves.get((p, name))
            if r is None:
                continue
            ax.plot(centers, r, "o-", color=col, linewidth=3, markersize=10, label=name)
        ax.set_xscale("log"); ax.set_yscale("log")
        style(ax, title=f"Puck {p}", xlabel="Distance from source (µm)",
              ylabel="P(clone bead)" if i == 0 else None)
        if i == 0:
            ax.legend(fontsize=FS_ANNOT, frameon=False, title="beads/clone",
                      title_fontsize=FS_ANNOT)
    fig.tight_layout(); savefig(fig, "Fig_Phase0_B_diffusion_kernel_by_band")

    # per-clone spread
    if len(pc):
        fig, axes = plt.subplots(1, 3, figsize=(30, 9), sharey=True)
        for i, p in enumerate(PUCKS):
            ax = axes[i]
            data, labels = [], []
            for name, _, _ in BANDS:
                s = pc[(pc.puck == p) & (pc.band == name)]
                if len(s):
                    data.append(s.frac_umi_within_50um.values); labels.append(name)
            if data:
                bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.6)
                for patch, col in zip(bp["boxes"], colors):
                    patch.set_facecolor(col); patch.set_alpha(0.8)
                for med in bp["medians"]:
                    med.set_color(C_DARK); med.set_linewidth(3)
            style(ax, title=f"Puck {p}", xlabel="Beads per clone",
                  ylabel="Clone UMIs within 50 µm of source" if i == 0 else None)
        fig.tight_layout(); savefig(fig, "Fig_Phase0_C_per_clone_spread")

    log("\n  A steep drop in the first shells means one cell bleeding onto")
    log("  neighbours, and the shell where it flattens is the correction scale.")
    log("  A flat curve means clone beads are placed independently of any source,")
    log("  pointing at index hopping or a public clone rather than diffusion.")
    log("  A bimodal per-clone distribution inside one band means two populations.")
    return pc


# ==============================================================================
# BEAT 7
# ==============================================================================
def beat7_compactness(tcr, adata, rows):
    banner("BEAT 7  clone compactness vs random labeling")
    log("  Radius of gyration per clone against the Phase 3 null: bead positions")
    log("  fixed, clone labels permuted among TCR-bearing beads.")

    def rg(a):
        c = a.mean(axis=0)
        return float(np.sqrt(((a - c) ** 2).sum(axis=1).mean()))

    for p in PUCKS:
        idx, xy, pos = puck_coords(adata, p)
        d = tcr[(tcr.puck == p) & (tcr.bead.isin(set(idx)))]
        beads = d["bead"].unique()
        bxy = xy[[pos[b] for b in beads]]

        sizes, obs = [], []
        for _, g in d.groupby("cloneId"):
            bs = g["bead"].unique()
            if len(bs) < 2:
                continue
            sizes.append(len(bs))
            obs.append(rg(xy[[pos[b] for b in bs]]))
        if not sizes:
            log(f"  puck {p}: no multi-bead clones"); continue
        sizes, obs = np.array(sizes), np.array(obs)

        null_med = []
        for _ in range(N_PERM_COMPACT):
            vals = [rg(bxy[RNG.choice(len(bxy), size=s, replace=False)]) for s in sizes]
            null_med.append(np.median(vals))
        null_med = np.array(null_med)

        om, nm = float(np.median(obs)), float(np.median(null_med))
        pval = (np.sum(null_med <= om) + 1) / (len(null_med) + 1)
        log(f"  puck {p}: {len(sizes):,} multi-bead clones")
        log(f"    observed median Rg {om:8.1f} um | null {nm:8.1f} um | "
            f"ratio {om/nm:.3f} | p = {pval:.4f}")
        rows.append({"beat": "compactness", "puck": p, "multi_bead_clones": len(sizes),
                     "obs_median_Rg_um": round(om, 1), "null_median_Rg_um": round(nm, 1),
                     "ratio": round(om/nm, 4), "p_empirical": round(pval, 5)})
    log("\n  A ratio well below 1 means clones are focal, which a real expanded")
    log("  clone and a diffusion halo both produce. BEAT 6 separates them.")


# ==============================================================================
# BEAT 8
# ==============================================================================
def beat8_chains(tcr, rows):
    banner("BEAT 8  chain composition and per-bead pairing")
    for p in PUCKS:
        d = tcr[tcr.puck == p]
        ch = d["chain"].value_counts()
        per = d.groupby("bead")["chain"].apply(set)
        both = sum(1 for s in per if {"TRA", "TRB"} <= s)
        ob = sum(1 for s in per if s == {"TRB"})
        oa = sum(1 for s in per if s == {"TRA"})
        pur = d["umi_purity"].dropna()
        log(f"  puck {p}: " + "; ".join(f"{k}:{v:,}" for k, v in ch.items()))
        log(f"    both chains {both:,} ({100*both/len(per):.1f}%) | "
            f"TRB only {ob:,} | TRA only {oa:,}")
        log(f"    UMI purity median {pur.median():.3f}, {100*(pur<0.9).mean():.1f}% below 0.9")
        rows.append({"beat": "chains", "puck": p, "TRA": int(ch.get("TRA", 0)),
                     "TRB": int(ch.get("TRB", 0)), "beads_both_chains": both,
                     "beads_TRB_only": ob, "beads_TRA_only": oa,
                     "pct_beads_paired": round(100*both/len(per), 2),
                     "median_umi_purity": round(float(pur.median()), 4)})
    log("\n  A low paired fraction means requiring both chains would cost most of")
    log("  the data. Whether that is acceptable is a Phase 1 decision.")


# ==============================================================================
# BEAT 9  clone layouts across size bands
# ==============================================================================
def beat9_layouts(tcr, adata, pc):
    banner("BEAT 9  spatial layout of representative clones across size bands")
    if not len(pc):
        log("  no per-clone geometry available, skipping"); return
    log("  One representative clone per band per puck, chosen as the MEDIAN-sized")
    log("  clone in that band so it is typical rather than extreme.")

    fig, axes = plt.subplots(len(PUCKS), len(BANDS), figsize=(13*len(BANDS), 13*len(PUCKS)))
    for r, p in enumerate(PUCKS):
        idx, xy, pos = puck_coords(adata, p)
        d = tcr[(tcr.puck == p) & (tcr.bead.isin(set(idx)))]
        for c, (name, _, _) in enumerate(BANDS):
            ax = axes[r, c] if len(PUCKS) > 1 else axes[c]
            s = pc[(pc.puck == p) & (pc.band == name)]
            ax.scatter(xy[:, 0], xy[:, 1], s=1, c=C_GRAY, alpha=0.35, linewidths=0)
            if not len(s):
                style(ax, title=f"Puck {p} | {name} | none")
                ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
                continue
            med = s.iloc[(s.n_beads - s.n_beads.median()).abs().argsort().iloc[0]]
            g = d[d.cloneId == med.cloneId]
            pb = g.groupby("bead").agg(umi=("umi", "nunique"))
            ci = [pos[b] for b in pb.index]
            ax.scatter(xy[ci, 0], xy[ci, 1], s=30 + 28*pb["umi"].values,
                       c=C_CORAL, alpha=0.85, linewidths=0)
            si = pos[med.source_bead]
            ax.scatter(xy[si, 0], xy[si, 1], s=420, marker="*", c=C_MUSTARD,
                       edgecolors=C_DARK, linewidths=2.5, zorder=5)
            style(ax, title=f"Puck {p} | {name} | {int(med.n_beads)} beads")
            ax.text(0.03, 0.97, f"Rg {med.Rg_um:.0f} µm\n"
                                f"{med.frac_umi_within_50um:.0%} UMI <50 µm",
                    transform=ax.transAxes, fontsize=FS_ANNOT, va="top")
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); savefig(fig, "Fig_Phase0_D_clone_layouts_by_band")
    log("  Star marks the source bead. Marker size scales with UMIs per bead.")


# ==============================================================================
# BEAT 10  singleton CDR3 adjacency
# ==============================================================================
def beat10_singletons(tcr, adata, rows):
    banner("BEAT 10  singleton CDR3 adjacency: error or rare clonotype?")
    log("  One-bead clones carry ~20% of all UMIs but cannot enter the kernel.")
    log("  If a singleton's CDR3 is one substitution from an abundant clone AND")
    log("  its bead sits inside that clone's territory, it is sequencing error")
    log("  split off from a real clone, not a rare T cell.")

    for p in PUCKS:
        idx, xy, pos = puck_coords(adata, p)
        d = tcr[(tcr.puck == p) & (tcr.bead.isin(set(idx)))]
        cl = d.groupby("cloneId").agg(n_beads=("bead", "nunique"))

        sing_ids = set(cl.index[cl.n_beads == 1])
        abun_ids = set(cl.index[cl.n_beads >= ABUNDANT_MIN_BEADS])
        if not sing_ids or not abun_ids:
            log(f"  puck {p}: insufficient clones"); continue

        # representative CDR3 per clone: the most common nucleotide sequence
        rep = (d.groupby(["cloneId", "chain", "nSeqCDR3"]).size()
                 .reset_index(name="n")
                 .sort_values("n", ascending=False)
                 .drop_duplicates("cloneId"))
        rep_map = {r.cloneId: (r.chain, str(r.nSeqCDR3)) for r in rep.itertuples()}

        # bucket by (chain, CDR3 length) so only comparable sequences are compared
        buckets = defaultdict(list)
        for cid in abun_ids:
            if cid in rep_map:
                ch, s = rep_map[cid]
                buckets[(ch, len(s))].append((cid, s))

        abun_beads = {cid: d.loc[d.cloneId == cid, "bead"].unique() for cid in abun_ids}
        sing_bead = d[d.cloneId.isin(sing_ids)].groupby("cloneId")["bead"].first()

        matched, dists, null_dists = 0, [], []
        tested = 0
        for cid in sing_ids:
            if cid not in rep_map:
                continue
            ch, s = rep_map[cid]
            cands = buckets.get((ch, len(s)))
            if not cands:
                continue
            tested += 1
            hit = next((acid for acid, aseq in cands if hamming_le1(s, aseq)), None)
            if hit is None:
                continue
            matched += 1
            sb = sing_bead.get(cid)
            if sb is None or sb not in pos:
                continue
            axy = xy[[pos[b] for b in abun_beads[hit] if b in pos]]
            if not len(axy):
                continue
            tree = cKDTree(axy)
            dists.append(float(tree.query(xy[pos[sb]])[0]))
            # null: a random TCR bead in this puck, same parent clone
            # v2: 20 null draws per matched singleton. v1 drew one, which at
            # n=45-66 matches was too noisy to interpret and made the median
            # ratio flip sign between pucks.
            rbs = d["bead"].sample(20, replace=True,
                                   random_state=int(RNG.integers(1e9))).values
            for rb in rbs:
                if rb in pos:
                    null_dists.append(float(tree.query(xy[pos[rb]])[0]))

        dists, null_dists = np.array(dists), np.array(null_dists)
        log(f"\n  puck {p}: {len(sing_ids):,} singletons, {len(abun_ids):,} abundant "
            f"(>= {ABUNDANT_MIN_BEADS} beads)")
        log(f"    comparable (chain + length bucket exists): {tested:,}")
        log(f"    within 1 substitution of an abundant clone: {matched:,} "
            f"({100*matched/max(1,tested):.1f}% of comparable, "
            f"{100*matched/len(sing_ids):.1f}% of all singletons)")
        if len(dists):
            log(f"    distance to that clone's nearest bead:")
            log(f"      observed  median {np.median(dists):8.1f} um  "
                f"p25 {np.percentile(dists,25):8.1f}  p75 {np.percentile(dists,75):8.1f}")
            if len(null_dists):
                log(f"      null      median {np.median(null_dists):8.1f} um  "
                    f"p25 {np.percentile(null_dists,25):8.1f}  "
                    f"p75 {np.percentile(null_dists,75):8.1f}")
                log(f"      ratio {np.median(dists)/max(1e-9, np.median(null_dists)):.3f}")
            for thr in (10, 25, 50, 100):
                log(f"      within {thr:3d} um of parent: {100*(dists<=thr).mean():5.1f}%"
                    + (f"   (null {100*(null_dists<=thr).mean():5.1f}%)" if len(null_dists) else ""))
        rows.append({"beat": "singletons", "puck": p, "singletons": len(sing_ids),
                     "abundant": len(abun_ids), "comparable": tested,
                     "matched_1nt": matched,
                     "pct_matched_of_comparable": round(100*matched/max(1, tested), 2),
                     "median_dist_to_parent_um": round(float(np.median(dists)), 1) if len(dists) else None,
                     "median_null_dist_um": round(float(np.median(null_dists)), 1) if len(null_dists) else None})
    log("\n  Observed well below null means error-derived singletons are landing")
    log("  inside their parent clone's territory, and much of the singleton UMI")
    log("  mass is error rather than rare clonotypes.")


# ==============================================================================
# MAIN
# ==============================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--validate-only", action="store_true")
    ap.add_argument("--skip-preflight", action="store_true",
                    help="NOT advised; the preflight gates every downstream join")
    ap.add_argument("--min-beads-kernel", type=int, default=3,
                    help="minimum beads for a clone to enter BEAT 6 (default 3; "
                         "2 gives one distance and no shape)")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)

    needed = [ADATA, NEOANTIGEN_BEADS] + [os.path.join(TCR_DIR, TCR_CSV[p]) for p in PUCKS]
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
    for c in (ANNOT_COL, "x_coord", "y_coord", "puck_id"):
        if c not in adata.obs.columns:
            sys.exit(f"missing obs column: {c}")

    if not args.skip_preflight:
        if not preflight(adata):
            log("\n" + "=" * 72)
            log("PREFLIGHT FAILED. Stopping before any TCR join.")
            log("Rerun Step05c -> Step06 -> Step07 and confirm the log reports")
            log("'all restored barcodes match the annotated bead set'.")
            log("=" * 72)
            sys.exit(1)
    else:
        log("\n  *** PREFLIGHT SKIPPED, results are not trustworthy ***")

    log("\nloading TCR tables")
    tcr = pd.concat([load_tcr(p) for p in PUCKS], ignore_index=True)
    log(f"  {len(tcr):,} rows")
    log(f"  mixcr_cloneId distinct values: {tcr['mixcr_cloneId'].nunique()} "
        "(1 means the column is vestigial)")
    log(f"  MIN_BEADS_FOR_KERNEL = {args.min_beads_kernel}")

    rows = []
    beat1_join(tcr, adata, rows)
    beat2_celltype(tcr, adata, rows)
    beat3_umi_error(tcr, rows)
    beat4_per_bead(tcr, adata, rows)
    beat5_clone_size(tcr, rows)
    pc = beat6_diffusion(tcr, adata, rows, args.min_beads_kernel)
    beat7_compactness(tcr, adata, rows)
    beat8_chains(tcr, rows)
    beat9_layouts(tcr, adata, pc)
    beat10_singletons(tcr, adata, rows)

    out = os.path.join(OUT, "phase0_summary.tsv")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    log(f"\nwrote {out}")

    banner("WHAT PHASE 1 NEEDS FROM THIS")
    log("  1. BEAT 2 enrichment: does TCR signal track T cells at all?")
    log("  2. BEAT 6 contrast per band: how far does signal travel, and do the")
    log("     giant clones behave like the small ones?")
    log("  3. BEAT 6 per-clone spread: is any band hiding two populations?")
    log("  4. BEAT 10 observed vs null: is the singleton UMI mass real or error?")
    log("  5. BEAT 8 paired fraction: what does requiring both chains cost?")
    log("  No thresholds are set here. Those follow from these numbers.")


if __name__ == "__main__":
    main()
