#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step05d_Signature_Refitting.py
=========================================================================
Semi-supervised COSMIC signature refitting on the per-bead genotypes from
Step05c, following the network paper's approach: a fixed core signature
set, plus additional HNSCC-associated signatures admitted by scree-plot
elbow detection, fitted per bead by NNLS.

METHOD
  1. 96-trinucleotide context matrix per bead (pyrimidine convention)
  2. Restrict COSMIC to the HNSCC-associated set
  3. Rank non-core candidates by explained variance when fitted alone
  4. Walk core -> core+1 -> core+2 ... recording Frobenius error
  5. Elbow on that curve sets the final signature count
  6. NNLS per bead against the selected set
  7. Reconstruction quality, per-bead weights, AnnData integration, figures

  SIGNATURE SETS (from the project config)
    core       SBS2, SBS13, SBS5
    candidates SBS1, SBS4, SBS7a, SBS7b, SBS16, SBS17a, SBS17b,
               SBS18, SBS29, SBS39, SBS40, SBS44

-------------------------------------------------------------------------
THREE DEPARTURES FROM THE ClusterCatcher signature_analysis.py
-------------------------------------------------------------------------
(a) ALT BASE FROM ALT_TRI, NOT ALT_expected. ClusterCatcher reads
    df['ALT_expected'], which at a multi-allelic site is "A,C" and
    produces a nonsense mutation class. Step05c wrote ALT_TRI using the
    actual Base_observed, so ref = REF_TRI[1] and alt = ALT_TRI[1] are
    self-consistent and multi-allelic-safe.

(b) NO BLANKET astype(str) ON adata.obs. ClusterCatcher coerces every
    non-numeric obs column to string before writing, which flattens every
    categorical in the annotated object including unified_annotation.
    AnnData writes categoricals natively; only object columns are coerced.

(c) SBS40 CONVENTION. COSMIC v3.4 split SBS40 into SBS40a/b/c. The
    requested set names SBS40. resolve_signature_names() detects which
    convention the COSMIC file uses and expands or collapses, reporting
    what it did rather than silently dropping the signature.

-------------------------------------------------------------------------
READ THIS BEFORE INTERPRETING SBS2
-------------------------------------------------------------------------
  The Step05b spectrum is 32.0% T>C in pyrimidine convention, roughly
  double the next class, and A>G / T>C are near-balanced (396 / 381),
  which is what strand-symmetric ADAR A-to-I editing looks like surviving
  SComatic's editing-site filter. That component was not removed.

  SBS16 is T>C dominant, so it will almost certainly rank first among the
  candidates and absorb much of that signal. SBS2 weight expressed as a
  fraction of total is therefore diluted by an artifact, not only by
  biology. The defensible comparison is epithelial versus the immune and
  stromal compartments, which carry the same editing background. That
  contrast is emitted as the per-cell-type companion fit.

  Roughly 430 beads carry any mutation, and most carry one. A 96-context
  NNLS fit on a single mutation is degenerate: all weight lands on the
  signature with the highest probability in that one context. MUT_THRESHOLD
  controls this. The script reports fit quality stratified by per-bead
  mutation count so the degenerate cells are visible rather than pooled
  into a mean.

Env: NETWORK or slide-TCR-seq (scanpy, pandas, numpy, scipy, matplotlib)
Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import glob
import time
import itertools
import collections

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import nnls
from scipy.stats import pearsonr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# =========================================================================
# CONFIGURATION
# =========================================================================

PROOT   = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTDIR  = f"{PROOT}/data/outputs/05_mutations"
SC_DIR  = f"{OUTDIR}/SComatic"
OUT_SC  = f"{SC_DIR}/SingleCell"

MUTATIONS = f"{OUT_SC}/FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv"
CALLABLE  = f"{OUT_SC}/CombinedCallableSites/complete_callable_sites.tsv"
H5AD      = f"{PROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"
META_FILE = f"{OUTDIR}/meta_unified_annotation.tsv"

OUT       = f"{OUTDIR}/06_signatures"
FIGDIR    = f"{OUT}/figures"

# COSMIC SBS reference matrix: rows = 96 contexts ("A[C>A]A"), columns = SBS.
# Confirmed on disk 2026-09-06. This is the SAME file the network paper used,
# which is what makes the signature weights here comparable to its results.
# v3.4 carries SBS40a/b/c rather than SBS40; resolve_signature_names() expands
# the requested SBS40 into all three.
COSMIC_FILE = ("/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/"
               "COSMIC_v3.4_SBS_GRCh38.txt")
COSMIC_SEARCH = [
    "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/*SBS*.txt",
    "/master/jlehle/WORKING/SC/ref/COSMIC/*SBS*.txt",
    "/master/jlehle/WORKING/**/COSMIC*SBS*.txt",
    "/master/jlehle/WORKING/**/*SBS_GRCh38*.txt",
]

# --- Signature sets ------------------------------------------------------
# Core is fixed and always enters the model. SBS2 (C>T at TCW) and SBS13
# (C>G at TCW) are the APOBEC pair; SBS5 is the flat clock-like background
# that absorbs unstructured counts and keeps the other two from soaking up
# noise. These three match the project config.
CORE_SIGNATURES = ["SBS2", "SBS13", "SBS5"]
HNSCC_SIGNATURES = [
    "SBS1", "SBS2", "SBS4", "SBS5", "SBS7a", "SBS7b", "SBS13",
    "SBS16", "SBS17a", "SBS17b", "SBS18", "SBS29", "SBS39",
    "SBS40", "SBS44",
]

USE_SCREE      = True
MAX_SIGNATURES = 15    # caps a pool of 17 after SBS40 -> SBS40a/b/c

# Minimum mutations per bead to enter the fit. A bead with 1 mutation has
# one non-zero context out of 96 and reconstructs perfectly under any
# signature carrying mass there, so its weights are an artifact of which
# context happened to be hit. 3 is the floor.
#
# Being straight about what 3 does and does not buy: it is not tied to the
# 3 core signatures, and it does not make the system determined. A bead
# with 3 mutations has at most 3 non-zero contexts, so once the scree
# admits more signatures the fit stays underdetermined and NNLS returns
# one of many equally good non-negative solutions. Per-bead weights above
# this threshold are more stable than at 1, not stable in absolute terms.
# per_bead_fit_quality.tsv and the per-cell-type pseudobulk companion are
# the honest reads.
MUT_THRESHOLD  = 3
ELBOW_METHOD   = "second_derivative"   # or "l_method"

# --- House figure style --------------------------------------------------
FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300
COL_SBS2 = "#ed6a5a"
COL_ALT  = "#F6D155"
COSMIC_COLORS = {"C>A": "#1EBFF0", "C>G": "#050708", "C>T": "#E62725",
                 "T>A": "#CBCACB", "T>C": "#A1CE63", "T>G": "#EDB6C2"}

# =========================================================================

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
CONTEXTS_96 = [f"{f}[{r}>{a}]{t}"
               for r in ("C", "T")
               for a in (("A", "G", "T") if r == "C" else ("A", "C", "G"))
               for f in "ACGT" for t in "ACGT"]
CONTEXTS_96.sort(key=lambda s: (CLASSES.index(s[2:5]), s[0], s[-1]))


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] [Step05d] {m}", flush=True)


def section(t):
    print("", flush=True)
    log("=" * 66)
    log(t)
    log("=" * 66)


def rc(s):
    return "".join(COMP[b] for b in reversed(s))


def save_fig(fig, name):
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIGDIR}/{name}.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    log(f"  figure: {name}.pdf / .png")


# =========================================================================
# 96-CONTEXT MATRIX
# =========================================================================

def build_matrix(path):
    """Per-bead 96-context counts, pyrimidine convention."""
    section("BUILDING 96-CONTEXT MATRIX")
    df = pd.read_csv(path, sep="\t")
    log(f"  {len(df):,} bead-variant rows, {df['CB'].nunique():,} beads")

    for c in ("CB", "REF_TRI", "ALT_TRI"):
        if c not in df.columns:
            sys.exit(f"FATAL: column '{c}' missing from {path}")

    n0 = len(df)
    df = df[df["REF_TRI"].str.len().eq(3) & df["ALT_TRI"].str.len().eq(3)]
    df = df[~df["REF_TRI"].str.contains("N") & ~df["ALT_TRI"].str.contains("N")]
    log(f"  dropped {n0 - len(df):,} rows with unresolved context")

    # Departure (a): ref/alt taken from the trinucleotides, not ALT_expected.
    ref = df["REF_TRI"].str[1]
    alt = df["ALT_TRI"].str[1]
    keep = ref.ne(alt)
    df, ref, alt = df[keep], ref[keep], alt[keep]

    tri, r2, a2 = [], [], []
    for t, rr, aa in zip(df["REF_TRI"], ref, alt):
        if rr in ("A", "G"):
            tri.append(rc(t)); r2.append(COMP[rr]); a2.append(COMP[aa])
        else:
            tri.append(t); r2.append(rr); a2.append(aa)
    ctx = pd.Series([f"{t[0]}[{r}>{a}]{t[2]}"
                     for t, r, a in zip(tri, r2, a2)], index=df.index)

    ok = ctx.isin(CONTEXTS_96)
    if (~ok).any():
        log(f"  dropped {int((~ok).sum()):,} rows with an invalid class")
    df, ctx = df[ok], ctx[ok]

    mat = (pd.crosstab(ctx, df["CB"])
             .reindex(index=CONTEXTS_96, fill_value=0)
             .astype(int))
    log(f"  matrix: {mat.shape[0]} contexts x {mat.shape[1]:,} beads, "
        f"{int(mat.values.sum()):,} mutations")

    per = mat.sum(axis=0)
    log("  beads by mutation count:")
    for t in (1, 2, 3, 5, 10, 20):
        log(f"    >= {t:>2}: {int((per >= t).sum()):,}")

    cb2ct = dict(zip(df["CB"], df.get("Cell_type", pd.Series(index=df.index))))
    return mat, cb2ct


# =========================================================================
# COSMIC
# =========================================================================

def resolve_signature_names(requested, available):
    """
    Departure (c): reconcile SBS40 against SBS40a/b/c, and SBS17/SBS7
    against their lettered variants, in whichever direction the file uses.
    """
    out, notes = [], []
    for sig in requested:
        if sig in available:
            out.append(sig)
            continue
        variants = sorted(s for s in available
                          if s.startswith(sig) and s[len(sig):].isalpha())
        if variants:
            out.extend(variants)
            notes.append(f"{sig} -> {', '.join(variants)}")
            continue
        if sig[-1].isalpha() and sig[:-1] in available:
            out.append(sig[:-1])
            notes.append(f"{sig} -> {sig[:-1]}")
            continue
        notes.append(f"{sig} -> NOT FOUND, dropped")
    for n in notes:
        log(f"    {n}")
    return list(dict.fromkeys(out))


def load_cosmic():
    section("LOADING COSMIC SIGNATURES")
    path = COSMIC_FILE
    if not os.path.exists(path):
        log(f"  {path} not found; searching ...")
        hits = []
        for pat in COSMIC_SEARCH:
            hits.extend(glob.glob(pat, recursive=True))
        hits = sorted(set(hits))
        if not hits:
            sys.exit("FATAL: no COSMIC SBS matrix found. Set COSMIC_FILE.")
        log(f"  candidates:\n    " + "\n    ".join(hits[:10]))
        path = hits[0]
        log(f"  using {path}")

    sep = "," if path.endswith(".csv") else "\t"
    cos = pd.read_csv(path, sep=sep, index_col=0)
    log(f"  {cos.shape[0]} rows x {cos.shape[1]} signatures")

    if not set(cos.index) >= set(CONTEXTS_96):
        # some releases use "A[C>A]A" vs "ACA>A"; try the transpose
        if set(cos.columns) >= set(CONTEXTS_96):
            cos = cos.T
            log("  transposed: contexts were on the columns")
        else:
            missing = [c for c in CONTEXTS_96 if c not in cos.index][:5]
            sys.exit(f"FATAL: COSMIC index is not the 96 contexts. "
                     f"Missing e.g. {missing}. Index sample: "
                     f"{list(cos.index[:5])}")
    cos = cos.loc[CONTEXTS_96]

    log("  resolving requested signature names:")
    hn = resolve_signature_names(HNSCC_SIGNATURES, set(cos.columns))
    core = resolve_signature_names(CORE_SIGNATURES, set(cos.columns))
    core = [c for c in core if c in hn]

    sigs = cos[hn].astype(float)
    sigs = sigs / sigs.sum(axis=0)          # each signature sums to 1
    log(f"  HNSCC set ({len(hn)}): {', '.join(hn)}")
    log(f"  core ({len(core)}): {', '.join(core)}")
    return sigs, core


# =========================================================================
# NNLS
# =========================================================================

def fit_nnls(mat, sigs):
    H = sigs.values
    X = mat.values.astype(float)
    W = np.zeros((H.shape[1], X.shape[1]))
    for i in range(X.shape[1]):
        W[:, i], _ = nnls(H, X[:, i])
    W = pd.DataFrame(W, index=sigs.columns, columns=mat.columns)
    R = pd.DataFrame(H @ W.values, index=mat.index, columns=mat.columns)
    return W, R


def frob(mat, recon):
    return float(np.linalg.norm(mat.values - recon.values, "fro"))


def per_cell_quality(mat, recon):
    X, Xr = mat.values.astype(float), recon.values
    cos, pear = [], []
    for i in range(X.shape[1]):
        a, b = X[:, i], Xr[:, i]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        cos.append(float(a @ b / (na * nb)) if na and nb else np.nan)
        if na and nb and np.ptp(a) and np.ptp(b):
            pear.append(pearsonr(a, b)[0])
        else:
            pear.append(np.nan)
    return np.array(cos), np.array(pear)


# =========================================================================
# SCREE SELECTION
# =========================================================================

def l_method(ns, errs):
    """Best two-line fit; the breakpoint is the elbow."""
    best, cut = np.inf, ns[0]
    for i in range(2, len(ns) - 1):
        e = 0.0
        for seg in (slice(0, i + 1), slice(i, len(ns))):
            x, y = np.array(ns[seg], float), np.array(errs[seg], float)
            if len(x) < 2:
                continue
            e += float(np.sum((y - np.polyval(np.polyfit(x, y, 1), x)) ** 2))
        if e < best:
            best, cut = e, ns[i]
    return cut


def select_signatures(mat, sigs, core):
    section("SIGNATURE SELECTION (scree)")
    cands = [s for s in sigs.columns if s not in core]
    Xn = float(np.linalg.norm(mat.values, "fro"))
    if Xn == 0:
        sys.exit("FATAL: mutation matrix is all zeros")

    log("  ranking candidates by explained variance, fitted alone:")
    scored = []
    for s in cands:
        _, R = fit_nnls(mat, sigs[[s]])
        ev = 1.0 - (frob(mat, R) / Xn) ** 2
        scored.append((s, ev))
    scored.sort(key=lambda x: -x[1])
    for s, ev in scored:
        log(f"    {s:<8} explained variance {ev:.4f}")

    ordered = core + [s for s, _ in scored]
    ns, errs, rows = [], [], []
    for n in range(len(core), min(MAX_SIGNATURES, len(ordered)) + 1):
        sel = ordered[:n]
        _, R = fit_nnls(mat, sigs[sel])
        e = frob(mat, R)
        ns.append(n); errs.append(e)
        rows.append({"n": n, "frobenius_error": e,
                     "relative_error": e / Xn,
                     "signatures": ",".join(sel)})
    scree = pd.DataFrame(rows)

    if len(ns) < 3:
        final_n = ns[-1]
        log(f"  too few points for an elbow; using all {final_n}")
    else:
        d2 = np.gradient(np.gradient(np.array(errs, float)))
        n_sd = ns[int(np.argmax(d2))]
        n_lm = l_method(ns, errs)
        log(f"  elbow, second derivative : {n_sd}")
        log(f"  elbow, L-method          : {n_lm}")
        final_n = n_sd if ELBOW_METHOD == "second_derivative" else n_lm
        if n_sd != n_lm:
            log("  NOTE: the two methods disagree. With few points the second")
            log("        derivative is noisy; L-method is the steadier read.")

    selected = ordered[:final_n]
    log(f"  selected {final_n}: {', '.join(selected)}")
    return selected, scree


# =========================================================================
# FIGURES
# =========================================================================

def plot_spectrum(counts, title, name):
    fig, ax = plt.subplots(figsize=(24, 7))
    cols = [COSMIC_COLORS[c[2:5]] for c in CONTEXTS_96]
    ax.bar(range(96), counts, color=cols, width=0.8)
    ax.set_xlim(-1, 96)
    ax.set_xticks(range(96))
    ax.set_xticklabels([f"{c[0]}{c[3]}{c[-1]}" for c in CONTEXTS_96],
                       rotation=90, fontsize=9, family="monospace")
    ax.set_ylabel("Mutations", fontsize=FS_LABEL)
    ax.set_title(title, fontsize=FS_TITLE)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    for i, cl in enumerate(CLASSES):
        ax.add_patch(plt.Rectangle((i * 16 - 0.4, 1.02), 15.8, 0.04,
                                   transform=ax.get_xaxis_transform(),
                                   color=COSMIC_COLORS[cl], clip_on=False))
        ax.text(i * 16 + 7.5, 1.07, cl, transform=ax.get_xaxis_transform(),
                ha="center", fontsize=FS_ANNOT)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save_fig(fig, name)


def plot_scree(scree, final_n):
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.plot(scree["n"], scree["frobenius_error"], "o-", lw=3, ms=12,
            color=COL_SBS2)
    ax.axvline(final_n, ls="--", lw=2.5, color="#555555")
    ax.text(final_n, ax.get_ylim()[1] * 0.95, f"  n = {final_n}",
            fontsize=FS_ANNOT, va="top")
    ax.set_xlabel("Signatures in model", fontsize=FS_LABEL)
    ax.set_ylabel("Frobenius error", fontsize=FS_LABEL)
    ax.set_title("Signature selection", fontsize=FS_TITLE)
    ax.tick_params(labelsize=FS_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save_fig(fig, "Step05d_scree")


def plot_weights(rel):
    mean = rel.mean(axis=1).sort_values()
    fig, ax = plt.subplots(figsize=(11, max(7, 0.6 * len(mean))))
    cols = [COL_SBS2 if s == "SBS2" else COL_ALT for s in mean.index]
    ax.barh(range(len(mean)), mean.values, color=cols)
    ax.set_yticks(range(len(mean)))
    ax.set_yticklabels(mean.index, fontsize=FS_TICK)
    ax.set_xlabel("Mean relative weight", fontsize=FS_LABEL)
    ax.set_title("Signature activity per bead", fontsize=FS_TITLE)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save_fig(fig, "Step05d_mean_weights")


def plot_quality(cos, per):
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    v = cos[~np.isnan(cos)]
    axes[0].hist(v, bins=40, color=COL_SBS2, alpha=0.85)
    axes[0].set_xlabel("Cosine similarity", fontsize=FS_LABEL)
    axes[0].set_ylabel("Beads", fontsize=FS_LABEL)
    axes[0].set_title("Reconstruction quality", fontsize=FS_TITLE)
    axes[1].scatter(per, cos, s=45, alpha=0.5, color=COL_SBS2)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Mutations per bead", fontsize=FS_LABEL)
    axes[1].set_ylabel("Cosine similarity", fontsize=FS_LABEL)
    axes[1].set_title("Quality vs burden", fontsize=FS_TITLE)
    for a in axes:
        a.tick_params(labelsize=FS_TICK)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    save_fig(fig, "Step05d_reconstruction_quality")


def plot_spatial(adata, rel):
    if "x_coord" not in adata.obs or "SBS2" not in rel.index:
        return
    pucks = list(adata.obs["puck_id"].astype(str).unique())
    sb = rel.loc["SBS2"]
    fig, axes = plt.subplots(1, len(pucks), figsize=(11 * len(pucks), 10))
    axes = np.atleast_1d(axes)
    for ax, p in zip(axes, pucks):
        m = adata.obs["puck_id"].astype(str) == p
        o = adata.obs[m]
        ax.scatter(o["x_coord"], o["y_coord"], s=1, color="#dddddd")
        hit = [b for b in o.index if b in sb.index]
        if hit:
            sub = o.loc[hit]
            sca = ax.scatter(sub["x_coord"], sub["y_coord"],
                             c=sb.loc[hit].values, s=60, cmap="plasma",
                             vmin=0, vmax=1, edgecolors="k", linewidths=0.4)
            plt.colorbar(sca, ax=ax).ax.tick_params(labelsize=FS_TICK)
        ax.set_title(f"{p}\nSBS2 weight ({len(hit)} beads)", fontsize=FS_TITLE)
        ax.set_aspect("equal"); ax.axis("off")
    save_fig(fig, "Step05d_SBS2_spatial")


# =========================================================================
# MAIN
# =========================================================================

def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIGDIR, exist_ok=True)

    section("STEP 05d: SIGNATURE REFITTING")
    log(f"  mutations : {MUTATIONS}")
    log(f"  core      : {', '.join(CORE_SIGNATURES)}")
    log(f"  threshold : {MUT_THRESHOLD}")

    if not os.path.exists(MUTATIONS):
        sys.exit(f"FATAL: not found: {MUTATIONS}. Run Step05c first.")

    mat, cb2ct = build_matrix(MUTATIONS)

    if os.path.exists(CALLABLE):
        try:
            cal = pd.read_csv(CALLABLE, sep="\t")
            good = set(cal.loc[cal["SitesPerCell"] > 0, "CB"])
            drop = [c for c in mat.columns if c not in good]
            if drop:
                log(f"  {len(drop):,} beads absent from callable sites, zeroed")
                mat = mat.drop(columns=drop)
        except Exception as e:
            log(f"  WARN: callable sites unusable ({e}); continuing")

    if MUT_THRESHOLD > 0:
        per = mat.sum(axis=0)
        n0 = mat.shape[1]
        mat = mat.loc[:, per >= MUT_THRESHOLD]
        log(f"  threshold {MUT_THRESHOLD}: {mat.shape[1]:,} of {n0:,} beads kept "
            f"({int(mat.values.sum()):,} mutations)")
    if mat.shape[1] == 0:
        sys.exit("FATAL: no beads left after filtering. Lower MUT_THRESHOLD.")

    # The scree ranking is computed on THIS matrix, so a small surviving set
    # makes the selected signature list itself unstable, not just the weights.
    if mat.shape[1] < 30:
        log("")
        log("  WARNING: fewer than 30 beads survived the threshold.")
        log("  Candidate ranking and the elbow are computed on this matrix, so")
        log("  the SELECTED SIGNATURE LIST is unstable here, not just the")
        log("  per-bead weights. Treat the selection as provisional and lean on")
        log("  the per-cell-type pseudobulk fit, which uses every mutation.")
        log("")
    elif mat.shape[1] < 100:
        log(f"  NOTE: {mat.shape[1]} beads is a thin basis for scree selection.")

    sigs, core = load_cosmic()

    if USE_SCREE:
        selected, scree = select_signatures(mat, sigs, core)
        scree.to_csv(f"{OUT}/Step05d_scree.tsv", sep="\t", index=False)
    else:
        selected, scree = list(sigs.columns), None
        log(f"  using all {len(selected)} HNSCC signatures")
    final = sigs[selected]

    section("FITTING")
    W, R = fit_nnls(mat, final)
    tot = W.sum(axis=0).replace(0, np.nan)
    rel = (W / tot).fillna(0.0)

    cos, pear = per_cell_quality(mat, R)
    per = mat.sum(axis=0).values
    Xn = float(np.linalg.norm(mat.values, "fro"))
    log(f"  beads fitted        : {mat.shape[1]:,}")
    log(f"  relative Frobenius  : {frob(mat, R)/Xn:.4f}")
    log(f"  mean cosine         : {np.nanmean(cos):.4f}")
    log(f"  median cosine       : {np.nanmedian(cos):.4f}")

    log("\n  quality by per-bead mutation count:")
    bins = [(1, 1), (2, 2), (3, 4), (5, 9), (10, 10**9)]
    for lo, hi in bins:
        m = (per >= lo) & (per <= hi)
        if m.sum():
            lbl = f"{lo}" if lo == hi else (f"{lo}+" if hi > 10**8 else f"{lo}-{hi}")
            log(f"    {lbl:>5} mut: {int(m.sum()):>5} beads, "
                f"cosine {np.nanmean(cos[m]):.3f}")
    log("    (beads with 1 mutation are degenerate by construction)")

    log("\n  mean relative weight:")
    for s, v in rel.mean(axis=1).sort_values(ascending=False).items():
        log(f"    {s:<8} {v:.4f}")

    # --- outputs: signatures as ROWS, beads as COLUMNS (needs .T on load) ---
    W.to_csv(f"{OUT}/signature_weights_per_cell.txt", sep="\t",
             float_format="%.6f")
    rel.to_csv(f"{OUT}/signature_weights_per_cell_relative.txt", sep="\t",
               float_format="%.6f")
    mat.to_csv(f"{OUT}/mutation_matrix_96contexts.txt", sep="\t")
    pd.DataFrame({"CB": mat.columns, "n_mutations": per,
                  "cosine": cos, "pearson": pear,
                  "cell_type": [cb2ct.get(c, "NA") for c in mat.columns]}
                 ).to_csv(f"{OUT}/per_bead_fit_quality.tsv",
                          sep="\t", index=False)

    # --- companion: pseudobulk per cell type (the stable reference) ---
    section("COMPANION: PER-CELL-TYPE PSEUDOBULK FIT")
    cts = pd.Series({c: cb2ct.get(c, "NA") for c in mat.columns})
    pb = pd.DataFrame({ct: mat.loc[:, cts[cts == ct].index].sum(axis=1)
                       for ct in sorted(set(cts)) if ct != "NA"})
    if pb.shape[1]:
        Wp, Rp = fit_nnls(pb, final)
        relp = (Wp / Wp.sum(axis=0).replace(0, np.nan)).fillna(0.0)
        relp.to_csv(f"{OUT}/signature_weights_per_celltype.txt", sep="\t",
                    float_format="%.6f")
        log("  relative weights per cell type:")
        log("    " + relp.round(3).to_string().replace("\n", "\n    "))
        log("\n  This is the stable comparison. Per-bead weights on 1 to 2")
        log("  mutations are dominated by which context happened to be hit.")

    # --- figures ---
    section("FIGURES")
    plot_spectrum(mat.sum(axis=1).values,
                  f"All beads (n={mat.shape[1]:,}, "
                  f"{int(mat.values.sum()):,} mutations)",
                  "Step05d_spectrum_all")
    for ct in sorted(set(cts)):
        if ct == "NA":
            continue
        cols = cts[cts == ct].index
        if len(cols) < 20:
            continue
        plot_spectrum(mat[cols].sum(axis=1).values,
                      f"{ct} (n={len(cols):,})",
                      f"Step05d_spectrum_{ct}")
    if scree is not None:
        plot_scree(scree, len(selected))
    plot_weights(rel)
    plot_quality(cos, per)

    # --- AnnData ---
    section("ANNDATA INTEGRATION")
    adata = sc.read_h5ad(H5AD)
    adata.obs["total_mutations"] = (
        mat.sum(axis=0).reindex(adata.obs_names).fillna(0).astype(int).values)
    for s in rel.index:
        adata.obs[s] = (rel.loc[s].reindex(adata.obs_names)
                        .fillna(0.0).astype(float).values)
    adata.obs["has_mutations"] = adata.obs["total_mutations"] > 0
    log(f"  beads with mutations in the object: "
        f"{int(adata.obs['has_mutations'].sum()):,}")

    # Departure (b): coerce object columns only; categoricals write natively.
    for c in adata.obs.columns:
        if adata.obs[c].dtype == object:
            adata.obs[c] = adata.obs[c].astype(str)

    plot_spatial(adata, rel)
    out_h5 = f"{OUT}/all_pucks_annotated_signatures.h5ad"
    adata.write_h5ad(out_h5, compression="gzip")
    log(f"  wrote {out_h5}")

    with open(f"{OUT}/Step05d_summary.txt", "w") as f:
        f.write("Step05d: semi-supervised COSMIC signature refitting\n")
        f.write(f"Completed: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Core: {', '.join(core)}\n")
        f.write(f"Selected ({len(selected)}): {', '.join(selected)}\n")
        f.write(f"Beads fitted: {mat.shape[1]:,}\n")
        f.write(f"Mutations: {int(mat.values.sum()):,}\n")
        f.write(f"Relative Frobenius error: {frob(mat, R)/Xn:.4f}\n")
        f.write(f"Mean cosine: {np.nanmean(cos):.4f}\n\n")
        f.write("Mean relative weight per bead:\n")
        for s, v in rel.mean(axis=1).sort_values(ascending=False).items():
            f.write(f"  {s:<8} {v:.4f}\n")
        f.write("\nCAVEAT: the input spectrum is ~32% T>C, consistent with\n")
        f.write("unannotated ADAR editing. SBS2 weight is diluted by that\n")
        f.write("component. Compare epithelial against immune and stromal\n")
        f.write("compartments rather than reading SBS2 in isolation.\n")

    section(f"STEP 05d COMPLETE in {(time.time()-t0)/60:.1f} min")
    log(f"  weights : {OUT}/signature_weights_per_cell.txt  (needs .T on load)")
    log(f"  object  : {out_h5}")
    log(f"  figures : {FIGDIR}/")


if __name__ == "__main__":
    main()
