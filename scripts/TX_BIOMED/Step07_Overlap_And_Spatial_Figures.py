#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step07_Overlap_And_Spatial_Figures.py   (v2)
=========================================================================
Compare the spatial neoantigen set against the network paper's, ask whether
neoantigen-carrying beads cluster in space, and draw the spatial figures.

-------------------------------------------------------------------------
CHANGES IN v2
-------------------------------------------------------------------------
  - Point sizes doubled in Fig A and raised throughout. Beads were too
    small to read on the tissue sections.
  - Bead metadata (puck, compartment, coordinates) is now joined directly
    from the h5ad. v1 relied on columns Step06 attached, which arrived
    all-NaN, so "bead-pairs per puck" printed nothing.
  - Locus-level overlap builds a location key from chrom+pos when the
    network ranking has no `location` column (v1 said "not comparable").
  - NEW: spatial clonality analysis and three new figures (F, G, H).

-------------------------------------------------------------------------
THE CLONALITY QUESTION AND ITS NULL
-------------------------------------------------------------------------
  If several beads in one small region all carry the SAME neoantigen
  mutation, the parsimonious reading is one expanded tumor clone rather
  than several independent mutational events. That is worth knowing, so
  the script tests it directly.

  THE CHOICE OF NULL IS THE WHOLE ANALYSIS. Neoantigen-carrying beads are
  by construction beads that had enough coverage to call a variant at 3
  ALT reads and 5x depth. Testing them against ALL epithelial beads would
  therefore mostly measure where sequencing depth was good, and would
  return "clustered" whether or not a clone exists.

  Two nulls are computed and both reported:
    vs epithelial         every epithelial bead in that puck. Answers
                          "is variant DETECTION spatially uneven?" Expect
                          yes; that is a coverage readout, not biology.
    vs mutation-carrying   only beads already carrying some somatic call.
                          Answers "given that a bead was callable, are
                          neoantigen carriers closer together than
                          chance?" THIS is the one to quote.

  For a specific mutation carried by k beads, the same logic applies to
  its pairwise distances, tested against random k-subsets of the
  mutation-carrying beads in the same puck.

  HONEST LIMIT ON POWER. Around 12 mutations are carried by 2 or more
  beads and the largest is carried by 3. A permutation test at k=2 has a
  floor on how small p can get, and nothing here survives multiple-testing
  correction across 411 mutations. Treat any hit as a candidate to inspect
  in the tissue, not as an established clone. Both raw and BH-adjusted p
  are reported and labelled.

  ACROSS-PUCK RECURRENCE IS NOT CLONALITY. The three pucks are three
  patients, so the same mutation appearing in two pucks is recurrence,
  not one expanded clone. `same_puck` separates the two cases.

-------------------------------------------------------------------------
OVERLAP LEVELS (unchanged from v1)
-------------------------------------------------------------------------
  locus    same base, same change. Strongest, rarest.
  mutation same amino acid substitution. The network paper's own unit.
  peptide  same presented 8-11mer. What a vaccine encodes.
  gene     same gene, possibly different mutation. Weakest.

  Only the gene level gets a hypergeometric test against the 19,885
  protein-coding symbols; the others have no enumerable universe.

  v1 result: mutation 0, peptide 0, gene 24 observed vs 9.6 expected,
  2.51-fold, p = 3.87e-05.

Env: NETWORK or slide-TCR-seq (scanpy, pandas, numpy, scipy, matplotlib)
Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import time
import collections

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Circle
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# =========================================================================
# CONFIGURATION
# =========================================================================

PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
NEO   = f"{PROOT}/data/outputs/07_neoantigen"
OUT   = f"{PROOT}/data/outputs/08_overlap_figures"
FIG   = f"{OUT}/figures"

SPATIAL_MUT  = f"{NEO}/epithelial_neoantigens_per_mutation.tsv"
SPATIAL_BEAD = f"{NEO}/neoantigens_per_bead.tsv"
BEAD_MAP     = (f"{PROOT}/data/outputs/05_mutations/SComatic/SingleCell/"
                "FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv")
H5AD         = f"{PROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"

NMF   = "/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_7"
NMF_RANKING = f"{NMF}/06_prevalence_ranking/neoantigen_prevalence_ranking_full.tsv"
NMF_GROUPS  = {"SBS2_HIGH": f"{NMF}/03_mhc_binding/SBS2_HIGH_neoantigens.tsv",
               "CNV_HIGH":  f"{NMF}/03_mhc_binding/CNV_HIGH_neoantigens.tsv"}

N_CODING_GENES = 19885
TOP_N          = 20

N_PERM         = 10000    # permutations for the clustering nulls
HOTSPOT_RADIUS = 150.0    # microns; beads within this join one hotspot
RANDOM_SEED    = 0

# --- House figure style. Point sizes doubled from v1. --------------------
FS_TITLE, FS_LABEL, FS_TICK, FS_ANNOT = 34, 30, 28, 28
DPI = 300

PT_CELLTYPE = 6.0      # was 1.5
PT_BG       = 4.8      # was 1.2
PT_MUT      = 165      # was 55
PT_NEO      = 390      # was 130
PT_NEO_BIG  = 630      # dedicated neoantigen figure

CT_COLORS = {
    "epithelial":    "#ed6a5a",
    "myeloid":       "#F6D155",
    "B_cell":        "#5B8FF9",
    "T_cell":        "#61DDAA",
    "fibroblast":    "#9B59B6",
    "endothelial":   "#F6903D",
    "smooth_muscle": "#7262FD",
    "mast":          "#78D3F8",
    "ambiguous":     "#C2C8D5",
}
C_RED     = "#D7263D"
C_TCW     = "#ed6a5a"
C_NONTCW  = "#9AA0A6"
C_SHARED  = "#9B59B6"
C_SPATIAL = "#ed6a5a"
C_NETWORK = "#F6D155"
C_BG      = "#E8E8E8"
C_EPI_BG  = "#F6C7C1"

_report = []
rng = np.random.default_rng(RANDOM_SEED)


def log(m=""):
    print(m, flush=True)
    _report.append(m)


def sep(t=""):
    log("")
    log("=" * 78)
    if t:
        log(f"  {t}")
        log("=" * 78)


def pick(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIG}/{name}.{ext}", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    log(f"  figure: {name}.pdf / .png")


def bh(pvals):
    p = np.asarray(pvals, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    adj = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        adj[i] = prev
    return adj


# =========================================================================
# OVERLAP
# =========================================================================

def load_network():
    if os.path.exists(NMF_RANKING):
        df = pd.read_csv(NMF_RANKING, sep="\t")
        log(f"  ranking: {len(df):,} neoantigen mutations")
        tc = pick(df, "tier")
        if tc:
            for k, v in df[tc].value_counts().items():
                log(f"    {k}: {v:,}")
        return df, tc
    log(f"  {NMF_RANKING} not found; using per-group binder sets")
    parts = []
    for g, p in NMF_GROUPS.items():
        if os.path.exists(p):
            d = pd.read_csv(p, sep="\t")
            d["group"] = g
            parts.append(d)
            log(f"    {g}: {len(d):,} rows")
    if not parts:
        sys.exit("FATAL: no network paper neoantigen files readable")
    return pd.concat(parts, ignore_index=True), None


def location_key(df):
    """`location` if present, else chrom:pos assembled from what exists."""
    c = pick(df, "location")
    if c:
        return df[c].astype(str)
    ch = pick(df, "chrom", "#CHROM", "CHROM")
    po = pick(df, "pos", "POS", "Start")
    if ch and po:
        return df[ch].astype(str) + ":" + df[po].astype(str)
    return None


def overlap(spatial, network, tier_col):
    res = {}
    sg, sh = pick(spatial, "gene"), pick(spatial, "hgvs_p")
    ng, nh = pick(network, "gene"), pick(network, "hgvs_p")
    sp = pick(spatial, "mut_peptide")
    np_ = pick(network, "mut_peptide", "mut_peptide_best", "peptide")

    sl, nl = location_key(spatial), location_key(network)
    if sl is not None and nl is not None:
        hits = set(sl) & set(nl)
        res["locus"] = sorted(hits)
        log(f"  locus (chrom:pos)      : {len(hits)} shared "
            f"(spatial {sl.nunique():,}, network {nl.nunique():,})")
    else:
        res["locus"] = []
        log("  locus                  : not comparable "
            "(no chrom/pos in the network table)")

    a = set(zip(spatial[sg].astype(str), spatial[sh].astype(str)))
    b = set(zip(network[ng].astype(str), network[nh].astype(str)))
    res["mutation"] = sorted(a & b)
    log(f"  mutation (gene,hgvs_p) : {len(a & b)} shared "
        f"(spatial {len(a):,}, network {len(b):,})")
    if (a & b) and tier_col:
        key = network.set_index([ng, nh])[tier_col].to_dict()
        for k, v in collections.Counter(
                key.get(h, "?") for h in (a & b)).most_common():
            log(f"      tier {k}: {v}")

    if sp and np_:
        A = set(spatial[sp].dropna().astype(str))
        B = set(network[np_].dropna().astype(str))
        res["peptide"] = sorted(A & B)
        log(f"  peptide (mut 8-11mer)  : {len(A & B)} shared "
            f"(spatial {len(A):,}, network {len(B):,})")
    else:
        res["peptide"] = []
        log("  peptide                : not comparable")

    A = set(spatial[sg].dropna().astype(str))
    B = set(network[ng].dropna().astype(str))
    hits = A & B
    res["gene"] = sorted(hits)
    k, K, n, N = len(hits), len(B), len(A), N_CODING_GENES
    exp = K * n / N
    p = hypergeom.sf(k - 1, N, K, n) if k else 1.0
    res["gene_stats"] = {"overlap": k, "spatial": n, "network": K,
                         "expected": exp,
                         "fold": k / exp if exp else np.nan, "p": p}
    log(f"  gene                   : {k} shared (spatial {n:,}, network {K:,})")
    log(f"      expected {exp:.1f}, fold {k/exp:.2f}, "
        f"hypergeometric p = {p:.3g}")
    if hits:
        log(f"      {', '.join(sorted(hits))}")
    return res


# =========================================================================
# SPATIAL CLONALITY
# =========================================================================

def mean_nn(xy):
    if len(xy) < 2:
        return np.nan
    d, _ = cKDTree(xy).query(xy, k=2)
    return float(d[:, 1].mean())


def perm_cluster(target_xy, pool_xy, n_perm=N_PERM):
    """One-sided: is the target's mean NN distance smaller than chance?"""
    k = len(target_xy)
    if k < 2 or len(pool_xy) <= k:
        return np.nan, np.nan, np.nan
    obs = mean_nn(target_xy)
    null = np.empty(n_perm)
    for i in range(n_perm):
        idx = rng.choice(len(pool_xy), size=k, replace=False)
        null[i] = mean_nn(pool_xy[idx])
    p = (np.sum(null <= obs) + 1) / (n_perm + 1)
    return obs, float(np.median(null)), float(p)


def analyze_clonality(pb, obs_df, mut_beads):
    rows = []
    for (g, h), grp in pb.groupby(["gene", "hgvs_p"]):
        cbs = [c for c in grp["CB"].unique() if c in obs_df.index]
        if len(cbs) < 2:
            continue
        sub = obs_df.loc[cbs]
        pucks = sorted(sub["puck_id"].unique())
        rec = {"gene": g, "hgvs_p": h, "n_beads": len(cbs),
               "n_pucks": len(pucks), "pucks": ",".join(pucks),
               "same_puck": len(pucks) == 1}
        if len(pucks) == 1:
            xy = sub[["x_coord", "y_coord"]].to_numpy(float)
            d = [float(np.hypot(*(xy[i] - xy[j])))
                 for i in range(len(xy)) for j in range(i + 1, len(xy))]
            rec["min_pair_dist_um"] = float(np.min(d))
            rec["mean_pair_dist_um"] = float(np.mean(d))
            pool_idx = [c for c in mut_beads if c in obs_df.index
                        and obs_df.at[c, "puck_id"] == pucks[0]]
            pool = obs_df.loc[pool_idx, ["x_coord", "y_coord"]].to_numpy(float)
            o, m, p = perm_cluster(xy, pool)
            rec.update({"obs_mean_nn_um": o, "null_median_nn_um": m,
                        "p_raw": p, "pool_size": len(pool)})
        else:
            rec.update({"min_pair_dist_um": np.nan, "mean_pair_dist_um": np.nan,
                        "obs_mean_nn_um": np.nan, "null_median_nn_um": np.nan,
                        "p_raw": np.nan, "pool_size": np.nan})
        rows.append(rec)

    df = pd.DataFrame(rows)
    if len(df):
        m = df["p_raw"].notna()
        df["p_bh"] = np.nan
        if m.any():
            df.loc[m, "p_bh"] = bh(df.loc[m, "p_raw"].to_numpy())
        df = df.sort_values(["same_puck", "n_beads", "min_pair_dist_um"],
                            ascending=[False, False, True])
    return df


def find_hotspots(pb, obs_df, radius=HOTSPOT_RADIUS):
    """Connected components of neoantigen beads within `radius`, per puck."""
    rows, hid = [], 0
    beads = [c for c in pb["CB"].unique() if c in obs_df.index]
    if not beads:
        return pd.DataFrame()
    sub_all = obs_df.loc[beads]
    for puck, sub in sub_all.groupby("puck_id"):
        xy = sub[["x_coord", "y_coord"]].to_numpy(float)
        if len(xy) < 2:
            continue
        pairs = cKDTree(xy).query_pairs(radius, output_type="ndarray")
        if len(pairs) == 0:
            continue
        adj = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                         shape=(len(xy), len(xy)))
        n_c, lab = connected_components(adj, directed=False)
        for c in range(n_c):
            members = list(sub.index[lab == c])
            if len(members) < 2:
                continue
            hid += 1
            mm = pb[pb["CB"].isin(members)]
            muts = sorted(set(zip(mm["gene"], mm["hgvs_p"])))
            shared = [f"{a} {b}" for a, b in muts
                      if mm[(mm["gene"] == a)
                            & (mm["hgvs_p"] == b)]["CB"].nunique() > 1]
            xs = sub.loc[members, "x_coord"].to_numpy(float)
            ys = sub.loc[members, "y_coord"].to_numpy(float)
            rows.append({
                "hotspot_id": hid, "puck_id": puck, "n_beads": len(members),
                "n_distinct_mutations": len(muts),
                "n_genes": mm["gene"].nunique(),
                "shared_mutations_within_hotspot": ";".join(shared) or "none",
                "clonal_candidate": bool(shared),
                "x_center": float(xs.mean()), "y_center": float(ys.mean()),
                "max_extent_um": float(np.hypot(np.ptp(xs), np.ptp(ys))),
                "genes": ",".join(sorted(mm["gene"].unique())),
                "beads": ",".join(members),
            })
    return pd.DataFrame(rows)


# =========================================================================
# FIGURES
# =========================================================================

def fig_celltypes(obs):
    pucks = sorted(obs["puck_id"].unique())
    fig, axes = plt.subplots(1, len(pucks), figsize=(13 * len(pucks), 13))
    axes = np.atleast_1d(axes)
    for ax, p in zip(axes, pucks):
        o = obs[obs["puck_id"] == p]
        for ct, col in CT_COLORS.items():
            s = o[o["unified_annotation"] == ct]
            if len(s):
                ax.scatter(s["x_coord"], s["y_coord"], s=PT_CELLTYPE, c=col,
                           linewidths=0, rasterized=True)
        ax.set_title(f"{p}\n{len(o):,} beads", fontsize=FS_TITLE)
        ax.set_aspect("equal")
        ax.axis("off")
    axes[-1].legend(
        handles=[Patch(facecolor=c, label=t.replace("_", " "))
                 for t, c in CT_COLORS.items()],
        fontsize=FS_ANNOT, loc="center left", bbox_to_anchor=(1.02, 0.5),
        frameon=False, markerscale=3)
    save(fig, "Fig_A_celltypes_per_puck")


def fig_spatial_mutations(obs, mut_beads, neo_beads):
    pucks = sorted(obs["puck_id"].unique())
    fig, axes = plt.subplots(1, len(pucks), figsize=(13 * len(pucks), 13))
    axes = np.atleast_1d(axes)
    for ax, p in zip(axes, pucks):
        o = obs[obs["puck_id"] == p]
        ax.scatter(o["x_coord"], o["y_coord"], s=PT_BG, c=C_BG,
                   linewidths=0, rasterized=True)
        ep = o[o["unified_annotation"] == "epithelial"]
        ax.scatter(ep["x_coord"], ep["y_coord"], s=PT_BG, c=C_EPI_BG,
                   linewidths=0, rasterized=True)
        m = o.loc[[b for b in mut_beads if b in o.index]]
        if len(m):
            ax.scatter(m["x_coord"], m["y_coord"], s=PT_MUT, c=C_NETWORK,
                       edgecolors="#7a6a1f", linewidths=0.7, zorder=3)
        n = o.loc[[b for b in neo_beads if b in o.index]]
        if len(n):
            ax.scatter(n["x_coord"], n["y_coord"], s=PT_NEO, c=C_SPATIAL,
                       edgecolors="black", linewidths=1.0, zorder=4)
        ax.set_title(f"{p}\n{len(m)} mutation beads, {len(n)} neoantigen beads",
                     fontsize=FS_TITLE)
        ax.set_aspect("equal")
        ax.axis("off")
    axes[-1].legend(handles=[
        Patch(facecolor=C_BG, label="all beads"),
        Patch(facecolor=C_EPI_BG, label="epithelial"),
        Patch(facecolor=C_NETWORK, label="carries a somatic call"),
        Patch(facecolor=C_SPATIAL, label="carries a predicted neoantigen")],
        fontsize=FS_ANNOT, loc="center left", bbox_to_anchor=(1.02, 0.5),
        frameon=False, markerscale=2)
    save(fig, "Fig_B_spatial_neoantigen_beads")


def fig_neoantigen_only(obs, neo_beads, hotspots):
    """NEW: where every neoantigen-carrying bead sits. Large red points."""
    pucks = sorted(obs["puck_id"].unique())
    fig, axes = plt.subplots(1, len(pucks), figsize=(14 * len(pucks), 14))
    axes = np.atleast_1d(axes)
    has_hs = len(hotspots) > 0 and "puck_id" in hotspots.columns
    for ax, p in zip(axes, pucks):
        o = obs[obs["puck_id"] == p]
        ax.scatter(o["x_coord"], o["y_coord"], s=PT_BG, c=C_BG,
                   linewidths=0, rasterized=True)
        ep = o[o["unified_annotation"] == "epithelial"]
        ax.scatter(ep["x_coord"], ep["y_coord"], s=PT_BG, c=C_EPI_BG,
                   linewidths=0, rasterized=True)
        n = o.loc[[b for b in neo_beads if b in o.index]]
        if len(n):
            ax.scatter(n["x_coord"], n["y_coord"], s=PT_NEO_BIG, c=C_RED,
                       edgecolors="black", linewidths=1.4, zorder=5)
        nh = 0
        if has_hs:
            hp = hotspots[hotspots["puck_id"] == p]
            nh = len(hp)
            for _, hs in hp.iterrows():
                ax.add_patch(Circle(
                    (hs["x_center"], hs["y_center"]),
                    max(hs["max_extent_um"], HOTSPOT_RADIUS) * 0.9,
                    fill=False, lw=3,
                    ls="-" if hs["clonal_candidate"] else "--",
                    ec="#2b2b2b", zorder=6))
        ax.set_title(f"{p}\n{len(n)} neoantigen beads, {nh} hotspots",
                     fontsize=FS_TITLE)
        ax.set_aspect("equal")
        ax.axis("off")
    axes[-1].legend(handles=[
        Patch(facecolor=C_EPI_BG, label="epithelial"),
        Patch(facecolor=C_RED, label="expresses a predicted neoantigen"),
        Patch(facecolor="white", edgecolor="#2b2b2b",
              label="hotspot (solid = shared mutation)")],
        fontsize=FS_ANNOT, loc="center left", bbox_to_anchor=(1.02, 0.5),
        frameon=False, markerscale=2)
    save(fig, "Fig_F_neoantigen_beads_only")


def fig_clonality(obs, pb, clon):
    """Carriers of one multi-bead mutation joined by a line."""
    if not len(clon):
        log("  no multi-bead neoantigens; skipping Fig G")
        return
    keys = list(zip(clon["gene"], clon["hgvs_p"]))
    cmap = plt.get_cmap("tab20")
    colors = {k: cmap(i % 20) for i, k in enumerate(keys)}
    pucks = sorted(obs["puck_id"].unique())
    fig, axes = plt.subplots(1, len(pucks), figsize=(13 * len(pucks), 13))
    axes = np.atleast_1d(axes)
    for ax, p in zip(axes, pucks):
        o = obs[obs["puck_id"] == p]
        ax.scatter(o["x_coord"], o["y_coord"], s=PT_BG, c=C_BG,
                   linewidths=0, rasterized=True)
        drawn = 0
        for (g, h) in keys:
            cbs = [c for c in pb[(pb["gene"] == g)
                                 & (pb["hgvs_p"] == h)]["CB"].unique()
                   if c in o.index]
            if len(cbs) < 2:
                continue
            s = o.loc[cbs]
            x = s["x_coord"].to_numpy(float)
            y = s["y_coord"].to_numpy(float)
            for i in range(len(x)):
                for j in range(i + 1, len(x)):
                    ax.plot([x[i], x[j]], [y[i], y[j]], lw=2.5,
                            color=colors[(g, h)], alpha=0.85, zorder=4)
            ax.scatter(x, y, s=PT_NEO, color=colors[(g, h)],
                       edgecolors="black", linewidths=1.1, zorder=5,
                       label=f"{g} {h}")
            drawn += 1
        ax.set_title(f"{p}\n{drawn} multi-bead neoantigens", fontsize=FS_TITLE)
        ax.set_aspect("equal")
        ax.axis("off")
        if drawn:
            ax.legend(fontsize=FS_ANNOT - 10, loc="upper right", frameon=False)
    save(fig, "Fig_G_multibead_neoantigen_clonality")


def fig_clustering(stats):
    if not stats:
        log("  no clustering statistics; skipping Fig H")
        return
    fig, ax = plt.subplots(figsize=(16, 10))
    labs = [f"{s['puck'].replace('Puck_211214_','P')}\nvs {s['pool']}"
            for s in stats]
    obs_v = [s["obs"] for s in stats]
    null_v = [s["null_median"] for s in stats]
    x = np.arange(len(labs))
    ax.bar(x - 0.2, null_v, 0.4, color="#C2C8D5", label="null median")
    ax.bar(x + 0.2, obs_v, 0.4, color=C_RED, label="observed")
    for i, s in enumerate(stats):
        ax.text(i + 0.2, s["obs"], f" p={s['p']:.3f}", ha="center",
                va="bottom", fontsize=FS_ANNOT - 10)
    ax.set_xticks(x)
    ax.set_xticklabels(labs, fontsize=FS_TICK - 8)
    ax.set_ylabel("Mean nearest-neighbour\ndistance (um)", fontsize=FS_LABEL)
    ax.set_title("Are neoantigen beads clustered?", fontsize=FS_TITLE)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    ax.legend(fontsize=FS_ANNOT, frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(0.5, -0.26, "lower = more clustered;  quote the "
            "'vs mutation-carrying' comparison",
            transform=ax.transAxes, ha="center", fontsize=FS_ANNOT - 10,
            color="#555555")
    save(fig, "Fig_H_spatial_clustering")


def fig_burden(obs, bead_ct):
    cts = [c for c in CT_COLORS if c in set(obs["unified_annotation"])]
    n_bead = obs["unified_annotation"].value_counts()
    n_mut = collections.Counter(bead_ct.values())
    fig, axes = plt.subplots(1, 2, figsize=(24, 10))
    axes[0].bar(range(len(cts)), [n_bead.get(c, 0) for c in cts],
                color=[CT_COLORS[c] for c in cts])
    axes[0].set_ylabel("Beads", fontsize=FS_LABEL)
    axes[0].set_title("Compartment size", fontsize=FS_TITLE)
    rate = [1000 * n_mut.get(c, 0) / max(n_bead.get(c, 1), 1) for c in cts]
    axes[1].bar(range(len(cts)), rate, color=[CT_COLORS[c] for c in cts])
    axes[1].set_ylabel("Mutation-carrying beads\nper 1,000 beads",
                       fontsize=FS_LABEL)
    axes[1].set_title("Detection rate by compartment", fontsize=FS_TITLE)
    for ax in axes:
        ax.set_xticks(range(len(cts)))
        ax.set_xticklabels([c.replace("_", " ") for c in cts],
                           rotation=45, ha="right", fontsize=FS_TICK)
        ax.tick_params(axis="y", labelsize=FS_TICK)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    save(fig, "Fig_C_burden_by_celltype")


def fig_candidates(mut, shared_genes):
    d = mut.sort_values(["n_beads", "delta_ic50"],
                        ascending=[False, False]).head(TOP_N).iloc[::-1]
    fig, ax = plt.subplots(figsize=(17, max(10, 0.75 * len(d))))
    y = np.arange(len(d))
    ax.barh(y, d["wt_ic50"].clip(upper=5e4), color="#DDDDDD", label="wild-type")
    ax.barh(y, d["mut_ic50"], height=0.55,
            color=[C_TCW if t else C_NONTCW for t in d["is_tcw_ct"]])
    ax.set_xscale("log")
    ax.axvline(500, ls="--", lw=2, color="#555555")
    ax.text(500, len(d), " 500 nM", fontsize=FS_ANNOT, va="top")
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{r['gene']} {r['hgvs_p']} ({int(r['n_beads'])})"
         + (" *" if r["gene"] in shared_genes else "")
         for _, r in d.iterrows()], fontsize=FS_TICK)
    ax.set_xlabel("MHC-I IC50 (nM, log scale)", fontsize=FS_LABEL)
    ax.set_title(f"Top {len(d)} neoantigens by bead support", fontsize=FS_TITLE)
    ax.tick_params(axis="x", labelsize=FS_TICK)
    ax.legend(handles=[Patch(facecolor=C_TCW, label="clean TCW C>T (APOBEC)"),
                       Patch(facecolor=C_NONTCW, label="other substitution"),
                       Patch(facecolor="#DDDDDD", label="wild-type")],
              fontsize=FS_ANNOT, frameon=False, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(1.0, -0.06,
            "(n) = beads with support;  * = gene shared with the network paper",
            transform=ax.transAxes, ha="right", fontsize=FS_ANNOT - 6,
            color="#555555")
    save(fig, "Fig_D_top_candidates")


def fig_overlap(res, n_spatial, n_network):
    fig, axes = plt.subplots(1, 2, figsize=(24, 11))
    ax = axes[0]
    k = len(res.get("mutation", []))
    for x, c, lab in ((-0.55, C_SPATIAL, f"spatial\n{n_spatial:,}"),
                      (0.55, C_NETWORK, f"network paper\n{n_network:,}")):
        ax.add_patch(Circle((x, 0), 1.0, alpha=0.55, color=c, lw=0))
        ax.text(x * 1.75, 0, lab, ha="center", va="center", fontsize=FS_ANNOT)
    ax.text(0, 0, str(k), ha="center", va="center", fontsize=FS_TITLE,
            fontweight="bold")
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Neoantigen mutations\n(gene, hgvs_p)", fontsize=FS_TITLE)

    ax = axes[1]
    levels = ["locus", "mutation", "peptide", "gene"]
    vals = [len(res.get(l, [])) for l in levels]
    ax.bar(range(len(levels)), vals, color=C_SHARED)
    for i, v in enumerate(vals):
        ax.text(i, v, f" {v}", ha="center", va="bottom", fontsize=FS_ANNOT)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels(levels, fontsize=FS_TICK)
    ax.set_ylabel("Shared with network paper", fontsize=FS_LABEL)
    ax.set_title("Overlap by comparison level", fontsize=FS_TITLE)
    ax.tick_params(axis="y", labelsize=FS_TICK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    gs = res.get("gene_stats", {})
    if gs:
        ax.text(0.5, -0.22,
                f"gene level: {gs['overlap']} observed vs {gs['expected']:.1f} "
                f"expected, {gs['fold']:.2f}-fold, p = {gs['p']:.2g}",
                transform=ax.transAxes, ha="center", fontsize=FS_ANNOT - 6,
                color="#555555")
    save(fig, "Fig_E_overlap")


# =========================================================================
# MAIN
# =========================================================================

def main():
    t0 = time.time()
    os.makedirs(FIG, exist_ok=True)

    sep("STEP 07 v2: OVERLAP, CLONALITY AND SPATIAL FIGURES")

    sep("STEP 1: load spatial neoantigens")
    mut = pd.read_csv(SPATIAL_MUT, sep="\t")
    log(f"  {len(mut):,} neoantigen mutations, "
        f"{int((mut['n_beads'] > 0).sum()):,} with bead support")
    log(f"  clean TCW C>T: {int(mut['is_tcw_ct'].sum())}")
    pb = pd.read_csv(SPATIAL_BEAD, sep="\t")
    log(f"  per-bead rows: {len(pb):,}, beads: {pb['CB'].nunique():,}")

    # v2: metadata straight from the h5ad, not Step06's columns.
    ad = sc.read_h5ad(H5AD)
    obs = ad.obs.copy()
    for c in ("x_coord", "y_coord"):
        obs[c] = pd.to_numeric(obs[c], errors="coerce")
    obs["puck_id"] = obs["puck_id"].astype(str)
    obs["unified_annotation"] = obs["unified_annotation"].astype(str)
    for c in ("puck_id", "unified_annotation", "x_coord", "y_coord"):
        pb[c] = pb["CB"].map(obs[c])
    miss = int(pb["puck_id"].isna().sum())
    if miss:
        # This is the single point of failure for every spatial panel. v1 and
        # v2 both dropped the rows and carried on, which made Fig F render an
        # empty tissue section and looked like "no neoantigen-carrying beads"
        # rather than "the barcodes did not match". Fail loudly instead.
        log(f"  {miss} of {len(pb)} bead rows did not join to the h5ad")
        log(f"    bead map CB : {pb['CB'].iloc[0]!r}")
        log(f"    h5ad obs    : {obs.index[0]!r}")
        if miss == len(pb):
            log("")
            log("  FATAL: NO bead row matched an obs_name. The CB strings in")
            log("  the Step05c bead map are not the h5ad obs_names, so every")
            log("  spatial panel would be empty. Compare the two reprs above:")
            log("  a truncation, an appended field, or stray whitespace will be")
            log("  visible immediately. Fix Step05c, do not paper over it here.")
            sys.exit(1)
        pb = pb.dropna(subset=["puck_id"])
    log(f"  coordinates: x [{obs['x_coord'].min():.0f}, {obs['x_coord'].max():.0f}], "
        f"y [{obs['y_coord'].min():.0f}, {obs['y_coord'].max():.0f}] microns")
    log("  neoantigen bead-pairs per puck:")
    for k, v in pb["puck_id"].value_counts().items():
        log(f"    {k}: {v}")

    sep("STEP 2: load network paper neoantigens")
    net, tier_col = load_network()

    sep("STEP 3: overlap")
    res = overlap(mut, net, tier_col)
    for lvl in ("locus", "peptide", "gene"):
        if res.get(lvl):
            pd.DataFrame({lvl: res[lvl]}).to_csv(
                f"{OUT}/shared_{lvl}.tsv", sep="\t", index=False)
    if res.get("mutation"):
        pd.DataFrame(res["mutation"], columns=["gene", "hgvs_p"]).to_csv(
            f"{OUT}/shared_neoantigen_mutations.tsv", sep="\t", index=False)
    pd.DataFrame([res["gene_stats"]]).to_csv(
        f"{OUT}/gene_overlap_stats.tsv", sep="\t", index=False)

    sep("STEP 4: spatial clonality")
    bm = pd.read_csv(BEAD_MAP, sep="\t")
    ctc = pick(bm, "Cell_type")
    bead_ct = dict(zip(bm["CB"], bm[ctc])) if ctc else {}
    mut_beads = set(bm["CB"])
    neo_beads = set(pb["CB"])
    log(f"  mutation-carrying beads   : {len(mut_beads):,}")
    log(f"  neoantigen-carrying beads : {len(neo_beads):,}")

    clon = analyze_clonality(pb, obs, mut_beads)
    if len(clon):
        clon.to_csv(f"{OUT}/neoantigen_clonality.tsv", sep="\t", index=False)
        same = clon[clon["same_puck"]]
        log(f"\n  neoantigens carried by >=2 beads: {len(clon)}")
        log(f"    within one puck : {len(same)}")
        log(f"    across pucks    : {len(clon) - len(same)} "
            "(different patients, so recurrence, not clonality)")
        if len(same):
            log("\n  same-puck multi-bead neoantigens:")
            for _, r in same.iterrows():
                log(f"    {r['gene']:<12} {r['hgvs_p']:<18} "
                    f"{int(r['n_beads'])} beads  {r['pucks']}  "
                    f"min pair {r['min_pair_dist_um']:>7.0f} um  "
                    f"p_raw {r['p_raw']:.3f}  p_BH {r['p_bh']:.3f}")
            log("\n    p_raw is a permutation p against random same-size sets of")
            log("    mutation-carrying beads in that puck. At k=2 with ~12 tests")
            log("    nothing survives correction; these are candidates to look")
            log("    at in the tissue, not established clones.")
    else:
        log("  no neoantigen is carried by 2 or more beads")

    hot = find_hotspots(pb, obs)
    if len(hot):
        hot.to_csv(f"{OUT}/neoantigen_hotspots.tsv", sep="\t", index=False)
        log(f"\n  hotspots (>=2 neoantigen beads within {HOTSPOT_RADIUS:.0f} um): "
            f"{len(hot)}")
        for _, h in hot.sort_values("n_beads", ascending=False).head(10).iterrows():
            tag = "CLONAL CANDIDATE" if h["clonal_candidate"] else "distinct mutations"
            log(f"    #{int(h['hotspot_id'])} {h['puck_id']}  "
                f"{int(h['n_beads'])} beads, {int(h['n_distinct_mutations'])} "
                f"mutations, extent {h['max_extent_um']:.0f} um  [{tag}]")
            if h["clonal_candidate"]:
                log(f"        shared: {h['shared_mutations_within_hotspot']}")
        n_clonal = int(hot["clonal_candidate"].sum())
        log(f"\n  hotspots where >=2 beads share the SAME mutation: {n_clonal}")
        if n_clonal == 0:
            log("    Every hotspot holds distinct mutations, so proximity here")
            log("    reflects where callable epithelium sits rather than one")
            log("    expanded clone carrying one neoantigen.")
    else:
        hot = pd.DataFrame()
        log(f"\n  no hotspots at {HOTSPOT_RADIUS:.0f} um")

    log("\n  global clustering test (mean nearest-neighbour distance):")
    stats = []
    for puck in sorted(obs["puck_id"].unique()):
        tgt_idx = [c for c in neo_beads if c in obs.index
                   and obs.at[c, "puck_id"] == puck]
        if len(tgt_idx) < 3:
            continue
        tgt = obs.loc[tgt_idx, ["x_coord", "y_coord"]].to_numpy(float)
        pools = (
            ("epithelial",
             [c for c in obs.index if obs.at[c, "puck_id"] == puck
              and obs.at[c, "unified_annotation"] == "epithelial"]),
            ("mutation-carrying",
             [c for c in mut_beads if c in obs.index
              and obs.at[c, "puck_id"] == puck]),
        )
        for pool_name, pool_idx in pools:
            pool = obs.loc[pool_idx, ["x_coord", "y_coord"]].to_numpy(float)
            o, m, p = perm_cluster(tgt, pool)
            if np.isnan(p):
                continue
            stats.append({"puck": puck, "pool": pool_name,
                          "n_target": len(tgt), "n_pool": len(pool),
                          "obs": o, "null_median": m, "p": p})
            log(f"    {puck} vs {pool_name:<18}: obs {o:>7.0f} um, "
                f"null {m:>7.0f} um, p = {p:.4f}")
    if stats:
        pd.DataFrame(stats).to_csv(f"{OUT}/spatial_clustering_stats.tsv",
                                   sep="\t", index=False)
        log("    'vs mutation-carrying' is the interpretable comparison;")
        log("    'vs epithelial' largely reports where coverage was adequate.")

    sep("STEP 5: figures")
    fig_celltypes(obs)
    fig_spatial_mutations(obs, mut_beads, neo_beads)
    fig_neoantigen_only(obs, neo_beads, hot)
    fig_clonality(obs, pb, clon)
    fig_clustering(stats)
    fig_burden(obs, bead_ct)
    fig_candidates(mut, set(res.get("gene", [])))
    fig_overlap(res, len(mut), len(net))

    (pb.groupby(["puck_id", "gene", "hgvs_p"]).size()
       .reset_index(name="beads")
       .sort_values("beads", ascending=False)
       .to_csv(f"{OUT}/neoantigen_beads_by_puck_gene.tsv", sep="\t", index=False))

    with open(f"{OUT}/step07_report.txt", "w") as f:
        f.write("\n".join(_report))

    sep(f"STEP 07 COMPLETE in {(time.time()-t0)/60:.1f} min")
    log(f"  tables  : {OUT}/")
    log(f"  figures : {FIG}/")


if __name__ == "__main__":
    main()
