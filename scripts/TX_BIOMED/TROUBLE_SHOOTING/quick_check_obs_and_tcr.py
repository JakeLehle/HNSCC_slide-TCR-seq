#!/usr/bin/env python3
"""
quick_check_obs_and_tcr.py
================================================================================
Fast, read-only. Run interactively, no SLURM needed. Answers three things before
Step08a is locked:

  PART 1  what is actually in adata.obs, so the coordinate and annotation column
          names are confirmed rather than assumed
  PART 2  the clone-size distribution, so MIN_BEADS_FOR_KERNEL is set from the
          data instead of picked
  PART 3  a preflight preview, so we know Step08a's hard gate will pass before
          submitting it

    conda activate slide-TCR-seq
    python quick_check_obs_and_tcr.py

Uses backed='r' so only obs is read, not the 99,341 x 19,822 matrix.
================================================================================
"""

import os
import numpy as np
import pandas as pd
import anndata as ad

PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
ADATA = f"{PROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"
NEO = f"{PROOT}/data/outputs/07_neoantigen/neoantigens_per_bead.tsv"
TCR_DIR = f"{PROOT}/data/inputs/tcr/processed"
PUCKS = ["29", "37", "40"]

pd.set_option("display.width", 200)


def hr(t):
    print("\n" + "=" * 72)
    print(t)
    print("=" * 72)


# ==============================================================================
hr("PART 1  adata.obs schema")
# ==============================================================================
adata = ad.read_h5ad(ADATA, backed="r")
print(f"shape: {adata.n_obs:,} beads x {adata.n_vars:,} genes")
print(f"\nobs columns ({len(adata.obs.columns)}):")
for c in adata.obs.columns:
    s = adata.obs[c]
    nuq = s.nunique(dropna=True)
    if nuq <= 12:
        vals = "; ".join(str(v) for v in sorted(s.dropna().unique(), key=str)[:12])
    else:
        vals = f"{nuq:,} distinct, e.g. " + "; ".join(str(v) for v in s.dropna().unique()[:3])
    print(f"  {c:<34} {str(s.dtype):<12} {vals[:95]}")

print(f"\nobsm keys: {list(adata.obsm.keys())}")
for k in adata.obsm.keys():
    print(f"  {k}: shape {adata.obsm[k].shape}")

print("\nwhat Step08a needs:")
for c in ("unified_annotation", "consensus_annotation", "x_coord", "y_coord"):
    print(f"  {c:<24} {'PRESENT' if c in adata.obs.columns else 'ABSENT'}")

if "x_coord" in adata.obs.columns:
    print(f"\ncoordinate ranges (microns):")
    print(f"  x: {adata.obs.x_coord.min():.1f} to {adata.obs.x_coord.max():.1f}")
    print(f"  y: {adata.obs.y_coord.min():.1f} to {adata.obs.y_coord.max():.1f}")
    # do obs columns agree with obsm['spatial']?
    if "spatial" in adata.obsm:
        sp = np.asarray(adata.obsm["spatial"])
        dx = np.abs(sp[:, 0] - adata.obs.x_coord.values).max()
        dy = np.abs(sp[:, 1] - adata.obs.y_coord.values).max()
        print(f"  max |obs - obsm['spatial']|: x {dx:.4f}, y {dy:.4f}  "
              f"({'identical' if max(dx, dy) < 1e-6 else 'DIFFER, pick one deliberately'})")

print("\nbeads per puck (from obs_names):")
pk = pd.Series([n.split("_", 1)[-1] for n in adata.obs_names]).value_counts()
for k, v in pk.items():
    print(f"  {k}: {v:,}")

if "unified_annotation" in adata.obs.columns:
    print("\nunified_annotation:")
    for k, v in adata.obs.unified_annotation.value_counts().items():
        print(f"  {k:<18} {v:>8,}")

obs_names = set(adata.obs_names)
obs_xy = adata.obs[["x_coord", "y_coord"]].copy() if "x_coord" in adata.obs.columns else None

# ==============================================================================
hr("PART 2  TCR clone-size distribution")
# ==============================================================================
print("The median beads-per-clone will be 1, because most clones are seen once.")
print("A threshold 'below the median' is therefore not available. What matters")
print("is how many clones survive each cutoff and how much signal they carry,")
print("so that is what is tabulated here.\n")

frames = []
for p in PUCKS:
    f = os.path.join(TCR_DIR, f"B59_{p}_hTCR_tcr.csv")
    d = pd.read_csv(f, sep=None, engine="python")
    d.columns = [c.strip() for c in d.columns]
    d["puck"] = p
    d["bead"] = d["bc"].astype(str) + f"-1_Puck_211214_{p}"
    frames.append(d)
tcr = pd.concat(frames, ignore_index=True)
print(f"total rows: {len(tcr):,}")

for p in PUCKS:
    d = tcr[tcr.puck == p]
    cl = d.groupby("cloneId").agg(n_beads=("bead", "nunique"),
                                  n_umi=("umi", "nunique"),
                                  n_reads=("n_reads", "sum"))
    tot_c, tot_u = len(cl), int(cl.n_umi.sum())
    print(f"\n--- puck {p}: {tot_c:,} clones, {tot_u:,} UMIs, "
          f"{d.bead.nunique():,} beads")
    print(f"    beads/clone: median {cl.n_beads.median():.0f}, "
          f"mean {cl.n_beads.mean():.2f}, max {cl.n_beads.max():,}")
    print(f"    {'cutoff':>7} {'clones':>9} {'% clones':>9} {'UMIs':>10} {'% UMIs':>8}")
    for k in (1, 2, 3, 4, 5, 8, 10, 20, 50):
        sub = cl[cl.n_beads >= k]
        if not len(sub):
            break
        print(f"    >={k:<5d} {len(sub):>9,} {100*len(sub)/tot_c:>8.2f}% "
              f"{int(sub.n_umi.sum()):>10,} {100*sub.n_umi.sum()/tot_u:>7.2f}%")

    print(f"    beads/clone percentiles: " + ", ".join(
        f"p{int(q*100)}={cl.n_beads.quantile(q):.0f}" for q in (0.5, 0.9, 0.99, 0.999)))

print("\n--- UMIs per bead (the other candidate threshold) ---")
for p in PUCKS:
    d = tcr[tcr.puck == p]
    per = d.groupby("bead")["umi"].nunique()
    print(f"  puck {p}: median {per.median():.0f}, mean {per.mean():.2f}, "
          f"max {per.max():,}, "
          + ", ".join(f"p{int(q*100)}={per.quantile(q):.0f}" for q in (0.9, 0.99, 0.999)))
    for k in (1, 2, 3, 5, 10):
        n = int((per >= k).sum())
        print(f"    beads with >= {k:2d} UMI: {n:>7,} ({100*n/len(per):5.2f}%)  "
              f"carrying {100*per[per >= k].sum()/per.sum():5.2f}% of UMIs")

# ==============================================================================
hr("PART 3  preflight preview")
# ==============================================================================
neo = pd.read_csv(NEO, sep="\t")
print(f"{NEO}")
print(f"  {len(neo):,} rows, columns: {list(neo.columns)}")

beads = neo["CB"].astype(str).unique()
print(f"  unique CB: {len(beads):,}")
miss = [b for b in beads if b not in obs_names]
print(f"  in adata.obs_names: {len(beads)-len(miss):,} / {len(beads):,}")
print(f"  -> preflight membership check would {'PASS' if not miss else 'FAIL'}")
for b in miss[:10]:
    print(f"       MISSING {b}")

print(f"  CB containing N: {sum(1 for b in beads if 'N' in b.split('-')[0])} "
      "(expected, corrected barcodes carry N)")
print(f"  CB with no '-1' suffix: {sum(1 for b in beads if '-' not in b)} (expect 0)")
print("  puck distribution:")
for k, v in neo.CB.map(lambda b: b.split('_', 1)[-1]).value_counts().items():
    print(f"    {k}: {v}")

# the stronger check: do the coordinates in this file agree with adata?
if obs_xy is not None:
    ok = [b for b in beads if b in obs_names]
    sub = neo[neo.CB.isin(ok)].drop_duplicates("CB").set_index("CB")
    j = sub.join(obs_xy, rsuffix="_adata")
    dx = (j.x_coord - j.x_coord_adata).abs()
    dy = (j.y_coord - j.y_coord_adata).abs()
    bad = int(((dx > 0.5) | (dy > 0.5)).sum())
    print(f"\n  coordinate agreement vs adata: max dx {dx.max():.3f}, dy {dy.max():.3f}")
    print(f"  beads disagreeing by > 0.5 um: {bad}")
    print("  -> a restored barcode landing on the WRONG real bead would show here,")
    print("     which pure membership testing cannot detect")

# cross-puck recurrence, the pooled-calling artifact
print("\n  mutations per puck count (pooled calling artifact check):")
mut = neo.assign(puck=neo.CB.map(lambda b: b.split('_', 1)[-1])) \
         .groupby(["gene", "hgvs_p"])["puck"].nunique().value_counts().sort_index()
for k, v in mut.items():
    print(f"    appearing in {k} puck(s): {v} mutations")
print("  Variants were called on the POOLED BAM, so one variant can be genotyped")
print("  against beads from several pucks. Cross-puck appearance is not")
print("  independent recurrence across patients.")

adata.file.close()
print("\ndone")
