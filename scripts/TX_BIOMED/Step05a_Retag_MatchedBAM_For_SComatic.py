#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step05a_Retag_MatchedBAM_For_SComatic.py
=========================================================================
Build one pooled, SComatic-ready BAM from Sophia Liu's three delivered
matched.bam files, plus the cell-type meta files SplitBamCellTypes.py needs.

WHY matched.bam AND NOT OUR RE-ALIGNMENT
  matched.bam merges BOTH transcriptome flowcells, including the deep 60nt
  H52J2DMXY run whose raw data was never delivered. Pooled it holds 357.4M
  reads at 69.0% MAPQ255, against ~62.7M at ~49% for our HLGH2BGXK-only
  re-alignment (Step01 v2): roughly 12x more usable reads on the same
  annotated beads. Our re-alignment is retained as an independent
  validation set, not discarded.

  Its alignments are also markedly cleaner, because Drop-seq uses stricter
  STAR defaults than the --outFilterMatchNminOverLread 0.3 we carried:

                        our re-alignment    matched.bam
      alignment <25bp        9.69%             0.89%
      alignment <30bp       18.24%             1.84%
      NM/length >10%        13.74%             0.79%

  TRADEOFF, to be stated in Methods: this is the Broad Drop-seq alignment
  (Drop-seq tools 2.4.0 + STAR, GRCh38.102 + GRCh38.102.gtf), not ours.
  Variant calls inherit their aligner parameters. This is exactly the
  reproducibility cost that Gap A (missing raw H52J2DMXY) imposes.

THREE TRANSFORMS APPLIED
  1. CONTIG RENAME. matched.bam uses Ensembl naming (1, 2, MT); our
     reference FASTA is chr-prefixed GENCODE. M5 checksums were compared
     across all 194 contigs: 193 match exactly, sequence-identical with
     only the label differing. Done by rewriting SQ names in the header
     dict. reference_id is a POSITIONAL index and contig order is
     unchanged, so alignment records need no edit. Verified after writing.

     The 194th contig is Y (Ensembl M5 ce3e3110..., ours b2b7e636...),
     almost certainly PAR masking, which Ensembl applies and GENCODE does
     not. Not chased: Y carries 96,433 of 121.7M reads (0.079%) and no
     Tier 1/2/3 neoantigen from the network paper sits on it. Excluded via
     EXCLUDE_CONTIGS rather than risk reference-base mismatches.

  2. nM TAG. SComatic's --max_nM filter reads STAR's nM tag, which the
     Drop-seq output does not carry: a tag census over 200,000 reads found
     NM on every read and nM on none. NM (edit distance, includes indels)
     is copied into nM. Slightly conservative, which is the safe direction.

  3. CB TAG. XB already holds the CORRECTED bead barcode with the -1
     suffix; XC holds the OBSERVED one (verified: 27,640/27,640 sampled XC
     values are in barcode_matching column 1, only 16,281 in column 2).
     CB is set to XB + "_" + puck_id, which is byte-identical to the h5ad
     obs_name, so the SComatic meta Index matches the BAM CB with no
     lookup table.

     The puck suffix is REQUIRED. Bead barcodes are 14nt with ~50k per
     puck, so pooling three patients without it would collide and silently
     undercount distinct beads in BaseCellCounter.

DUPLICATES ARE DELIBERATELY RETAINED
  matched.bam contains PCR duplicates and SComatic counts reads, not UMI
  families, so a duplicate family carrying the same error contributes
  multiple apparent supporting reads. The beta-binomial model assumes
  independence, so this inflates confidence at low-support sites. Retained
  anyway: the network paper ran SComatic on Cell Ranger BAMs with the same
  property, and deduplicating here would break that comparison. State it
  in Methods rather than fixing it in code.

POOLING THREE PATIENTS
  The three pucks are three different patients. A patient-private somatic
  variant is therefore diluted ~3x in the pooled cell-type counts while
  the comparison cell types stay at zero. Cell-type specificity survives;
  sensitivity drops. Set PER_PUCK=True to emit three separate BAMs and
  meta files instead if the pooled call comes back thin.

OUTPUT (all under OUTDIR)
  pooled.matched.retagged.bam(.bai)   three pucks, CB = h5ad obs_name
  meta_unified_annotation.tsv         Index / Cell_type  <- PRIMARY
  meta_consensus_annotation.tsv       Index / Cell_type  <- alternative
  celltype_read_counts.tsv            reads per cell type, BOTH columns
  retag_stats.tsv                     per-puck read accounting
  Step05a_summary.txt                 human-readable summary

Env: slide-TCR-seq (pysam, scanpy, pandas, samtools)
Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import time
import subprocess
import collections
import multiprocessing as mp

import pysam
import pandas as pd
import scanpy as sc

# =========================================================================
# CONFIGURATION
# =========================================================================

PROOT  = "/master/jlehle/WORKING/slide-TCR-seq-working"
INPUTS = f"{PROOT}/data/inputs/fastq"
H5AD   = f"{PROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"
OUTDIR = f"{PROOT}/data/outputs/05_mutations"

PUCK_DIRS = {
    "Puck_211214_29": "2022-01-28_Puck_211214_29",
    "Puck_211214_37": "2022-01-28_Puck_211214_37",
    "Puck_211214_40": "2022-01-28_Puck_211214_40",
}

# --- Annotation columns -------------------------------------------------
# PRIMARY first. unified_annotation (9 categories, Step04f per-bead marker
# argmax) is primary because consensus_annotation scatters epithelial beads
# across 15 separate categories (epithelial 46,837 plus 14 ambiguous:*
# variants), and SComatic flags variants seen in multiple cell types. A true
# tumor variant would appear in up to 15 cell types and be filtered out.
ANNOT_COLS = ["unified_annotation", "consensus_annotation"]

# Categories below this many beads collapse into COLLAPSE_LABEL, which stays
# in the meta file so those reads are EXCLUDED from the real groups rather
# than leaking into them. Near-noop for unified_annotation (smallest real
# category is mast at 724); matters for consensus_annotation.
MIN_BEADS      = 500
COLLAPSE_LABEL = "low_n_other"

# --- CB format ----------------------------------------------------------
# EVERY SComatic script truncates the barcode at the first hyphen:
#     SplitBamCellTypes.py:74, BaseCellCounter.py:243,
#     SitesPerCell.py:177, SingleCellGenotype.py:165 and :204
# all run `barcode = barcode.split("-")[0]` (standard 10x "-1" handling).
#
# Our CB is the h5ad obs_name, "TATATGATTCTCGT-1_Puck_211214_29", so that
# split discards BOTH the suffix and the puck tag and leaves the bare 14 nt
# barcode. Two consequences:
#   1. Output labels lose the puck. Step05c recovers them from this BAM.
#   2. Worse, reads from beads sharing a bare barcode ACROSS PUCKS are
#      pooled during SplitBam and BaseCellCounter. That is cross-patient
#      merging, and no downstream relabelling can undo it.
#
# Measured here: 13 colliding barcodes, 26 of 99,341 beads (0.03%), none
# carrying a called variant. Immaterial on this pilot. It will not stay
# immaterial at full depth with more beads and more calls.
#
# HYPHEN_FREE_CB writes the CB and the meta Index as
#     "TATATGATTCTCGT.1_Puck_211214_29"
# which survives split("-")[0] intact: no truncation, no pooling. Step05c
# converts the dot back to a hyphen to recover the obs_name, and still
# handles the legacy hyphenated form, so both work.
#
# Takes effect only on a full Step05a + Step05b rerun. There is no reason
# to redo 3 hours of counting for 26 beads on the current data.
HYPHEN_FREE_CB = True


def to_cb(obs_name):
    """h5ad obs_name -> the CB written into the BAM and the meta Index."""
    return obs_name.replace("-", ".") if HYPHEN_FREE_CB else obs_name

# --- Read filters -------------------------------------------------------
MIN_MAPQ      = 255    # STAR unique. Verified: MAPQ255 count == NH==1 count.
MAX_NM        = 5      # SComatic --max_nM default, applied here on NM
MIN_ALIGN_LEN = 25     # costs 0.89% on this BAM; blocks 15bp-anchor variants
MAX_NM_RATE   = 0.10   # costs 0.79%; catches short-but-divergent alignments
EXCLUDE_CONTIGS = ("chrY",)   # post-rename names. () to include.

# --- Contig rename: Ensembl -> GENCODE. Scaffolds already share names. ---
RENAME = {str(i): f"chr{i}" for i in range(1, 23)}
RENAME.update({"X": "chrX", "Y": "chrY", "MT": "chrM"})

# --- Execution ----------------------------------------------------------
PER_PUCK    = False    # True = 3 separate BAMs + meta files, no pooling
N_PROCS     = 3        # one per puck; each is a single streaming pass
THREADS     = 16       # samtools merge / index
KEEP_TEMPS  = False    # keep tmp_{puck}.retagged.bam after merge

# =========================================================================


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] [Step05a] {msg}", flush=True)


def sanitize(s):
    """SComatic cell-type labels: alphanumeric and underscore only."""
    out = "".join(c if (c.isalnum() or c == "_") else "_" for c in str(s))
    while "__" in out:
        out = out.replace("__", "_")
    return out.strip("_")


def build_renamed_header(bam):
    """
    Relabel SQ names. Order and lengths are preserved so reference_id
    stays valid without touching alignment records.
    """
    h = bam.header.to_dict()
    old_sq = list(h["SQ"])
    n = 0
    for sq in h["SQ"]:
        if sq["SN"] in RENAME:
            sq["SN"] = RENAME[sq["SN"]]
            n += 1

    # Order/length invariant: this is what makes the positional index safe.
    assert len(h["SQ"]) == len(old_sq), "contig count changed"
    for a, b in zip(h["SQ"], old_sq):
        assert a["LN"] == b["LN"], f"length changed for {a['SN']}"

    h.setdefault("CO", []).append(
        "Step05a: SQ relabelled Ensembl->GENCODE (M5-verified, 193/194); "
        "CB set from XB+puck_id; nM copied from NM; "
        f"filters MAPQ>={MIN_MAPQ} NM<={MAX_NM} alnlen>={MIN_ALIGN_LEN} "
        f"NM/len<={MAX_NM_RATE}; excluded {','.join(EXCLUDE_CONTIGS) or 'none'}"
    )
    return pysam.AlignmentHeader.from_dict(h), n


def retag_puck(task):
    """
    Stream one matched.bam, apply the three transforms, write a temp BAM.
    Annotation-agnostic: membership is tested against the h5ad bead set,
    which is identical for every annotation column. Returns per-bead read
    counts so the parent can tabulate per cell type for ANY column.
    """
    puck, subdir, beads = task
    src = f"{INPUTS}/{subdir}/{puck}.matched.bam"
    dst = f"{OUTDIR}/tmp_{puck}.retagged.bam"

    if not os.path.exists(src):
        return puck, None, None, f"FATAL: {src} not found"

    if os.path.exists(dst) and os.path.getsize(dst) > 0:
        log(f"  {puck}: {os.path.basename(dst)} exists, skipping retag")
        return puck, dst, None, None

    bam = pysam.AlignmentFile(src, "rb", check_sq=False)
    hdr, n_renamed = build_renamed_header(bam)
    out = pysam.AlignmentFile(dst, "wb", header=hdr)

    # tids of contigs to drop, resolved against the ORIGINAL names
    drop_tids = set()
    for old, new in RENAME.items():
        if new in EXCLUDE_CONTIGS:
            tid = bam.get_tid(old)
            if tid >= 0:
                drop_tids.add(tid)

    s = collections.Counter()
    per_bead = collections.Counter()
    xb_checked = 0

    for r in bam.fetch(until_eof=True):
        s["total"] += 1

        if r.is_unmapped:
            s["unmapped"] += 1
            continue
        if r.reference_id in drop_tids:
            s["excluded_contig"] += 1
            continue
        if r.mapping_quality < MIN_MAPQ:
            s["low_mapq"] += 1
            continue

        nm = r.get_tag("NM") if r.has_tag("NM") else 0
        if nm > MAX_NM:
            s["high_nm"] += 1
            continue

        alen = r.query_alignment_length or 0
        if MIN_ALIGN_LEN and alen < MIN_ALIGN_LEN:
            s["short_aln"] += 1
            continue
        if MAX_NM_RATE and alen and (nm / alen) > MAX_NM_RATE:
            s["high_nm_rate"] += 1
            continue

        if not r.has_tag("XB"):
            s["no_xb"] += 1
            continue
        xb = r.get_tag("XB")

        # Sanity-check XB format on the first reads that reach this point.
        if xb_checked < 100:
            xb_checked += 1
            if not xb.endswith("-1"):
                out.close(); bam.close()
                return puck, None, None, (
                    f"FATAL: XB '{xb}' lacks the -1 suffix. Expected the "
                    "corrected barcode. Has the tag convention changed?")

        name = f"{xb}_{puck}"
        if name not in beads:
            s["not_annotated"] += 1
            continue
        r.set_tag("CB", to_cb(name), value_type="Z")
        r.set_tag("nM", int(nm), value_type="i")
        out.write(r)
        s["kept"] += 1
        per_bead[name] += 1

    out.close()
    bam.close()

    if s["kept"] == 0:
        return puck, None, None, f"FATAL: {puck} kept zero reads"

    # Read back and confirm the relabel actually landed on disk.
    chk = pysam.AlignmentFile(dst, "rb", check_sq=False)
    refs = set(chk.references)
    chk.close()
    if "chr1" not in refs:
        return puck, None, None, (
            f"FATAL: {dst} header not relabelled (no chr1 in references)")

    s["contigs_renamed"] = n_renamed
    return puck, dst, (dict(s), dict(per_bead)), None


def write_meta_files(obs):
    """
    One meta file per annotation column. Returns {col: {obs_name: label}}
    so the parent can tabulate reads per cell type for each.
    """
    maps = {}
    lines = []
    for col in ANNOT_COLS:
        vc = obs[col].value_counts()
        keep = set(vc[vc >= MIN_BEADS].index)
        labels = obs[col].astype(str).map(
            lambda v: sanitize(v) if v in keep else COLLAPSE_LABEL)

        path = f"{OUTDIR}/meta_{col}.tsv"
        # Index must match the BAM CB byte for byte, so it goes through the
        # same to_cb() transform. Recover the obs_name with .replace(".","-").
        pd.DataFrame({"Index": [to_cb(i) for i in obs.index],
                      "Cell_type": labels.values}).to_csv(
            path, sep="\t", index=False)
        maps[col] = dict(zip(obs.index, labels.values))

        n_keep = len(keep)
        n_coll = len(vc) - n_keep
        pct = 100.0 * vc[vc >= MIN_BEADS].sum() / vc.sum()
        msg = (f"{col}: {n_keep} categories kept covering {pct:.2f}% of beads, "
               f"{n_coll} collapsed into {COLLAPSE_LABEL}")
        log(msg)
        lines.append(msg)
        for ct, n in pd.Series(labels).value_counts().items():
            log(f"    {ct}: {n:,} beads")
    return maps, lines


def main():
    t0 = time.time()
    os.makedirs(OUTDIR, exist_ok=True)

    log("=" * 70)
    log("STEP 05a: RETAG matched.bam FOR SComatic")
    log("=" * 70)
    log(f"  primary annotation : {ANNOT_COLS[0]}")
    log(f"  filters            : MAPQ>={MIN_MAPQ}, NM<={MAX_NM}, "
        f"alnlen>={MIN_ALIGN_LEN}, NM/len<={MAX_NM_RATE}")
    log(f"  excluded contigs   : {', '.join(EXCLUDE_CONTIGS) or 'none'}")
    log(f"  mode               : {'per-puck' if PER_PUCK else 'pooled'}")

    # --- annotation ------------------------------------------------------
    log(f"Loading {H5AD}")
    adata = sc.read_h5ad(H5AD)
    obs = adata.obs
    log(f"  {adata.n_obs:,} beads, {adata.n_vars:,} genes")

    for c in ANNOT_COLS + ["puck_id"]:
        if c not in obs.columns:
            sys.exit(f"FATAL: column '{c}' not in adata.obs")

    beads = set(obs.index)
    cb_beads = {to_cb(b) for b in beads}
    if HYPHEN_FREE_CB:
        log("  CB format: hyphen-free "
            f"(e.g. {to_cb(obs.index[0])}) so SComatic cannot truncate it")
    for puck in PUCK_DIRS:
        n = int((obs["puck_id"] == puck).sum())
        log(f"  {puck}: {n:,} annotated beads")
        if n == 0:
            sys.exit(f"FATAL: no beads for {puck}; check puck_id values")

    log("Writing meta files ...")
    label_maps, meta_lines = write_meta_files(obs)

    # --- retag in parallel, one process per puck -------------------------
    tasks = [(p, d, beads) for p, d in PUCK_DIRS.items()]
    log(f"Retagging {len(tasks)} pucks with {min(N_PROCS, len(tasks))} processes ...")
    log("  (streaming pass over ~357M reads; expect 1-3 hours)")

    with mp.Pool(min(N_PROCS, len(tasks))) as pool:
        results = pool.map(retag_puck, tasks)

    parts, all_stats, per_bead_all = [], {}, collections.Counter()
    for puck, path, payload, err in results:
        if err:
            sys.exit(err)
        parts.append(path)
        if payload is None:
            log(f"  {puck}: resumed from existing temp BAM, stats unavailable")
            continue
        s, per_bead = payload
        all_stats[puck] = s
        per_bead_all.update(per_bead)
        kept, total = s.get("kept", 0), max(s.get("total", 1), 1)
        log(f"  {puck}: {total:,} records -> {kept:,} kept ({100*kept/total:.1f}%)")
        for k in ("unmapped", "excluded_contig", "low_mapq", "high_nm",
                  "short_aln", "high_nm_rate", "no_xb", "not_annotated"):
            if s.get(k):
                log(f"      dropped {k}: {s[k]:,}")

    if all_stats:
        pd.DataFrame(all_stats).fillna(0).astype("int64").to_csv(
            f"{OUTDIR}/retag_stats.tsv", sep="\t")
        log(f"  wrote {OUTDIR}/retag_stats.tsv")

    # --- reads per cell type, for EVERY annotation column ----------------
    ct_lines = []
    if per_bead_all:
        rows = []
        for col, lab in label_maps.items():
            agg = collections.Counter()
            beads_hit = collections.Counter()
            for bead, n in per_bead_all.items():
                ct = lab.get(bead)
                if ct is None:
                    continue
                agg[ct] += n
                beads_hit[ct] += 1
            for ct, n in agg.items():
                rows.append({"annotation_column": col, "cell_type": ct,
                             "reads": n, "beads_with_reads": beads_hit[ct],
                             "reads_per_bead": round(n / beads_hit[ct], 1)})
        df = pd.DataFrame(rows).sort_values(
            ["annotation_column", "reads"], ascending=[True, False])
        df.to_csv(f"{OUTDIR}/celltype_read_counts.tsv", sep="\t", index=False)
        log(f"  wrote {OUTDIR}/celltype_read_counts.tsv")

        log("Reads per cell type (PRIMARY column):")
        sub = df[df["annotation_column"] == ANNOT_COLS[0]]
        for _, r in sub.iterrows():
            line = (f"    {r['cell_type']}: {r['reads']:,} reads across "
                    f"{r['beads_with_reads']:,} beads "
                    f"({r['reads_per_bead']}/bead)")
            log(line)
            ct_lines.append(line)

    # --- merge and index -------------------------------------------------
    finals = []
    if PER_PUCK:
        for puck, path in zip(PUCK_DIRS, parts):
            final = f"{OUTDIR}/{puck}.matched.retagged.bam"
            if not os.path.exists(final):
                os.rename(path, final)
            subprocess.run(["samtools", "index", "-@", str(THREADS), final],
                           check=True)
            finals.append(final)
            log(f"  {final}")
    else:
        pooled = f"{OUTDIR}/pooled.matched.retagged.bam"
        if not os.path.exists(pooled):
            log("Merging three pucks ...")
            subprocess.run(
                ["samtools", "merge", "-@", str(THREADS), "-f", pooled, *parts],
                check=True)
        log("Indexing ...")
        subprocess.run(["samtools", "index", "-@", str(THREADS), pooled],
                       check=True)
        finals.append(pooled)

        if not KEEP_TEMPS:
            for p in parts:
                if os.path.exists(p) and p != pooled:
                    os.remove(p)
            log("  removed per-puck temp BAMs")

    # --- final verification ---------------------------------------------
    log("Verifying output ...")
    chk = pysam.AlignmentFile(finals[0], "rb")
    refs = list(chk.references)
    assert "chr1" in refs, "FATAL: pooled BAM is not chr-prefixed"
    n_cb, n_nm, n_bad = 0, 0, 0
    for i, r in enumerate(chk.fetch(until_eof=True)):
        if i >= 10000:
            break
        if r.has_tag("CB"):
            n_cb += 1
            if r.get_tag("CB") not in cb_beads:
                n_bad += 1
        if r.has_tag("nM"):
            n_nm += 1
    chk.close()
    log(f"  contigs: {len(refs)} (chr1 present)")
    log(f"  first 10,000 reads: CB on {n_cb}, nM on {n_nm}, "
        f"CB not in h5ad on {n_bad}")
    if n_cb < 10000 or n_nm < 10000 or n_bad > 0:
        sys.exit("FATAL: output verification failed")

    size_gb = os.path.getsize(finals[0]) / 1e9
    log(f"  {finals[0]} ({size_gb:.1f} GB)")

    # --- summary file ----------------------------------------------------
    with open(f"{OUTDIR}/Step05a_summary.txt", "w") as f:
        f.write("Step05a: retag matched.bam for SComatic\n")
        f.write(f"Completed: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Primary annotation column: {ANNOT_COLS[0]}\n")
        f.write(f"Filters: MAPQ>={MIN_MAPQ} NM<={MAX_NM} "
                f"alnlen>={MIN_ALIGN_LEN} NM/len<={MAX_NM_RATE}\n")
        f.write(f"Excluded contigs: {', '.join(EXCLUDE_CONTIGS) or 'none'}\n")
        f.write(f"Duplicates: RETAINED (see header docstring)\n\n")
        for line in meta_lines:
            f.write(line + "\n")
        f.write("\nReads per cell type (primary column):\n")
        for line in ct_lines:
            f.write(line + "\n")
        f.write(f"\nOutput: {finals[0]}\n")

    log("=" * 70)
    log(f"STEP 05a COMPLETE in {(time.time()-t0)/60:.1f} min")
    log("=" * 70)
    log("Read celltype_read_counts.tsv before Step05b: it shows whether each")
    log("cell type has the depth to clear min_cov 5 / min_cells 5.")
    log("")
    log("Next:")
    log(f"  python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py \\")
    log(f"    --bam {finals[0]} \\")
    log(f"    --meta {OUTDIR}/meta_{ANNOT_COLS[0]}.tsv \\")
    log(f"    --id pooled --max_nM {MAX_NM} --max_NH 1 "
        f"--min_MQ {MIN_MAPQ} --outdir {OUTDIR}/SComatic")


if __name__ == "__main__":
    main()
