#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step05b_SComatic_CellType_Variants.py
=========================================================================
Cell-type-level somatic variant calling on the pooled, retagged Slide-seq
BAM produced by Step05a.

PIPELINE
  1. SplitBamCellTypes.py          pooled BAM  -> one BAM per cell type
  2. BaseCellCounter.py            per-cell-type base counts (parallel)
  3. MergeBaseCellCounts.py        one merged count table
  4. BaseCellCalling.step1.py      beta-binomial calling
  5. BaseCellCalling.step2.py      PoN + RNA-editing filtering
  6. bedtools intersect            mappable regions, PASS only
  7. GetAllCallableSites.py        callable sites per cell type
  8. TrinucleotideContextBackground.py   background for signature refitting

  Stops at cell-type-level variants. Per-bead genotyping
  (SingleCellGenotype.py) is deliberately NOT run here; see PER-BEAD below.

INPUT (from Step05a)
  pooled.matched.retagged.bam   224,548,629 reads, CB = h5ad obs_name
  meta_unified_annotation.tsv   99,341 beads, 9 cell types

  Verified depth per cell type (Step05a celltype_read_counts.tsv):
      epithelial     135,993,940 reads / 50,678 beads / 2683.5 per bead
      myeloid         51,323,978 /  25,062 / 2047.9
      B_cell          10,274,314 /   6,056 / 1696.6
      fibroblast      10,183,348 /   5,538 / 1838.8
      T_cell           8,778,648 /   5,184 / 1693.4
      endothelial      3,195,426 /   2,881 / 1109.1
      smooth_muscle    2,986,644 /   2,364 / 1263.4
      mast             1,325,400 /     724 / 1830.7
      ambiguous          486,931 /     854 /  570.2
  All nine clear min_cov 5 comfortably.

PARALLELISM
  BaseCellCounter parallelizes internally by genomic chunk (--nprocs), so
  the useful unit is cores-per-cell-type, not cell-types-at-once. Read
  counts span 280x between epithelial and ambiguous, so a flat allocation
  would leave epithelial running for hours after everything else finished.
  nprocs is therefore allocated PROPORTIONAL to each cell type's read
  count (floor NPROC_MIN, cap NPROC_MAX) and all cell types launch at
  once. Small ones finish early and return their cores to the scheduler.

  This differs from the ClusterCatcher version, which used --nprocs 1 per
  cell type and parallelized across cell types with a ThreadPool. That is
  the right shape when many samples each have modest depth; here there is
  one pooled sample with one dominant cell type.

PARAMETERS MATCH THE NETWORK PAPER
  Everything is at SComatic defaults except --min_bq 30, which is what
  the network paper used. Do not change these casually: the point of the
  exercise is comparing this variant set against the network paper's
  Tier 1/2/3 neoantigen lists, and that comparison is only clean if the
  calling parameters match.

  ONE CAVEAT WORTH A SENSITIVITY RUN. matched.bam retains PCR duplicates
  (deliberately; see Step05a). At ~2,684 reads/bead against ~660 median
  UMIs/bead, duplication is roughly 4x, so MIN_AC_READS=3 can be met by a
  single molecule sequenced three times. The duplicate-robust filter is
  MIN_AC_CELLS=2, which requires two distinct beads. If reviewers push on
  this, rerun phases 4 onward with MIN_AC_READS=5; phases 1 to 3 are
  unaffected and the checkpoints will skip them.

POOLING THREE PATIENTS
  The three pucks are three different patients. A patient-private somatic
  variant is diluted ~3x in pooled cell-type counts while comparison cell
  types stay at zero, so cell-type specificity survives but sensitivity
  drops. Germline variants are filtered correctly, appearing across all
  cell types once pooled. If the pooled call comes back thin, rerun
  Step05a with PER_PUCK=True and point this script at each puck.

THE ambiguous CELL TYPE
  854 beads / 486,931 reads are labelled 'ambiguous' and kept as their own
  cell type so their reads are EXCLUDED from the real compartments rather
  than leaking in. Residual risk: some ambiguous beads are epithelial, so
  a true tumor variant could surface there too and trip SComatic's
  multiple-cell-type flag. At 570 reads/bead the odds of clearing
  min_ac_cells 2 AND min_ac_reads 3 there are low, so this is accepted
  rather than engineered around. Watch the Cell_types column in the
  step1 output; if epithelial variants are routinely co-flagged in
  ambiguous, drop those beads from the meta file and rerun.

PER-BEAD GENOTYPING (not run here)
  SingleCellGenotype.py was scoped out when we expected ~65 reads/bead.
  The actual figure is 2,684, so per-bead genotyping is now plausible and
  is what per-cell COSMIC signature refitting would require. That is a
  separate decision and a separate script (Step05c), not a silent
  addition here.

Env: slide-TCR-seq (pysam, pandas) + bedtools
Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import sys
import glob
import time
import shutil
import subprocess
import collections

import pandas as pd

# =========================================================================
# CONFIGURATION
# =========================================================================

PROOT   = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTDIR  = f"{PROOT}/data/outputs/05_mutations"
SCOMATIC = "/master/jlehle/WORKING/SComatic"

POOLED_BAM = f"{OUTDIR}/pooled.matched.retagged.bam"
META_FILE  = f"{OUTDIR}/meta_unified_annotation.tsv"
READ_COUNTS = f"{OUTDIR}/celltype_read_counts.tsv"
ANNOT_COL   = "unified_annotation"

GENOME_FA  = f"{PROOT}/data/inputs/ref/GRCh38/GRCh38.primary_assembly.genome.fa"
PON_FILE   = f"{SCOMATIC}/PoNs/PoN.scRNAseq.hg38.tsv"
EDIT_SITES = f"{SCOMATIC}/RNAediting/AllEditingSites.hg38.txt"
BED_FILE   = f"{SCOMATIC}/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed"

SAMPLE_ID = "pooled"

# --- Working directories -------------------------------------------------
SC_DIR      = f"{OUTDIR}/SComatic"
SPLIT_DIR   = f"{SC_DIR}/SplitBam"
COUNTS_DIR  = f"{SC_DIR}/BaseCellCounts"
MERGED_DIR  = f"{SC_DIR}/MergedCounts"
CALL_DIR    = f"{SC_DIR}/VariantCalling"
FILT_DIR    = f"{SC_DIR}/FilteredVariants"
CALLABLE_DIR = f"{SC_DIR}/CellTypeCallableSites"
TRINUC_DIR  = f"{SC_DIR}/TrinucleotideBackground"
TMP_ROOT    = f"{SC_DIR}/tmp"
CKPT_DIR    = f"{SC_DIR}/checkpoints"

# --- SComatic parameters (network paper settings) ------------------------
MIN_BQ        = 30     # network paper value; SComatic default is 20
MIN_MQ        = 255    # STAR unique
MAX_NM        = 5
MAX_NH        = 1
MIN_AC_READS  = 3      # see duplicate caveat in the docstring
MIN_AC_CELLS  = 2      # the duplicate-robust filter
MAX_COV       = 150    # GetAllCallableSites
MIN_CELL_TYPES = 2     # GetAllCallableSites

# --- Parallelism ---------------------------------------------------------
TOTAL_CORES = int(os.environ.get("SLURM_CPUS_PER_TASK", 80))
NPROC_MIN   = 2
NPROC_MAX   = 48

KEEP_TMP = False

# =========================================================================


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] [Step05b] {msg}", flush=True)


def section(title):
    print("", flush=True)
    log("=" * 66)
    log(title)
    log("=" * 66)


def ckpt_done(name):
    return os.path.exists(f"{CKPT_DIR}/{name}.done")


def ckpt_set(name):
    os.makedirs(CKPT_DIR, exist_ok=True)
    open(f"{CKPT_DIR}/{name}.done", "w").close()
    log(f"  checkpoint set: {name}")


def run(cmd, label):
    """Run a subprocess, streaming stderr on failure."""
    log(f"  $ {' '.join(str(c) for c in cmd[:6])} ...")
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       text=True)
    if r.returncode != 0:
        log(f"  FAILED ({label}), exit {r.returncode}")
        for line in (r.stderr or "").strip().splitlines()[-25:]:
            log(f"    | {line}")
        return False
    return True


def allocate_nprocs():
    """
    Cores per cell type, proportional to read count. Falls back to an even
    split if celltype_read_counts.tsv is unreadable.
    """
    try:
        df = pd.read_csv(READ_COUNTS, sep="\t")
        df = df[df["annotation_column"] == ANNOT_COL]
        counts = dict(zip(df["cell_type"], df["reads"]))
    except Exception as e:
        log(f"  WARN: could not read {READ_COUNTS} ({e}); using even split")
        return None, {}

    total = sum(counts.values())
    alloc = {}
    for ct, n in counts.items():
        p = int(round(TOTAL_CORES * n / total))
        alloc[ct] = max(NPROC_MIN, min(NPROC_MAX, p))
    return alloc, counts


# =========================================================================
# PHASE 1: SplitBamCellTypes
# =========================================================================

def phase1_split():
    section("PHASE 1: SplitBamCellTypes")
    if ckpt_done("phase1_split"):
        log("  already done, skipping")
        return True

    os.makedirs(SPLIT_DIR, exist_ok=True)
    ok = run([
        "python", f"{SCOMATIC}/scripts/SplitBam/SplitBamCellTypes.py",
        "--bam", POOLED_BAM,
        "--meta", META_FILE,
        "--id", SAMPLE_ID,
        "--min_MQ", str(MIN_MQ),
        "--max_nM", str(MAX_NM),
        "--max_NH", str(MAX_NH),
        "--outdir", SPLIT_DIR,
    ], "SplitBamCellTypes")
    if not ok:
        return False

    bams = sorted(glob.glob(f"{SPLIT_DIR}/{SAMPLE_ID}.*.bam"))
    if not bams:
        log("  FATAL: no split BAMs produced")
        return False

    log(f"  {len(bams)} cell-type BAMs:")
    for b in bams:
        log(f"    {os.path.basename(b)}  ({os.path.getsize(b)/1e9:.2f} GB)")
    ckpt_set("phase1_split")
    return True


# =========================================================================
# PHASE 2: BaseCellCounter
# =========================================================================

def phase2_count():
    section("PHASE 2: BaseCellCounter")
    bams = sorted(glob.glob(f"{SPLIT_DIR}/{SAMPLE_ID}.*.bam"))
    if not bams:
        log("  FATAL: no split BAMs found")
        return False

    os.makedirs(COUNTS_DIR, exist_ok=True)
    alloc, counts = allocate_nprocs()

    # Launch largest first so the long pole starts immediately.
    def ct_of(b):
        return os.path.basename(b)[len(SAMPLE_ID) + 1:-4]

    bams.sort(key=lambda b: -counts.get(ct_of(b), 0))

    procs = []
    for b in bams:
        ct = ct_of(b)
        if ckpt_done(f"phase2_count_{ct}"):
            log(f"  {ct}: already counted, skipping")
            continue
        n = alloc.get(ct, max(NPROC_MIN, TOTAL_CORES // len(bams))) if alloc \
            else max(NPROC_MIN, TOTAL_CORES // len(bams))
        tmp = f"{TMP_ROOT}/count_{ct}"
        os.makedirs(tmp, exist_ok=True)
        cmd = [
            "python", f"{SCOMATIC}/scripts/BaseCellCounter/BaseCellCounter.py",
            "--bam", b,
            "--ref", GENOME_FA,
            "--chrom", "all",
            "--out_folder", COUNTS_DIR,
            "--min_bq", str(MIN_BQ),
            "--min_mq", str(MIN_MQ),
            "--tmp_dir", tmp,
            "--nprocs", str(n),
        ]
        log(f"  launching {ct}: {counts.get(ct, 0):,} reads, nprocs={n}")
        procs.append((ct, tmp, subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)))

    if not procs:
        log("  all cell types already counted")
        return True

    log(f"  {len(procs)} BaseCellCounter jobs running "
        f"({sum(alloc.get(c, 0) for c, _, _ in procs)} cores allocated)")

    failed = []
    for ct, tmp, p in procs:
        _, err = p.communicate()
        if p.returncode != 0:
            log(f"  FAILED: {ct} (exit {p.returncode})")
            for line in (err or "").strip().splitlines()[-20:]:
                log(f"    | {line}")
            failed.append(ct)
        else:
            log(f"  done: {ct}")
            ckpt_set(f"phase2_count_{ct}")
        if not KEEP_TMP and os.path.exists(tmp):
            shutil.rmtree(tmp, ignore_errors=True)

    if failed:
        log(f"  FATAL: {len(failed)} cell type(s) failed: {', '.join(failed)}")
        return False

    tsvs = glob.glob(f"{COUNTS_DIR}/*.tsv")
    log(f"  {len(tsvs)} count tables written")
    return len(tsvs) > 0


# =========================================================================
# PHASE 3 to 6: merge, call, filter
# =========================================================================

def phase3_merge():
    section("PHASE 3: MergeBaseCellCounts")
    out = f"{MERGED_DIR}/{SAMPLE_ID}.BaseCellCounts.AllCellTypes.tsv"
    if ckpt_done("phase3_merge") and os.path.getsize(out) > 0:
        log("  already done, skipping")
        return True
    os.makedirs(MERGED_DIR, exist_ok=True)
    if not run([
        "python", f"{SCOMATIC}/scripts/MergeCounts/MergeBaseCellCounts.py",
        "--tsv_folder", COUNTS_DIR,
        "--outfile", out,
    ], "MergeBaseCellCounts"):
        return False
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        log(f"  FATAL: merged file missing or empty: {out}")
        return False
    log(f"  {out} ({os.path.getsize(out)/1e9:.2f} GB)")
    ckpt_set("phase3_merge")
    return True


def phase4_call():
    section("PHASE 4: BaseCellCalling step1 and step2")
    os.makedirs(CALL_DIR, exist_ok=True)
    merged = f"{MERGED_DIR}/{SAMPLE_ID}.BaseCellCounts.AllCellTypes.tsv"
    prefix = f"{CALL_DIR}/{SAMPLE_ID}"
    s1 = f"{prefix}.calling.step1.tsv"
    s2 = f"{prefix}.calling.step2.tsv"

    if not ckpt_done("phase4_step1"):
        if not run([
            "python", f"{SCOMATIC}/scripts/BaseCellCalling/BaseCellCalling.step1.py",
            "--infile", merged,
            "--outfile", prefix,
            "--ref", GENOME_FA,
            "--min_ac_reads", str(MIN_AC_READS),
            "--min_ac_cells", str(MIN_AC_CELLS),
        ], "BaseCellCalling.step1"):
            return False
        if not os.path.exists(s1) or os.path.getsize(s1) == 0:
            log(f"  FATAL: {s1} missing or empty")
            return False
        ckpt_set("phase4_step1")
    else:
        log("  step1 already done, skipping")

    if not ckpt_done("phase4_step2"):
        if not run([
            "python", f"{SCOMATIC}/scripts/BaseCellCalling/BaseCellCalling.step2.py",
            "--infile", s1,
            "--outfile", prefix,
            "--editing", EDIT_SITES,
            "--pon", PON_FILE,
        ], "BaseCellCalling.step2"):
            return False
        if not os.path.exists(s2) or os.path.getsize(s2) == 0:
            log(f"  FATAL: {s2} missing or empty")
            return False
        ckpt_set("phase4_step2")
    else:
        log("  step2 already done, skipping")

    for f in (s1, s2):
        n = sum(1 for l in open(f) if not l.startswith("#"))
        log(f"  {os.path.basename(f)}: {n:,} rows")
    return True


def phase5_bed():
    section("PHASE 5: BED filtering (mappable regions, PASS only)")
    os.makedirs(FILT_DIR, exist_ok=True)
    s2 = f"{CALL_DIR}/{SAMPLE_ID}.calling.step2.tsv"
    out = f"{FILT_DIR}/{SAMPLE_ID}.calling.filtered.tsv"
    if ckpt_done("phase5_bed") and os.path.exists(out):
        log("  already done, skipping")
        return True

    cmd = (f"bedtools intersect -header -a {s2} -b {BED_FILE} | "
           f"awk '$1 ~ /^#/ || $6 == \"PASS\"' > {out}")
    log(f"  $ {cmd}")
    r = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        log(f"  FAILED: {r.stderr[-500:]}")
        return False

    n = sum(1 for l in open(out) if not l.startswith("#"))
    log(f"  {n:,} PASS variants in mappable regions")
    if n == 0:
        log("  WARN: zero variants survived. Check the Cell_types and FILTER")
        log("        columns of the step2 output before assuming a null result.")
    ckpt_set("phase5_bed")
    return True


def phase6_callable():
    section("PHASE 6: GetAllCallableSites")
    if ckpt_done("phase6_callable"):
        log("  already done, skipping")
        return True
    os.makedirs(CALLABLE_DIR, exist_ok=True)
    if not run([
        "python", f"{SCOMATIC}/scripts/GetCallableSites/GetAllCallableSites.py",
        "--infile", f"{CALL_DIR}/{SAMPLE_ID}.calling.step1.tsv",
        "--outfile", f"{CALLABLE_DIR}/{SAMPLE_ID}",
        "--max_cov", str(MAX_COV),
        "--min_cell_types", str(MIN_CELL_TYPES),
    ], "GetAllCallableSites"):
        return False
    ckpt_set("phase6_callable")
    return True


def phase7_trinuc():
    section("PHASE 7: TrinucleotideContextBackground")
    if ckpt_done("phase7_trinuc"):
        log("  already done, skipping")
        return True
    os.makedirs(TRINUC_DIR, exist_ok=True)
    listfile = f"{TRINUC_DIR}/step1_files.txt"
    with open(listfile, "w") as f:
        f.write(f"{CALL_DIR}/{SAMPLE_ID}.calling.step1.tsv\n")
    if not run([
        "python",
        f"{SCOMATIC}/scripts/TrinucleotideBackground/TrinucleotideContextBackground.py",
        "--in_tsv", listfile,
        "--out_file", f"{TRINUC_DIR}/{SAMPLE_ID}.trinucleotide_background.tsv",
    ], "TrinucleotideContextBackground"):
        return False
    ckpt_set("phase7_trinuc")
    return True


# =========================================================================
# SUMMARY
# =========================================================================

def summarize():
    section("SUMMARY")
    filt = f"{FILT_DIR}/{SAMPLE_ID}.calling.filtered.tsv"
    lines = []

    if os.path.exists(filt):
        hdr, rows = None, []
        for line in open(filt):
            if line.startswith("##"):
                continue
            if line.startswith("#") and hdr is None:
                hdr = line.lstrip("#").rstrip("\n").split("\t")
                continue
            rows.append(line.rstrip("\n").split("\t"))
        if hdr and rows:
            df = pd.DataFrame(rows, columns=hdr[:len(rows[0])])
            lines.append(f"Total PASS variants: {len(df):,}")
            log(lines[-1])

            for col in ("Cell_types", "Cell_type_Filter", "CellTypes"):
                if col in df.columns:
                    log(f"\nVariants per cell type ({col}):")
                    c = collections.Counter()
                    for v in df[col]:
                        for ct in str(v).replace("|", ",").split(","):
                            ct = ct.strip()
                            if ct and ct not in (".", "nan"):
                                c[ct] += 1
                    for ct, n in c.most_common():
                        line = f"    {ct}: {n:,}"
                        log(line); lines.append(line)
                    break

            for col in ("REF", "ALT"):
                if col not in df.columns:
                    break
            else:
                log("\nSubstitution spectrum:")
                sub = collections.Counter(
                    f"{r}>{a}" for r, a in zip(df["REF"], df["ALT"])
                    if len(str(r)) == 1 and len(str(a)) == 1)
                tot = sum(sub.values()) or 1
                for k, n in sub.most_common(12):
                    line = f"    {k}: {n:,} ({100*n/tot:.1f}%)"
                    log(line); lines.append(line)
                ct_frac = 100 * (sub.get("C>T", 0) + sub.get("G>A", 0)) / tot
                line = (f"    C>T plus G>A: {ct_frac:.1f}%  "
                        "(APOBEC SBS2 is C>T at TCW; see Step05c)")
                log(line); lines.append(line)
        else:
            log("No variants in the filtered output.")
    else:
        log(f"Filtered file not found: {filt}")

    with open(f"{SC_DIR}/Step05b_summary.txt", "w") as f:
        f.write("Step05b: SComatic cell-type-level variant calling\n")
        f.write(f"Completed: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Input BAM:  {POOLED_BAM}\n")
        f.write(f"Meta:       {META_FILE} ({ANNOT_COL})\n")
        f.write(f"Parameters: min_bq {MIN_BQ}, min_mq {MIN_MQ}, "
                f"max_nM {MAX_NM}, max_NH {MAX_NH}, "
                f"min_ac_reads {MIN_AC_READS}, min_ac_cells {MIN_AC_CELLS}\n")
        f.write("Duplicates: RETAINED (see docstring)\n")
        f.write("Pooled across 3 patients; chrY excluded upstream.\n\n")
        for l in lines:
            f.write(l + "\n")
    log(f"\nWrote {SC_DIR}/Step05b_summary.txt")


def main():
    t0 = time.time()
    for d in (SC_DIR, SPLIT_DIR, COUNTS_DIR, MERGED_DIR, CALL_DIR, FILT_DIR,
              CALLABLE_DIR, TRINUC_DIR, TMP_ROOT, CKPT_DIR):
        os.makedirs(d, exist_ok=True)

    section("STEP 05b: SComatic CELL-TYPE VARIANT CALLING")
    log(f"  BAM        : {POOLED_BAM}")
    log(f"  meta       : {META_FILE}  ({ANNOT_COL})")
    log(f"  cores      : {TOTAL_CORES}")
    log(f"  min_bq {MIN_BQ}  min_mq {MIN_MQ}  max_nM {MAX_NM}  max_NH {MAX_NH}")
    log(f"  min_ac_reads {MIN_AC_READS}  min_ac_cells {MIN_AC_CELLS}")

    for f in (POOLED_BAM, META_FILE, GENOME_FA, PON_FILE, EDIT_SITES, BED_FILE):
        if not os.path.exists(f):
            sys.exit(f"FATAL: not found: {f}")

    for name, fn in [("phase1", phase1_split), ("phase2", phase2_count),
                     ("phase3", phase3_merge),  ("phase4", phase4_call),
                     ("phase5", phase5_bed),    ("phase6", phase6_callable),
                     ("phase7", phase7_trinuc)]:
        if not fn():
            sys.exit(f"FATAL: {name} failed after {(time.time()-t0)/60:.1f} min")

    summarize()

    if not KEEP_TMP and os.path.exists(TMP_ROOT):
        shutil.rmtree(TMP_ROOT, ignore_errors=True)

    section(f"STEP 05b COMPLETE in {(time.time()-t0)/60:.1f} min")
    log(f"  variants : {FILT_DIR}/{SAMPLE_ID}.calling.filtered.tsv")
    log(f"  callable : {CALLABLE_DIR}/")
    log(f"  trinuc   : {TRINUC_DIR}/")
    log("")
    log("Next: SnpEff annotation and neoantigen prediction, then overlap")
    log("against the network paper Tier 1/2/3 lists by gene + AA substitution.")


if __name__ == "__main__":
    main()
