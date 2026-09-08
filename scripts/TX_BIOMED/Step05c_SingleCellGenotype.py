#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step05c_SingleCellGenotype.py
=========================================================================
Per-bead genotyping of the Step05b variant set, producing the table that
signature_analysis.py consumes and that the neoantigen work needs to know
WHICH beads carry which variants.

PIPELINE
  1. SingleCellGenotype.py per cell-type BAM (parallel, proportional nprocs)
  2. Filter and annotate: observed base must match ALT, minimum ALT reads,
     minimum depth, then add REF_TRI / ALT_TRI FROM THE GENOME
  3. Feasibility report: mutations-per-bead distribution
  4. SitesPerCell.py -> complete_callable_sites.tsv  (optional, slow)

OUTPUT (05_mutations/SComatic/SingleCell/)
  {cell_type}.single_cell_genotype.tsv            raw, per cell type
  FilteredSingleCellAlleles/
      {cell_type}.single_cell_genotype.filtered.tsv
      all_cell.single_cell_genotype.filtered.tsv  <- feed to signature_analysis.py
  mutations_per_bead.tsv                          per-bead counts + cell type
  Step05c_feasibility.txt                         convergence readout
  CombinedCallableSites/complete_callable_sites.tsv   (if enabled)

-------------------------------------------------------------------------
THREE PLACES THE ClusterCatcher LOGIC WOULD BREAK HERE
-------------------------------------------------------------------------
(a) CB SUFFIX. ClusterCatcher does
        df['CB'] = df['CB'] + f"-1-{sample_id}"
    because its BAM carried raw 10x barcodes that needed a sample tag.
    Our CB is ALREADY the full h5ad obs_name ("BARCODE-1_Puck_211214_29"),
    written by Step05a. Appending anything orphans every bead from the
    AnnData join. NO SUFFIX IS APPLIED HERE. Do not add one.

(b) TRINUCLEOTIDE CONTEXT. ClusterCatcher builds REF_TRI from the step2
    Up_context / Down_context fields, taking index [1] of each. The
    network paper established that SComatic's context fields are
    genome-correct in only 23 of 355 cases. This script fetches the
    trinucleotide from GRCh38 with pysam instead, and reports the
    disagreement rate against SComatic's fields so the effect is visible
    rather than assumed.

    SComatic's Start column may be 0-based or 1-based depending on
    version. Rather than guess, detect_coord_base() tests both offsets
    against REF on the first 200 variants and uses whichever agrees.
    (Measured: 1-based, matching the network paper. Offset 0.)

(c) SComatic TRUNCATES THE CB AT THE FIRST HYPHEN.  <-- v2 FIX
    SingleCellGenotype.py applies standard 10x "-1" suffix handling, so
        'TATATGATTCTCGT-1_Puck_211214_29'  (31 chars, the obs_name)
    comes back as
        'TATATGATTCTCGT'                   (14 chars, the bare barcode)
    losing both the suffix and the puck tag. v1 wrote those 14-character
    strings straight into the output, and every downstream join to the
    h5ad silently matched zero rows: Step05d reported "435 beads absent
    from callable sites", and Step07 rendered an empty tissue section that
    looked like "no neoantigen-carrying beads" rather than "the barcodes
    do not match".

    RE-APPENDING THE PUCK IS NOT A VALID FIX. Slide-seq bead barcodes are
    14 nt from one combinatorial library, so the same barcode can exist on
    more than one puck and a bare barcode may map to several obs_names.

    restore_cb() resolves it exactly instead: for each surviving row it
    queries the pooled retagged BAM at that variant position and keeps
    full CB tags whose pre-hyphen prefix matches. Only beads with reads at
    that position are candidates, which in practice makes the answer
    unique. Rows that stay ambiguous or resolve to nothing are dropped and
    counted, never guessed.

-------------------------------------------------------------------------
WILL SIGNATURES CONVERGE? (read Step05c_feasibility.txt)
-------------------------------------------------------------------------
  Step05b produced 1,584 epithelial variants across 50,678 epithelial
  beads, which is 0.031 variants per bead. Each variant cleared
  min_ac_cells 2, so it sits in at least 2 beads; expect on the order of
  3,000 to 8,000 bead-variant pairs, most beads carrying exactly one
  mutation.

  A 96-context NNLS fit on a bead with one mutation is degenerate: all
  weight lands on whichever signature has the highest probability in that
  single context. The lever is signature_analysis.py's --mutation-threshold,
  which drops sparse cells before fitting. Phase 3 reports how many beads
  survive thresholds of 1, 3, 5, 10 and 20 so the choice is made on the
  actual distribution rather than an argument.

  If too few beads survive, the per-bead genotype table is still the
  deliverable that matters: it says which beads carry which variants, and
  that is what the neoantigen and spatial colocalization work needs.

-------------------------------------------------------------------------
WHAT THE VARIANT SET LOOKS LIKE (Step05b, for context)
-------------------------------------------------------------------------
  2,431 PASS variants. Collapsed to pyrimidine convention the spectrum is
  T>C 32.0%, T>G 16.5%, C>T 15.2%, C>G 13.7%, C>A 11.3%, T>A 11.3%.
  The T>C excess is consistent with unannotated ADAR A-to-I editing
  surviving SComatic's editing-site filter (A>G 396 and T>C 381 are nearly
  balanced, as strand-symmetric editing predicts). This is not corrected
  here. Any signature refit on this input will assign substantial weight
  to T>C-heavy signatures (SBS16, SBS5) for that reason, so interpret SBS2
  weight relative to that, not in isolation.

Env: slide-TCR-seq (pysam, pandas, numpy)
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

import numpy as np
import pandas as pd
import pysam

# =========================================================================
# CONFIGURATION
# =========================================================================

PROOT    = "/master/jlehle/WORKING/slide-TCR-seq-working"
OUTDIR   = f"{PROOT}/data/outputs/05_mutations"
SC_DIR   = f"{OUTDIR}/SComatic"
SCOMATIC = "/master/jlehle/WORKING/SComatic"

SPLIT_DIR = f"{SC_DIR}/SplitBam"
CALL_DIR  = f"{SC_DIR}/VariantCalling"
FILT_DIR  = f"{SC_DIR}/FilteredVariants"
OUT_SC    = f"{SC_DIR}/SingleCell"
OUT_FILT  = f"{OUT_SC}/FilteredSingleCellAlleles"
CALLABLE  = f"{OUT_SC}/CombinedCallableSites"
TMP_ROOT  = f"{OUT_SC}/tmp"
CKPT_DIR  = f"{OUT_SC}/checkpoints"

META_FILE   = f"{OUTDIR}/meta_unified_annotation.tsv"
READ_COUNTS = f"{OUTDIR}/celltype_read_counts.tsv"
ANNOT_COL   = "unified_annotation"
GENOME_FA   = f"{PROOT}/data/inputs/ref/GRCh38/GRCh38.primary_assembly.genome.fa"
SAMPLE_ID   = "pooled"

# Source of truth for full CB tags. SingleCellGenotype.py truncates CB at the
# first hyphen, so the full obs_name is recovered from this BAM. See (c).
POOLED_BAM  = f"{OUTDIR}/pooled.matched.retagged.bam"

# Genotype against the BED-filtered PASS set (2,431 sites) rather than the
# full step2 output (790,232). Same format, ~325x fewer sites, and it is the
# analysis set. Flip to the step2 file if SingleCellGenotype rejects it.
VARIANT_FILE = f"{FILT_DIR}/{SAMPLE_ID}.calling.filtered.tsv"
VARIANT_FALLBACK = f"{CALL_DIR}/{SAMPLE_ID}.calling.step2.tsv"

# --- Per-bead filters (ClusterCatcher values) ---------------------------
# NOTE on MIN_ALT_READS: matched.bam retains PCR duplicates at roughly 4x
# (2,684 reads/bead vs ~660 median UMIs/bead), so 3 ALT reads can be one
# molecule sequenced three times. There is no per-bead equivalent of
# min_ac_cells to fall back on. Raise to 5 for a stricter pass.
MIN_ALT_READS   = 3
MIN_TOTAL_DEPTH = 5

# --- Execution ----------------------------------------------------------
TOTAL_CORES = int(os.environ.get("SLURM_CPUS_PER_TASK", 80))
NPROC_MIN, NPROC_MAX = 2, 32
RUN_SITES_PER_CELL = True   # slow: iterates the 3.8 GB step1 file per BAM
KEEP_TMP = False

# =========================================================================

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] [Step05c] {m}", flush=True)


def section(t):
    print("", flush=True)
    log("=" * 66)
    log(t)
    log("=" * 66)


def ckpt_done(n):
    return os.path.exists(f"{CKPT_DIR}/{n}.done")


def ckpt_set(n):
    os.makedirs(CKPT_DIR, exist_ok=True)
    open(f"{CKPT_DIR}/{n}.done", "w").close()
    log(f"  checkpoint set: {n}")


def read_scomatic_tsv(path):
    """Parse a SComatic output file: '##' comments, one '#' header, rows."""
    hdr, rows = None, []
    for line in open(path):
        if line.startswith("##"):
            continue
        if line.startswith("#") and hdr is None:
            hdr = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if line.strip():
            rows.append(line.rstrip("\n").split("\t"))
    if hdr is None or not rows:
        return None
    w = len(rows[0])
    return pd.DataFrame([r[:w] for r in rows], columns=hdr[:w])


def allocate_nprocs(cell_types):
    try:
        df = pd.read_csv(READ_COUNTS, sep="\t")
        df = df[df["annotation_column"] == ANNOT_COL]
        counts = dict(zip(df["cell_type"], df["reads"]))
    except Exception as e:
        log(f"  WARN: {READ_COUNTS} unreadable ({e}); even split")
        return {ct: max(NPROC_MIN, TOTAL_CORES // len(cell_types))
                for ct in cell_types}, {}
    total = sum(counts.get(ct, 1) for ct in cell_types) or 1
    return ({ct: max(NPROC_MIN, min(NPROC_MAX,
             int(round(TOTAL_CORES * counts.get(ct, 1) / total))))
             for ct in cell_types}, counts)


def ct_of(bam):
    return os.path.basename(bam)[len(SAMPLE_ID) + 1:-4]


# =========================================================================
# PHASE 1: SingleCellGenotype
# =========================================================================

def phase1_genotype():
    section("PHASE 1: SingleCellGenotype")

    bams = sorted(glob.glob(f"{SPLIT_DIR}/{SAMPLE_ID}.*.bam"))
    if not bams:
        log(f"  FATAL: no split BAMs in {SPLIT_DIR}. Run Step05b first.")
        return False

    vfile = VARIANT_FILE
    if not os.path.exists(vfile) or os.path.getsize(vfile) == 0:
        log(f"  WARN: {vfile} missing; falling back to step2")
        vfile = VARIANT_FALLBACK
    n_var = sum(1 for l in open(vfile) if not l.startswith("#"))
    log(f"  variants: {n_var:,} from {os.path.basename(vfile)}")

    os.makedirs(OUT_SC, exist_ok=True)
    alloc, counts = allocate_nprocs([ct_of(b) for b in bams])
    bams.sort(key=lambda b: -counts.get(ct_of(b), 0))

    procs = []
    for b in bams:
        ct = ct_of(b)
        out = f"{OUT_SC}/{ct}.single_cell_genotype.tsv"
        if ckpt_done(f"phase1_gt_{ct}"):
            log(f"  {ct}: already genotyped, skipping")
            continue
        tmp = f"{TMP_ROOT}/gt_{ct}"
        os.makedirs(tmp, exist_ok=True)
        cmd = [
            "python",
            f"{SCOMATIC}/scripts/SingleCellGenotype/SingleCellGenotype.py",
            "--bam", b,
            "--infile", vfile,
            "--meta", META_FILE,
            "--ref", GENOME_FA,
            "--outfile", out,
            "--tmp_dir", tmp,
            "--nprocs", str(alloc.get(ct, NPROC_MIN)),
        ]
        log(f"  launching {ct}: nprocs={alloc.get(ct, NPROC_MIN)}")
        procs.append((ct, tmp, subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)))

    if not procs:
        log("  all cell types already genotyped")
        return True

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
            ckpt_set(f"phase1_gt_{ct}")
        if not KEEP_TMP:
            shutil.rmtree(tmp, ignore_errors=True)

    if failed:
        log(f"  FATAL: failed for {', '.join(failed)}")
        return False
    return True


# =========================================================================
# PHASE 2: filter and annotate with GENOME trinucleotide context
# =========================================================================

def detect_coord_base(df, fa, chrom_col, pos_col, ref_col):
    """
    SComatic's Start may be 0-based or 1-based. Test both on up to 200
    variants and return the offset whose fetched base agrees with REF.
    """
    trials = {0: 0, 1: 0}
    n = 0
    for _, r in df.head(200).iterrows():
        ref = str(r[ref_col])
        if len(ref) != 1:
            continue
        try:
            p = int(r[pos_col])
        except (TypeError, ValueError):
            continue
        n += 1
        for off in (0, 1):
            try:
                b = fa.fetch(str(r[chrom_col]), p - off, p - off + 1).upper()
            except Exception:
                continue
            if b == ref:
                trials[off] += 1
    best = max(trials, key=trials.get)
    log(f"  coordinate probe on {n} variants: "
        f"0-based {trials[0]}, 1-based {trials[1]} -> using offset {best}")
    if n and trials[best] / n < 0.9:
        log("  WARN: neither offset agrees with REF on >90% of variants.")
        log("        Check that the BAM and FASTA are the same assembly.")
    return best


def pick(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None


def restore_cb(df, chrom_c, pos_c, cb_c, offset, bam):
    """
    Recover the full obs_name for each row. See docstring item (c).

    SingleCellGenotype.py truncates CB at the first hyphen, so a row carries
    only the 14 nt bare barcode. Bead barcodes are not unique across pucks,
    so the puck cannot simply be re-appended. Instead, query the pooled BAM
    at the variant position: only beads with reads there are candidates, and
    those reads still carry the untruncated CB.

    Returns (restored_series, stats). Rows that resolve to zero or more than
    one candidate get None and are dropped by the caller.
    """
    out, stats = [], collections.Counter()
    cache = {}
    for _, r in df.iterrows():
        bare = str(r[cb_c])
        chrom = str(r[chrom_c])
        try:
            p = int(r[pos_c]) - offset          # -> 1-based
        except (TypeError, ValueError):
            out.append(None); stats["bad_position"] += 1; continue

        key = (chrom, p)
        if key not in cache:
            seen = collections.defaultdict(set)
            try:
                for read in bam.fetch(chrom, p - 1, p):
                    if read.has_tag("CB"):
                        full = read.get_tag("CB")
                        seen[full.split("-")[0]].add(full)
            except Exception:
                seen = {}
            cache[key] = seen
        cands = cache[key].get(bare, set())

        if len(cands) == 1:
            out.append(next(iter(cands))); stats["resolved"] += 1
        elif len(cands) > 1:
            out.append(None); stats["ambiguous"] += 1
        else:
            out.append(None); stats["no_read_at_site"] += 1
    return pd.Series(out, index=df.index), stats


def phase2_filter_annotate():
    section("PHASE 2: filter and annotate (genome trinucleotide context)")
    os.makedirs(OUT_FILT, exist_ok=True)

    files = sorted(glob.glob(f"{OUT_SC}/*.single_cell_genotype.tsv"))
    if not files:
        log("  FATAL: no genotype files found")
        return False

    if not os.path.exists(POOLED_BAM):
        log(f"  FATAL: {POOLED_BAM} not found. It is required to recover the")
        log("  full CB tags that SingleCellGenotype.py truncated. See (c).")
        return False
    if not os.path.exists(POOLED_BAM + ".bai"):
        log(f"  FATAL: {POOLED_BAM}.bai missing; CB recovery needs random access")
        return False
    pooled = pysam.AlignmentFile(POOLED_BAM, "rb")

    fa = pysam.FastaFile(GENOME_FA)
    stats = collections.Counter()
    cb_stats = collections.Counter()
    parts, offset = [], None

    for path in files:
        ct = os.path.basename(path).split(".")[0]
        df = read_scomatic_tsv(path)
        if df is None or df.empty:
            log(f"  {ct}: empty, skipping")
            continue

        chrom_c = pick(df, "#CHROM", "CHROM", "Chr")
        pos_c   = pick(df, "Start", "POS", "Pos")
        ref_c   = pick(df, "REF", "Ref")
        alt_c   = pick(df, "ALT_expected", "ALT")
        obs_c   = pick(df, "Base_observed", "ALT_observed", "Base")
        nrd_c   = pick(df, "Num_reads", "ALT_reads", "Nreads")
        dep_c   = pick(df, "Total_depth", "Depth", "DP")
        cb_c    = pick(df, "CB", "Index", "Cell")

        missing = [n for n, c in [("chrom", chrom_c), ("pos", pos_c),
                                  ("REF", ref_c), ("ALT_expected", alt_c),
                                  ("Base_observed", obs_c), ("Num_reads", nrd_c),
                                  ("Total_depth", dep_c), ("CB", cb_c)]
                   if c is None]
        if missing:
            log(f"  FATAL: {ct} missing columns {missing}")
            log(f"    available: {list(df.columns)}")
            return False

        if offset is None:
            offset = detect_coord_base(df, fa, chrom_c, pos_c, ref_c)

        n0 = len(df)
        stats["rows"] += n0

        def alt_set(v):
            if pd.isna(v) or str(v) in (".", ""):
                return set()
            return {b.strip() for b in str(v).replace("|", ",").split(",")
                    if b.strip()}

        df["_alts"] = df[alt_c].map(alt_set)
        m_base = df.apply(lambda r: r[obs_c] in r["_alts"], axis=1)
        stats["drop_base_mismatch"] += int((~m_base).sum())

        nrd = pd.to_numeric(df[nrd_c], errors="coerce").fillna(0)
        dep = pd.to_numeric(df[dep_c], errors="coerce").fillna(0)
        m_rd = nrd >= MIN_ALT_READS
        m_dp = dep >= MIN_TOTAL_DEPTH
        stats["drop_alt_reads"] += int((~m_rd & m_base).sum())
        stats["drop_depth"] += int((~m_dp & m_base & m_rd).sum())

        keep = df[m_base & m_rd & m_dp].drop(columns=["_alts"]).copy()
        if keep.empty:
            log(f"  {ct}: 0 of {n0:,} rows survived")
            continue

        # --- (c) recover the full obs_name from the pooled BAM -----------
        bare_len = keep[cb_c].astype(str).str.len().unique()
        restored, cbs = restore_cb(keep, chrom_c, pos_c, cb_c, offset, pooled)
        cb_stats.update(cbs)
        n_before = len(keep)
        keep[cb_c] = restored
        keep = keep[keep[cb_c].notna()].copy()
        # If Step05a wrote hyphen-free CBs (HYPHEN_FREE_CB), turn the dot back
        # into a hyphen to recover the obs_name. Legacy hyphenated CBs have no
        # dot, so this is a no-op for them and both formats work.
        keep[cb_c] = keep[cb_c].astype(str).str.replace(".", "-", regex=False)
        log(f"  {ct}: CB recovery {cbs['resolved']}/{n_before} resolved "
            f"(bare len {list(bare_len)}, ambiguous {cbs['ambiguous']}, "
            f"no read at site {cbs['no_read_at_site']})")
        if keep.empty:
            log(f"  {ct}: nothing left after CB recovery")
            continue

        # --- trinucleotide context FROM THE GENOME ---
        ref_tri, alt_tri, bad = [], [], 0
        sc_up = pick(keep, "Up_context", "UP_context")
        sc_dn = pick(keep, "Down_context", "DOWN_context")
        disagree = 0

        for _, r in keep.iterrows():
            try:
                p = int(r[pos_c]) - offset
                tri = fa.fetch(str(r[chrom_c]), p - 1, p + 2).upper()
            except Exception:
                tri = ""
            if len(tri) != 3 or tri[1] != str(r[ref_c]):
                ref_tri.append("N" * 3)
                alt_tri.append("N" * 3)
                bad += 1
                continue
            ref_tri.append(tri)
            alt_tri.append(tri[0] + str(r[obs_c]) + tri[2])
            if sc_up and sc_dn:
                u, d = str(r[sc_up]), str(r[sc_dn])
                if len(u) >= 2 and len(d) >= 2:
                    if (u[1] + str(r[ref_c]) + d[1]) != tri:
                        disagree += 1

        keep["REF_TRI"] = ref_tri
        keep["ALT_TRI"] = alt_tri
        stats["bad_context"] += bad
        stats["context_disagree"] += disagree

        # signature_analysis.py requires these exact names.
        if ref_c != "REF":
            keep = keep.rename(columns={ref_c: "REF"})
        if alt_c != "ALT_expected":
            keep = keep.rename(columns={alt_c: "ALT_expected"})
        if cb_c != "CB":
            keep = keep.rename(columns={cb_c: "CB"})
        keep["Cell_type"] = ct

        # NO CB SUFFIX. CB is already the h5ad obs_name. See docstring (a).

        out = f"{OUT_FILT}/{ct}.single_cell_genotype.filtered.tsv"
        keep.to_csv(out, sep="\t", index=False)
        parts.append(keep)
        stats["kept"] += len(keep)
        log(f"  {ct}: {len(keep):,} of {n0:,} rows kept, "
            f"{keep['CB'].nunique():,} beads")

    fa.close()
    pooled.close()
    if not parts:
        log("  FATAL: nothing survived filtering in any cell type")
        return False

    combined = pd.concat(parts, ignore_index=True)
    combined.to_csv(f"{OUT_FILT}/all_cell.single_cell_genotype.filtered.tsv",
                    sep="\t", index=False)

    # --- (c) the join that v1 got wrong: verify it, do not assume it -----
    meta_idx = set(pd.read_csv(META_FILE, sep="\t")["Index"].astype(str))
    joined = combined["CB"].astype(str).isin(meta_idx)
    log("")
    log("CB recovery:")
    log(f"  resolved         : {cb_stats['resolved']:,}")
    log(f"  ambiguous        : {cb_stats['ambiguous']:,} "
        "(same bare barcode on >1 puck at that site; dropped)")
    log(f"  no read at site  : {cb_stats['no_read_at_site']:,} (dropped)")
    log(f"  join to meta     : {int(joined.sum()):,} / {len(combined):,}")
    if not joined.all():
        log("")
        log("  FATAL: restored CB values still do not match the meta Index.")
        log(f"    restored : {combined['CB'].iloc[0]!r}")
        log(f"    expected : {sorted(meta_idx)[0]!r}")
        log("  Do not proceed; every downstream spatial join depends on this.")
        return False
    log("  all restored barcodes match the annotated bead set")

    log("")
    log("Filtering summary:")
    log(f"  rows in            : {stats['rows']:,}")
    log(f"  base mismatch      : {stats['drop_base_mismatch']:,}")
    log(f"  ALT reads < {MIN_ALT_READS}      : {stats['drop_alt_reads']:,}")
    log(f"  depth < {MIN_TOTAL_DEPTH}          : {stats['drop_depth']:,}")
    log(f"  kept               : {stats['kept']:,}")
    log(f"  unresolved context : {stats['bad_context']:,}")
    if stats["context_disagree"]:
        pct = 100 * stats["context_disagree"] / max(stats["kept"], 1)
        log(f"  SComatic context disagrees with genome: "
            f"{stats['context_disagree']:,} ({pct:.1f}%)")
        log("    (genome-derived context is used; this is the reason why)")
    ckpt_set("phase2_filter")
    return True


# =========================================================================
# PHASE 3: feasibility
# =========================================================================

def phase3_feasibility():
    section("PHASE 3: per-bead mutation burden")
    path = f"{OUT_FILT}/all_cell.single_cell_genotype.filtered.tsv"
    df = pd.read_csv(path, sep="\t")
    per = df.groupby("CB").size().sort_values(ascending=False)

    meta = pd.read_csv(META_FILE, sep="\t")
    ct_of_bead = dict(zip(meta["Index"], meta["Cell_type"]))
    n_beads_total = len(meta)

    tbl = pd.DataFrame({"CB": per.index, "n_mutations": per.values})
    tbl["cell_type"] = tbl["CB"].map(ct_of_bead)
    tbl.to_csv(f"{OUT_SC}/mutations_per_bead.tsv", sep="\t", index=False)

    lines = []

    def emit(s):
        log(s)
        lines.append(s)

    emit(f"bead-variant pairs        : {len(df):,}")
    emit(f"beads with >=1 mutation   : {len(per):,} of {n_beads_total:,} "
         f"({100*len(per)/n_beads_total:.2f}%)")
    emit(f"median mutations per bead : {int(per.median())}")
    emit(f"max mutations in one bead : {int(per.max())}")
    emit("")
    emit("Beads surviving --mutation-threshold:")
    for t in (1, 2, 3, 5, 10, 20):
        n = int((per >= t).sum())
        emit(f"    >= {t:>2}: {n:,} beads")
    emit("")
    emit("Beads with >=1 mutation, by cell type:")
    for ct, n in tbl["cell_type"].value_counts().items():
        emit(f"    {ct}: {n:,}")
    emit("")
    emit("A 96-context NNLS fit on a bead with one mutation is degenerate:")
    emit("all weight lands on whichever signature is highest in that single")
    emit("context. Use the table above to pick --mutation-threshold. If few")
    emit("beads clear 5, report signature weights per CELL TYPE and use the")
    emit("per-bead table for neoantigen and spatial work instead.")

    with open(f"{OUT_SC}/Step05c_feasibility.txt", "w") as f:
        f.write("Step05c per-bead mutation burden\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Filters: ALT reads >= {MIN_ALT_READS}, "
                f"depth >= {MIN_TOTAL_DEPTH}\n\n")
        f.write("\n".join(lines) + "\n")
    log(f"\nwrote {OUT_SC}/Step05c_feasibility.txt")
    ckpt_set("phase3_feasibility")
    return True


# =========================================================================
# PHASE 4: SitesPerCell (optional)
# =========================================================================

def phase4_sites_per_cell():
    section("PHASE 4: SitesPerCell (callable sites per bead)")
    if not RUN_SITES_PER_CELL:
        log("  disabled by config, skipping")
        return True
    if ckpt_done("phase4_sites"):
        log("  already done, skipping")
        return True

    step1 = f"{CALL_DIR}/{SAMPLE_ID}.calling.step1.tsv"
    if not os.path.exists(step1):
        log(f"  WARN: {step1} missing, skipping")
        return True

    os.makedirs(CALLABLE, exist_ok=True)
    spc_dir = f"{OUT_SC}/UniqueCellCallableSites"
    os.makedirs(spc_dir, exist_ok=True)
    bams = sorted(glob.glob(f"{SPLIT_DIR}/{SAMPLE_ID}.*.bam"))
    alloc, counts = allocate_nprocs([ct_of(b) for b in bams])
    bams.sort(key=lambda b: -counts.get(ct_of(b), 0))

    log(f"  iterating a {os.path.getsize(step1)/1e9:.1f} GB step1 file per BAM;")
    log("  this is the long pole. Set RUN_SITES_PER_CELL=False to skip.")

    procs = []
    for b in bams:
        ct = ct_of(b)
        tmp = f"{TMP_ROOT}/spc_{ct}"
        os.makedirs(tmp, exist_ok=True)
        cmd = ["python", f"{SCOMATIC}/scripts/SitesPerCell/SitesPerCell.py",
               "--bam", b, "--infile", step1, "--ref", GENOME_FA,
               "--out_folder", spc_dir, "--tmp_dir", tmp,
               "--nprocs", str(alloc.get(ct, NPROC_MIN))]
        procs.append((ct, tmp, subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)))

    for ct, tmp, p in procs:
        _, err = p.communicate()
        if p.returncode != 0:
            log(f"  WARN: SitesPerCell failed for {ct} (exit {p.returncode})")
            for line in (err or "").strip().splitlines()[-10:]:
                log(f"    | {line}")
        else:
            log(f"  done: {ct}")
        if not KEEP_TMP:
            shutil.rmtree(tmp, ignore_errors=True)

    meta = pd.read_csv(META_FILE, sep="\t")
    allc = pd.DataFrame({"CB": meta["Index"], "SitesPerCell": 0})
    got = {}
    for f in glob.glob(f"{spc_dir}/*.tsv"):
        if os.path.getsize(f) < 20:
            continue
        sep = "," if "," in open(f).readline() else "\t"
        d = pd.read_csv(f, sep=sep)
        if "CB" in d.columns and "SitesPerCell" in d.columns:
            got.update(dict(zip(d["CB"], d["SitesPerCell"])))
    allc["SitesPerCell"] = allc["CB"].map(got).fillna(0).astype(int)
    out = f"{CALLABLE}/complete_callable_sites.tsv"
    allc.to_csv(out, sep="\t", index=False)
    n = int((allc["SitesPerCell"] > 0).sum())
    log(f"  {out}: {n:,} of {len(allc):,} beads with callable sites")
    ckpt_set("phase4_sites")
    return True


def main():
    t0 = time.time()
    for d in (OUT_SC, OUT_FILT, CALLABLE, TMP_ROOT, CKPT_DIR):
        os.makedirs(d, exist_ok=True)

    section("STEP 05c: PER-BEAD GENOTYPING")
    log(f"  variants : {VARIANT_FILE}")
    log(f"  meta     : {META_FILE}")
    log(f"  filters  : ALT reads >= {MIN_ALT_READS}, depth >= {MIN_TOTAL_DEPTH}")
    log(f"  cores    : {TOTAL_CORES}")

    for f in (META_FILE, GENOME_FA, f"{GENOME_FA}.fai"):
        if not os.path.exists(f):
            sys.exit(f"FATAL: not found: {f}")

    for name, fn in [("phase1", phase1_genotype),
                     ("phase2", phase2_filter_annotate),
                     ("phase3", phase3_feasibility),
                     ("phase4", phase4_sites_per_cell)]:
        if not fn():
            sys.exit(f"FATAL: {name} failed after {(time.time()-t0)/60:.1f} min")

    if not KEEP_TMP:
        shutil.rmtree(TMP_ROOT, ignore_errors=True)

    section(f"STEP 05c COMPLETE in {(time.time()-t0)/60:.1f} min")
    log(f"  genotypes : {OUT_FILT}/all_cell.single_cell_genotype.filtered.tsv")
    log(f"  per bead  : {OUT_SC}/mutations_per_bead.tsv")
    log(f"  readout   : {OUT_SC}/Step05c_feasibility.txt")
    log("")
    log("Next: signature_analysis.py, e.g.")
    log(f"  --mutations {OUT_FILT}/all_cell.single_cell_genotype.filtered.tsv \\")
    log(f"  --adata {PROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad \\")
    log(f"  --callable-sites {CALLABLE}/complete_callable_sites.tsv \\")
    log("  --hnscc-only --use-scree --core-signatures SBS1 SBS2 SBS5 SBS13 \\")
    log("  --mutation-threshold <pick from Step05c_feasibility.txt>")


if __name__ == "__main__":
    main()
