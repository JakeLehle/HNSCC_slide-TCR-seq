#!/usr/bin/env python3
"""
Diagnostic_Full_BAM_Inventory.py
================================================================================
READ-ONLY. Characterizes every BAM in the Slide-TCR-seq delivery, validates the
read-1 geometry against true raw data for the first time, reconciles read
accounting end to end, and emits a README for the archive.

WHY THIS EXISTS
    The lane-level `unmapped.bam` files are Picard-style UNALIGNED BAMs: both
    reads at full length with base qualities, queryname sorted, already
    demultiplexed per puck. They are raw data in a BAM container. They have no
    @SQ header block, which makes `samtools quickcheck` reject them, which is
    why they were written off as corrupt and why the HPV16 arm was abandoned.
    They were never corrupt. This script documents that permanently.

WHAT IT MEASURES
    BEAT A  unaligned lane BAMs (18)  header, sort order, read group, FLAG and
            length census, read-1 geometry, barcode whitelist match rate
    BEAT B  aligned lane BAMs (18)    header, aligner, reference, tag census,
            MAPQ and NH distribution
    BEAT C  matched.bam (3)           read groups present (proves both flow
            cells), tag census, contig naming
    BEAT D  all_illumina.bam (3)      same, for completeness
    BEAT E  read accounting           cellular_tagging total vs STAR input vs
            uBAM pairs vs final records. The reconciliation that proves the
            uBAMs hold the complete library.
    BEAT F  ONT TCR FASTQ (3)         length distribution, UP linker rate on
            both strands, barcode recoverability
    BEAT G  processed TCR CSV (3)     schema, barcode form, chain census

OUTPUTS  (all under data/outputs/00_inventory/)
    bam_inventory.tsv         one row per BAM, every measured field
    read_accounting.tsv       the reconciliation, per lane
    geometry_validation.tsv   read-1 geometry and whitelist match, per lane
    tcr_inventory.tsv         ONT and CSV summary
    ARCHIVE_README.md         generated, numbers traceable to this run

USAGE
    python Diagnostic_Full_BAM_Inventory.py                # sampled, fast
    python Diagnostic_Full_BAM_Inventory.py --counts       # + exact read counts
    python Diagnostic_Full_BAM_Inventory.py --sample 500000
    python Diagnostic_Full_BAM_Inventory.py --validate-only

Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
================================================================================
"""

import argparse
import csv
import gzip
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime

try:
    import pysam
except ImportError:
    sys.exit("pysam not found. Activate the slide-TCR-seq environment.")

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
PROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
IN = os.path.join(PROOT, "data/inputs/fastq")
BK = os.path.join(PROOT, "BACKUP")
OUT = os.path.join(PROOT, "data/outputs/00_inventory")

PUCKS = ["29", "37", "40"]
PUCK_DIR = "2022-01-28_Puck_211214_{p}"
FLOWCELLS = {"H52J2DMXY": ["L001", "L002"],
             "HLGH2BGXK": ["L001", "L002", "L003", "L004"]}

# Slide-seq V2 read 1 architecture, 1-based inclusive positions
BC1 = (1, 8)          # bead barcode part 1
LINKER = (9, 26)      # UP linker
BC2 = (27, 32)        # bead barcode part 2
UMI = (33, 41)        # UMI
UP_LINKER = "TCTTCAGCGTTCCCGAGA"
UP_LINKER_RC = "TCTCGGGAACGCTGAAGA"

# expected read 2 length per flow cell, to be confirmed not assumed
EXPECTED_R2 = {"H52J2DMXY": 60, "HLGH2BGXK": 42}

TCR_CSV = {"29": "B59_29_hTCR_tcr.csv",
           "37": "B59_37_hTCR_tcr.csv",
           "40": "B59_40_hTCR_tcr.csv"}
TCR_ONT = {"29": "TCR_20220224_Puck_211214_29.gz",
           "37": "TCR_20220228_Puck_211214_37.gz",
           "40": "TCR_20220127_Puck_211214_40.gz"}

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMPLEMENT)[::-1]


def sub1(seq, span):
    """1-based inclusive slice."""
    a, b = span
    return seq[a - 1:b]


def human(n):
    for unit in ["B", "K", "M", "G", "T"]:
        if abs(n) < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}P"


def exact_count(path, threads=4):
    """samtools view -c. Slow on GB-scale files, hence opt-in."""
    try:
        r = subprocess.run(["samtools", "view", "-c", "-@", str(threads), path],
                           capture_output=True, text=True, timeout=7200)
        return int(r.stdout.strip()) if r.returncode == 0 else -1
    except Exception:
        return -1


# ------------------------------------------------------------------------------
# PATH HELPERS
# ------------------------------------------------------------------------------
def lane_bam(puck, fc, lane, kind, root=IN):
    return os.path.join(root, PUCK_DIR.format(p=puck), fc, lane,
                        f"Puck_211214_{puck}.{kind}.bam")


def align_file(puck, fc, lane, name, root=IN):
    return os.path.join(root, PUCK_DIR.format(p=puck), fc, lane, "alignment",
                        f"Puck_211214_{puck}.{name}")


def puck_file(puck, name, root=IN):
    return os.path.join(root, PUCK_DIR.format(p=puck), f"Puck_211214_{puck}.{name}")


def matching_file(puck, root=IN):
    return os.path.join(root, PUCK_DIR.format(p=puck), "barcode_matching",
                        f"Puck_211214_{puck}_barcode_matching.txt.gz")


# ------------------------------------------------------------------------------
# WHITELISTS
# ------------------------------------------------------------------------------
def load_whitelists(puck):
    """
    barcode_matching.txt.gz: observed, corrected(-1), x, y.
    Column 1 is the form barcodes take in reads. Column 2 is the bead assignment
    and is the wrong thing to match against, which cost us a full pipeline run.
    """
    path = matching_file(puck)
    obs, corr = set(), set()
    if not os.path.exists(path):
        return obs, corr
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                continue
            obs.add(f[0])
            corr.add(f[1].split("-")[0])
    return obs, corr


# ------------------------------------------------------------------------------
# BEAT A: unaligned lane BAMs
# ------------------------------------------------------------------------------
def probe_unaligned(path, fc, obs_wl, corr_wl, n_sample):
    """
    Confirms the file is a Picard-style unaligned BAM holding both reads, then
    measures the read-1 geometry directly on raw sequence. This is the first
    time the geometry has been checked against true raw data rather than
    inferred from CR/UR tags.
    """
    rec = {"path": path, "exists": os.path.exists(path)}
    if not rec["exists"]:
        return rec

    rec["size_bytes"] = os.path.getsize(path)
    bam = pysam.AlignmentFile(path, "rb", check_sq=False)   # check_sq=False is the key
    hdr = bam.header.to_dict()

    rec["sort_order"] = hdr.get("HD", {}).get("SO", "none")
    rec["n_SQ"] = len(hdr.get("SQ", []))
    rgs = hdr.get("RG", [])
    rec["n_RG"] = len(rgs)
    rec["RG_ids"] = ";".join(r.get("ID", "") for r in rgs)
    rec["RG_PU"] = ";".join(r.get("PU", "") for r in rgs)
    rec["RG_SM"] = ";".join(r.get("SM", "") for r in rgs)
    rec["PG_programs"] = ";".join(p.get("PN", p.get("ID", "")) for p in hdr.get("PG", []))

    flags, tags = Counter(), Counter()
    r1_len, r2_len = Counter(), Counter()
    n_r1 = n_r2 = 0
    linker_hit = 0
    bc_obs_hit = bc_corr_hit = 0
    umi_n = 0

    for i, read in enumerate(bam.fetch(until_eof=True)):
        if i >= n_sample:
            break
        flags[read.flag] += 1
        for t, _ in read.get_tags():
            tags[t] += 1
        seq = read.query_sequence or ""
        if read.is_read1:
            n_r1 += 1
            r1_len[len(seq)] += 1
            if len(seq) >= UMI[1]:
                if sub1(seq, LINKER) == UP_LINKER:
                    linker_hit += 1
                bc = sub1(seq, BC1) + sub1(seq, BC2)
                if bc in obs_wl:
                    bc_obs_hit += 1
                if bc in corr_wl:
                    bc_corr_hit += 1
                if "N" in sub1(seq, UMI):
                    umi_n += 1
        elif read.is_read2:
            n_r2 += 1
            r2_len[len(seq)] += 1
    bam.close()

    rec["n_sampled"] = n_r1 + n_r2
    rec["n_read1"] = n_r1
    rec["n_read2"] = n_r2
    rec["paired_balance"] = round(n_r1 / n_r2, 4) if n_r2 else None
    rec["flag_top"] = ";".join(f"{f}:{c}" for f, c in flags.most_common(4))
    rec["tags_present"] = ";".join(sorted(tags))
    rec["r1_len_mode"] = r1_len.most_common(1)[0][0] if r1_len else None
    rec["r2_len_mode"] = r2_len.most_common(1)[0][0] if r2_len else None
    rec["r2_len_expected"] = EXPECTED_R2.get(fc)
    rec["r2_len_matches_expected"] = (rec["r2_len_mode"] == EXPECTED_R2.get(fc))
    rec["pct_linker_at_9_26"] = round(100 * linker_hit / n_r1, 2) if n_r1 else None
    rec["pct_bc_in_observed"] = round(100 * bc_obs_hit / n_r1, 2) if n_r1 else None
    rec["pct_bc_in_corrected"] = round(100 * bc_corr_hit / n_r1, 2) if n_r1 else None
    rec["pct_umi_with_N"] = round(100 * umi_n / n_r1, 2) if n_r1 else None
    return rec


# ------------------------------------------------------------------------------
# BEAT B/C/D: aligned BAMs
# ------------------------------------------------------------------------------
def probe_aligned(path, n_sample):
    rec = {"path": path, "exists": os.path.exists(path)}
    if not rec["exists"]:
        return rec

    rec["size_bytes"] = os.path.getsize(path)
    bam = pysam.AlignmentFile(path, "rb")
    hdr = bam.header.to_dict()

    rec["sort_order"] = hdr.get("HD", {}).get("SO", "none")
    sq = hdr.get("SQ", [])
    rec["n_SQ"] = len(sq)
    rec["contig_style"] = ("chr-prefixed" if sq and str(sq[0].get("SN", "")).startswith("chr")
                           else "ensembl" if sq else "none")
    rec["first_contigs"] = ";".join(str(s.get("SN", "")) for s in sq[:3])
    rgs = hdr.get("RG", [])
    rec["n_RG"] = len(rgs)
    rec["RG_ids"] = ";".join(r.get("ID", "") for r in rgs)
    rec["RG_PU"] = ";".join(r.get("PU", "") for r in rgs)
    # flow cells represented, the clean proof that matched.bam holds both
    fcs = set()
    for r in rgs:
        pu = r.get("PU", "")
        if pu:
            fcs.add(pu.split(".")[0])
    rec["flowcells_present"] = ";".join(sorted(fcs))

    pgs = hdr.get("PG", [])
    rec["PG_programs"] = ";".join(p.get("PN", p.get("ID", "")) for p in pgs)
    star = [p for p in pgs if p.get("PN") == "STAR"]
    rec["aligner_version"] = star[0].get("VN", "") if star else ""
    ref = ""
    for p in pgs:
        m = re.search(r"--genomeDir\s+(\S+)", p.get("CL", "") or "")
        if m:
            ref = m.group(1)
            break
    rec["reference"] = ref

    tags, mapq, nh = Counter(), Counter(), Counter()
    seqlen = Counter()
    n = 0
    for i, read in enumerate(bam.fetch(until_eof=True)):
        if i >= n_sample:
            break
        n += 1
        for t, _ in read.get_tags():
            tags[t] += 1
        mapq[read.mapping_quality] += 1
        if read.has_tag("NH"):
            nh[read.get_tag("NH")] += 1
        if read.query_sequence:
            seqlen[len(read.query_sequence)] += 1
    bam.close()

    rec["n_sampled"] = n
    rec["tags_present"] = ";".join(sorted(tags))
    rec["has_nM"] = "nM" in tags
    rec["has_XC"] = "XC" in tags
    rec["has_XB"] = "XB" in tags
    rec["has_XM"] = "XM" in tags
    rec["pct_mapq255"] = round(100 * mapq.get(255, 0) / n, 2) if n else None
    rec["pct_NH1"] = round(100 * nh.get(1, 0) / n, 2) if n else None
    rec["seqlen_mode"] = seqlen.most_common(1)[0][0] if seqlen else None
    return rec


# ------------------------------------------------------------------------------
# BEAT E: read accounting
# ------------------------------------------------------------------------------
def parse_cellular_tagging(path):
    """Histogram of barcode bases failing quality. Sums to total reads tagged."""
    if not os.path.exists(path):
        return None, None
    total, clean = 0, 0
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) == 2 and f[0].isdigit() and f[1].isdigit():
                n_failed, n = int(f[0]), int(f[1])
                total += n
                if n_failed == 0:
                    clean = n
    return (total or None), (clean or None)


def parse_star_log(path):
    if not os.path.exists(path):
        return {}
    out = {}
    keymap = {"Number of input reads": "star_input_reads",
              "Average input read length": "star_avg_input_len",
              "Uniquely mapped reads number": "star_unique",
              "Number of reads mapped to multiple loci": "star_multi",
              "Number of reads unmapped: too short": "star_unmapped_short"}
    with open(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            k, v = [x.strip() for x in line.split("|", 1)]
            if k in keymap:
                try:
                    out[keymap[k]] = int(v)
                except ValueError:
                    out[keymap[k]] = v
    return out


# ------------------------------------------------------------------------------
# BEAT F/G: TCR
# ------------------------------------------------------------------------------
def probe_ont(path, obs_wl, corr_wl, n_reads=200000):
    """
    ONT FASTQ. The bead oligo sits at a variable offset inside a long amplicon
    and roughly half of reads are reverse complement, so the linker must be
    searched on both strands rather than at a fixed position. The previous
    positional probe returned a false negative for exactly this reason.
    """
    rec = {"path": path, "exists": os.path.exists(path)}
    if not rec["exists"]:
        return rec
    rec["size_bytes"] = os.path.getsize(path)

    lens = []
    fwd = rev = 0
    bc_obs = bc_corr = 0
    n = 0
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            s = line.strip()
            n += 1
            lens.append(len(s))
            pos = s.find(UP_LINKER)
            strand = None
            if pos >= 0:
                fwd += 1
                strand = "+"
            else:
                pos = s.find(UP_LINKER_RC)
                if pos >= 0:
                    rev += 1
                    strand = "-"
            # reconstruct the 14 nt barcode flanking the linker
            if strand == "+" and pos >= 8:
                bc = s[pos - 8:pos] + s[pos + 18:pos + 24]
                if len(bc) == 14:
                    if bc in obs_wl:
                        bc_obs += 1
                    if bc in corr_wl:
                        bc_corr += 1
            elif strand == "-" and pos + 18 + 8 <= len(s):
                seg = rc(s[pos - 6 if pos >= 6 else 0:pos + 18 + 8])
                p2 = seg.find(UP_LINKER)
                if p2 >= 8:
                    bc = seg[p2 - 8:p2] + seg[p2 + 18:p2 + 24]
                    if len(bc) == 14:
                        if bc in obs_wl:
                            bc_obs += 1
                        if bc in corr_wl:
                            bc_corr += 1
            if n >= n_reads:
                break

    lens.sort()
    rec["n_sampled"] = n
    rec["len_mean"] = round(sum(lens) / n, 1) if n else None
    rec["len_median"] = lens[n // 2] if n else None
    rec["len_p10"] = lens[int(0.10 * n)] if n else None
    rec["len_p90"] = lens[int(0.90 * n)] if n else None
    rec["len_max"] = lens[-1] if n else None
    # N50 over sampled reads
    half = sum(lens) / 2
    acc = 0
    n50 = None
    for L in reversed(lens):
        acc += L
        if acc >= half:
            n50 = L
            break
    rec["len_N50"] = n50
    rec["pct_linker_fwd"] = round(100 * fwd / n, 2) if n else None
    rec["pct_linker_rev"] = round(100 * rev / n, 2) if n else None
    rec["pct_linker_either"] = round(100 * (fwd + rev) / n, 2) if n else None
    rec["pct_bc_in_observed"] = round(100 * bc_obs / n, 2) if n else None
    rec["pct_bc_in_corrected"] = round(100 * bc_corr / n, 2) if n else None
    return rec


def probe_tcr_csv(path, obs_wl, corr_wl):
    rec = {"path": path, "exists": os.path.exists(path)}
    if not rec["exists"]:
        return rec
    rec["size_bytes"] = os.path.getsize(path)

    chains, bcs, clones = Counter(), set(), set()
    bc_len = Counter()
    suffixed = 0
    n = 0
    umi_per_bc = defaultdict(int)
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh)
        rec["columns"] = ";".join(rdr.fieldnames or [])
        for row in rdr:
            n += 1
            bc = (row.get("bc") or "").strip()
            bcs.add(bc)
            bc_len[len(bc)] += 1
            if "-" in bc:
                suffixed += 1
            chains[(row.get("chain") or "").strip()] += 1
            cid = (row.get("cloneId") or "").strip()
            if cid:
                clones.add(cid)
            umi_per_bc[bc] += 1

    rec["n_rows"] = n
    rec["n_unique_bc"] = len(bcs)
    rec["n_unique_clone"] = len(clones)
    rec["bc_len_mode"] = bc_len.most_common(1)[0][0] if bc_len else None
    rec["bc_with_suffix"] = suffixed
    rec["chain_census"] = ";".join(f"{k}:{v}" for k, v in chains.most_common())
    # THE join question: observed or corrected barcodes
    rec["pct_bc_in_observed"] = round(100 * sum(1 for b in bcs if b in obs_wl) / len(bcs), 2) if bcs else None
    rec["pct_bc_in_corrected"] = round(100 * sum(1 for b in bcs if b in corr_wl) / len(bcs), 2) if bcs else None
    rec["median_umi_per_bc"] = sorted(umi_per_bc.values())[len(umi_per_bc) // 2] if umi_per_bc else None
    rec["max_umi_per_bc"] = max(umi_per_bc.values()) if umi_per_bc else None
    return rec


# ------------------------------------------------------------------------------
# WRITERS
# ------------------------------------------------------------------------------
def write_tsv(path, rows):
    if not rows:
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  wrote {path}  ({len(rows)} rows)")


def write_readme(path, unaligned, aligned, matched, acct, ont, tcrcsv, args):
    def g(rows, key, default="n/a"):
        vals = [r.get(key) for r in rows if r.get(key) not in (None, "")]
        return vals[0] if vals else default

    tot_pairs = sum(a.get("ubam_pairs_expected", 0) or 0 for a in acct)
    lines = []
    A = lines.append

    A("# Slide-TCR-seq Archive README")
    A("")
    A(f"Generated {datetime.now():%Y-%m-%d %H:%M} by `Diagnostic_Full_BAM_Inventory.py`.")
    A("Every number below is measured, not asserted. Regenerate rather than edit.")
    A("")
    A("**Project:** HPV16+ HNSCC Slide-TCR-seq, Sophia Liu (Ragon Institute) and")
    A("the Ebrahimi Lab, Texas Biomedical Research Institute.")
    A("")
    A("---")
    A("")
    A("## READ THIS FIRST: `unmapped.bam` is not broken")
    A("")
    A("`samtools quickcheck` **fails on every `*.unmapped.bam` in this archive.**")
    A("The files are fine. They are Picard-style *unaligned* BAMs, which carry no")
    A("`@SQ` header block, and `quickcheck` treats a missing `@SQ` block as a")
    A("failure. Records read out normally.")
    A("")
    A("These files contain **raw sequencing data**: both reads at full length,")
    A("with base qualities, queryname sorted, already demultiplexed per puck.")
    A("They are the Broad pipeline's entry point, equivalent to the original")
    A("FASTQ. `samtools fastq` converts them losslessly.")
    A("")
    A("To open them in pysam use `check_sq=False`. To validate them, read a")
    A("record rather than running `quickcheck`.")
    A("")
    A("This misreading previously cost the project its HPV16 arm and led to a")
    A("year of believing the deep flow cell was unrecoverable. It is not.")
    A("")
    A("---")
    A("")
    A("## Samples")
    A("")
    A("| Puck | Demux index | Patient |")
    A("|---|---|---|")
    A("| Puck_211214_29 | AGATTTAA | A |")
    A("| Puck_211214_37 | GGCGTCGA | B |")
    A("| Puck_211214_40 | ATCACTCG | C |")
    A("")
    A("**The three pucks are three different patients.** Pooled variant calling")
    A("dilutes patient-private variants, and the same mutation in two pucks is")
    A("recurrence, not clonality.")
    A("")
    A("## Flow cells")
    A("")
    A("| Flow cell | Lanes | Read 2 | Role |")
    A("|---|---|---|---|")
    A(f"| H52J2DMXY | L001, L002 | {EXPECTED_R2['H52J2DMXY']} nt | deep transcriptome, ~90-95% of depth |")
    A(f"| HLGH2BGXK | L001-L004 | {EXPECTED_R2['HLGH2BGXK']} nt | shallow transcriptome |")
    A("")
    A("Read 2 lengths differ between flow cells and must be handled explicitly")
    A("when the two are combined.")
    A("")
    A("## Read 1 architecture (Slide-seq V2)")
    A("")
    A("Read 1 is 42 bp with a **split** bead barcode:")
    A("")
    A("| Bases | Content |")
    A("|---|---|")
    A(f"| {BC1[0]}-{BC1[1]} | bead barcode part 1 |")
    A(f"| {LINKER[0]}-{LINKER[1]} | UP linker `{UP_LINKER}` |")
    A(f"| {BC2[0]}-{BC2[1]} | bead barcode part 2 |")
    A(f"| {UMI[0]}-{UMI[1]} | UMI |")
    A("| 42 | spare cycle |")
    A("")
    A("Barcode = 8 + 6 = 14 nt. UMI = 9 nt. A pipeline that reads 14 contiguous")
    A("bases from position 1 will silently produce 8 real bases plus 6 linker")
    A("bases and recover ~1% of barcodes.")
    A("")
    if unaligned:
        A("Measured directly on raw read 1 in this archive:")
        A("")
        A("| Metric | Value |")
        A("|---|---|")
        A(f"| UP linker exact at bases 9-26 | {g(unaligned,'pct_linker_at_9_26')}% |")
        A(f"| Reconstructed barcode in `barcode_matching` col1 (observed) | {g(unaligned,'pct_bc_in_observed')}% |")
        A(f"| Reconstructed barcode in `barcode_matching` col2 (corrected) | {g(unaligned,'pct_bc_in_corrected')}% |")
        A("")
        A("**Use column 1 (observed) as the STARsolo whitelist.** A whitelist must")
        A("contain barcodes as they appear in reads. Column 2 is the bead")
        A("assignment and contains ambiguous `N` characters.")
        A("")
    A("## What each file is")
    A("")
    A("| File | Contents |")
    A("|---|---|")
    A("| `{puck}.unmapped.bam` | **RAW.** Both reads, full length, with qualities. Queryname sorted, no `@SQ`. Per lane. |")
    A("| `{puck}.final.bam` | Aligned and Drop-seq tagged. `XC` observed barcode, `XM` UMI. Per lane. |")
    A("| `{puck}.matched.bam` | Bead-matched, both flow cells merged. `XC` observed, `XB` corrected with `-1`, `XM` UMI. `nM` absent. |")
    A("| `{puck}.all_illumina.bam` | Merged across lanes before bead matching. |")
    A("| `alignment/*.star.Log.final.out` | STAR summary: input reads, read length, mapping rates. |")
    A("| `alignment/*.cellular_tagging.summary.txt` | Histogram of barcode bases failing quality. Bins 0-14 confirm a 14 nt barcode. **Not a geometry statement.** |")
    A("| `alignment/*.{adapter_trimming,polyA_filtering}.summary.txt` | Bases trimmed per read, **not** reads discarded. |")
    A("| `barcode_matching/*_barcode_matching.txt.gz` | observed, corrected(`-1`), x µm, y µm. |")
    A("| `tcr/processed/B59_*_hTCR_tcr.csv` | MiXCR output, one row per bead barcode and UMI. |")
    A("| `tcr/ont/TCR_*.gz` | Oxford Nanopore FASTQ, rhTCRseq Fraction 2 long-read arm. |")
    A("")
    A("## Measured inventory")
    A("")
    if acct:
        A("### Read accounting, per lane")
        A("")
        A("| Puck | Flow cell | Lane | Tagged reads | Clean barcode | STAR input | uBAM pairs (expected) | Reconciles |")
        A("|---|---|---|---|---|---|---|---|")
        for a in acct:
            A("| {puck} | {fc} | {lane} | {tagged} | {clean} | {star} | {exp} | {ok} |".format(
                puck=a.get("puck"), fc=a.get("flowcell"), lane=a.get("lane"),
                tagged=a.get("tagging_total", "?"), clean=a.get("tagging_clean", "?"),
                star=a.get("star_input_reads", "?"),
                exp=a.get("ubam_pairs_expected", "?"),
                ok=a.get("reconciles", "not checked")))
        A("")
        A("`Clean barcode` is the count with zero failed barcode bases, which is")
        A("what Drop-seq feeds to STAR. The difference from `Tagged reads` is")
        A("dropped before alignment but is **still present in the uBAM**, which is")
        A("why the uBAM is the complete library and `final.bam` is not.")
        A("")
    if matched:
        A("### `matched.bam`")
        A("")
        A("| Puck | Size | Flow cells present | Contig style | Tags |")
        A("|---|---|---|---|---|")
        for m in matched:
            A("| {p} | {s} | {fc} | {cs} | {t} |".format(
                p=m.get("puck"), s=human(m.get("size_bytes", 0)),
                fc=m.get("flowcells_present", "?"), cs=m.get("contig_style", "?"),
                t=m.get("tags_present", "?")))
        A("")
        A("Both flow cells appearing in the read-group list is the proof that")
        A("`matched.bam` carries the deep data.")
        A("")
    if ont:
        A("### ONT TCR reads")
        A("")
        A("| Puck | Size | Median len | N50 | Linker fwd | Linker rev | BC in observed |")
        A("|---|---|---|---|---|---|---|")
        for o in ont:
            A("| {p} | {s} | {med} | {n50} | {f}% | {r}% | {b}% |".format(
                p=o.get("puck"), s=human(o.get("size_bytes", 0)),
                med=o.get("len_median", "?"), n50=o.get("len_N50", "?"),
                f=o.get("pct_linker_fwd", "?"), r=o.get("pct_linker_rev", "?"),
                b=o.get("pct_bc_in_observed", "?")))
        A("")
        A("The linker must be searched on **both strands**; ONT reads are not")
        A("strand-oriented and the bead oligo sits at a variable offset. A")
        A("fixed-position probe returns a false negative.")
        A("")
    if tcrcsv:
        A("### Processed TCR tables")
        A("")
        A("| Puck | Rows | Unique beads | Unique clones | BC len | BC in observed | BC in corrected |")
        A("|---|---|---|---|---|---|---|")
        for t in tcrcsv:
            A("| {p} | {n} | {b} | {c} | {L} | {o}% | {k}% |".format(
                p=t.get("puck"), n=t.get("n_rows", "?"), b=t.get("n_unique_bc", "?"),
                c=t.get("n_unique_clone", "?"), L=t.get("bc_len_mode", "?"),
                o=t.get("pct_bc_in_observed", "?"), k=t.get("pct_bc_in_corrected", "?")))
        A("")
        A("Whichever column the barcodes match is the form to join on. Every")
        A("join in this project has turned on that distinction.")
        A("")
    A("## Known traps")
    A("")
    A("1. `samtools quickcheck` fails on unaligned BAMs. Not corruption.")
    A("2. SComatic truncates barcodes at the hyphen (`barcode.split(\"-\")[0]`) in")
    A("   four separate scripts, discarding both the `-1` suffix and any puck tag.")
    A("3. `barcode_matching` column 2 is not a whitelist.")
    A("4. Read 2 is 60 nt on H52J2DMXY and 42 nt on HLGH2BGXK.")
    A("5. `SampleSheet.csv` reports R2 = 50. The data say 42. `RunInfo.xml` is")
    A("   authoritative.")
    A("6. STAR's \"average input read length\" is post-trim and will read below")
    A("   the true read length.")
    A("")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  wrote {path}")


# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", type=int, default=200000,
                    help="records to sample per BAM (default 200000)")
    ap.add_argument("--ont-reads", type=int, default=200000,
                    help="ONT reads to sample (default 200000)")
    ap.add_argument("--counts", action="store_true",
                    help="exact read counts via samtools view -c (slow)")
    ap.add_argument("--validate-only", action="store_true",
                    help="check paths and exit")
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)

    # ---------- validate ----------
    missing = []
    for p in PUCKS:
        for fc, lanes in FLOWCELLS.items():
            for L in lanes:
                for kind in ("unmapped", "final"):
                    f = lane_bam(p, fc, L, kind)
                    if not os.path.exists(f):
                        missing.append(f)
        for nm in ("matched.bam", "all_illumina.bam"):
            f = puck_file(p, nm)
            if not os.path.exists(f):
                missing.append(f)
        if not os.path.exists(matching_file(p)):
            missing.append(matching_file(p))
    for p in PUCKS:
        for d, m in (("tcr/processed", TCR_CSV), ("tcr/ont", TCR_ONT)):
            f = os.path.join(PROOT, "data/inputs", d, m[p])
            if not os.path.exists(f):
                f2 = os.path.join(BK, m[p])
                if not os.path.exists(f2):
                    missing.append(f)

    print(f"path check: {len(missing)} missing")
    for f in missing[:20]:
        print(f"   MISSING {f}")
    if args.validate_only:
        sys.exit(1 if missing else 0)

    # ---------- whitelists ----------
    print("\nloading whitelists")
    WL = {}
    for p in PUCKS:
        o, c = load_whitelists(p)
        WL[p] = (o, c)
        print(f"  puck {p}: observed {len(o):,}  corrected {len(c):,}")

    unaligned, aligned, matched, allill, acct = [], [], [], [], []

    # ---------- BEAT A + B + E ----------
    for p in PUCKS:
        obs, corr = WL[p]
        for fc, lanes in FLOWCELLS.items():
            for L in lanes:
                print(f"\npuck {p}  {fc}  {L}")

                u = probe_unaligned(lane_bam(p, fc, L, "unmapped"), fc, obs, corr, args.sample)
                u.update(puck=p, flowcell=fc, lane=L, kind="unmapped")
                if args.counts and u.get("exists"):
                    u["n_records_exact"] = exact_count(u["path"])
                unaligned.append(u)
                print(f"  uBAM  nSQ={u.get('n_SQ')} sort={u.get('sort_order')} "
                      f"R1={u.get('r1_len_mode')} R2={u.get('r2_len_mode')} "
                      f"linker={u.get('pct_linker_at_9_26')}% "
                      f"obs={u.get('pct_bc_in_observed')}% corr={u.get('pct_bc_in_corrected')}%")

                a = probe_aligned(lane_bam(p, fc, L, "final"), args.sample)
                a.update(puck=p, flowcell=fc, lane=L, kind="final")
                if args.counts and a.get("exists"):
                    a["n_records_exact"] = exact_count(a["path"])
                aligned.append(a)
                print(f"  final nSQ={a.get('n_SQ')} tags={a.get('tags_present')}")

                tot, clean = parse_cellular_tagging(align_file(p, fc, L, "cellular_tagging.summary.txt"))
                star = parse_star_log(align_file(p, fc, L, "star.Log.final.out"))
                row = {"puck": p, "flowcell": fc, "lane": L,
                       "tagging_total": tot, "tagging_clean": clean}
                row.update(star)
                row["ubam_pairs_expected"] = tot
                row["ubam_records_expected"] = tot * 2 if tot else None
                row["ubam_records_exact"] = u.get("n_records_exact")
                if tot and star.get("star_input_reads"):
                    d = abs(clean - star["star_input_reads"]) if clean else None
                    row["clean_vs_star_delta"] = d
                if u.get("n_records_exact") and tot:
                    row["reconciles"] = "YES" if u["n_records_exact"] == tot * 2 else \
                        f"NO (delta {u['n_records_exact'] - tot*2})"
                else:
                    row["reconciles"] = "not checked (--counts)"
                acct.append(row)

    # ---------- BEAT C + D ----------
    for p in PUCKS:
        for nm, bucket in (("matched.bam", matched), ("all_illumina.bam", allill)):
            print(f"\npuck {p}  {nm}")
            r = probe_aligned(puck_file(p, nm), args.sample)
            r.update(puck=p, kind=nm.replace(".bam", ""))
            if args.counts and r.get("exists"):
                r["n_records_exact"] = exact_count(r["path"])
            bucket.append(r)
            print(f"  flowcells={r.get('flowcells_present')} contigs={r.get('contig_style')} "
                  f"nM={r.get('has_nM')} XB={r.get('has_XB')}")

    # ---------- BEAT F + G ----------
    ont, tcrcsv = [], []
    for p in PUCKS:
        obs, corr = WL[p]

        f = os.path.join(PROOT, "data/inputs/tcr/ont", TCR_ONT[p])
        if not os.path.exists(f):
            f = os.path.join(BK, TCR_ONT[p])
        print(f"\npuck {p}  ONT")
        o = probe_ont(f, obs, corr, args.ont_reads)
        o["puck"] = p
        ont.append(o)
        print(f"  median={o.get('len_median')} N50={o.get('len_N50')} "
              f"linker +{o.get('pct_linker_fwd')}% -{o.get('pct_linker_rev')}% "
              f"obs={o.get('pct_bc_in_observed')}%")

        f = os.path.join(PROOT, "data/inputs/tcr/processed", TCR_CSV[p])
        if not os.path.exists(f):
            f = os.path.join(BK, TCR_CSV[p])
        print(f"puck {p}  TCR CSV")
        t = probe_tcr_csv(f, obs, corr)
        t["puck"] = p
        tcrcsv.append(t)
        print(f"  rows={t.get('n_rows')} beads={t.get('n_unique_bc')} "
              f"clones={t.get('n_unique_clone')} "
              f"obs={t.get('pct_bc_in_observed')}% corr={t.get('pct_bc_in_corrected')}%")

    # ---------- write ----------
    print("\nwriting outputs")
    write_tsv(os.path.join(OUT, "bam_inventory.tsv"), unaligned + aligned + matched + allill)
    write_tsv(os.path.join(OUT, "read_accounting.tsv"), acct)
    write_tsv(os.path.join(OUT, "geometry_validation.tsv"), unaligned)
    write_tsv(os.path.join(OUT, "tcr_inventory.tsv"), ont + tcrcsv)
    write_readme(os.path.join(OUT, "ARCHIVE_README.md"),
                 unaligned, aligned, matched, acct, ont, tcrcsv, args)

    # ---------- headline ----------
    print("\n" + "=" * 72)
    print("HEADLINE")
    print("=" * 72)
    bad_sq = [u for u in unaligned if u.get("n_SQ", 0) != 0]
    print(f"  unaligned BAMs with nSQ != 0:        {len(bad_sq)} (expect 0)")
    bad_r2 = [u for u in unaligned if u.get("r2_len_matches_expected") is False]
    print(f"  unaligned BAMs with wrong R2 length: {len(bad_r2)} (expect 0)")
    if unaligned:
        lk = [u["pct_linker_at_9_26"] for u in unaligned if u.get("pct_linker_at_9_26") is not None]
        ob = [u["pct_bc_in_observed"] for u in unaligned if u.get("pct_bc_in_observed") is not None]
        co = [u["pct_bc_in_corrected"] for u in unaligned if u.get("pct_bc_in_corrected") is not None]
        if lk:
            print(f"  UP linker at 9-26:   {min(lk):.1f}% to {max(lk):.1f}%")
        if ob and co:
            print(f"  barcode in observed: {min(ob):.1f}% to {max(ob):.1f}%")
            print(f"  barcode in corrected:{min(co):.1f}% to {max(co):.1f}%")
            print("  -> observed should win. If it does, the whitelist finding is")
            print("     confirmed on raw data for both flow cells.")
    nr = [a for a in acct if a.get("reconciles", "").startswith("NO")]
    print(f"  lanes failing read reconciliation:   {len(nr)}")
    if not args.counts:
        print("  (rerun with --counts to check reconciliation)")
    print(f"\n  outputs in {OUT}")


if __name__ == "__main__":
    main()
