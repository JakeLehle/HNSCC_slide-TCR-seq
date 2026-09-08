#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step06_Neoantigen_Prediction.py
=========================================================================
Predict MHC-I neoantigens from the spatial SComatic variant set, and map
each one back to the individual beads that carry it.

Ports the network paper's NEOANTIGEN pipeline (Step01 + Step02 + Step03)
into a single script, because the spatial variant set is ~2,400 rows
rather than ~10,000 and does not need three stages or an env switch.

PIPELINE
  1. Load Step05b PASS variants (cell-type resolved) + Step05c per-bead calls
  2. Detect the SComatic coordinate convention against GRCh38, then write VCFs
  3. SnpEff annotate (GRCh38.p14, -canon), parse ANN, keep protein-altering
  4. Non-epithelial subtraction: drop anything also called outside epithelium
  5. Missense -> mutant/wild-type peptides from the Ensembl r115 proteome
  6. MHCflurry Class1AffinityPredictor over the network paper's 10-allele panel
  7. Map neoantigens back to beads, cell types, pucks and spatial coordinates

-------------------------------------------------------------------------
DECISIONS CARRIED FROM THE NETWORK PAPER (do not change casually)
-------------------------------------------------------------------------
COORDINATES. SComatic's `Start` is ALREADY 1-based. VCF POS = Start, with
  NO +1. The network paper's Step01 docstring records that an earlier +1
  shifted every variant one base 3', made SnpEff annotate the wrong base,
  and emitted WARNING_REF_DOES_NOT_MATCH_GENOME on nearly every record
  (SBS2 3296, CNV 2713, NORMAL 3438). detect_offset() below re-derives the
  convention from the genome and aborts if agreement is poor, so this
  cannot silently regress.

HLA PANEL. The same 10 alleles the network paper used, roughly 80%
  population coverage. You asked for ~90%; that would need a wider panel,
  and a wider panel breaks the overlap comparison, because a peptide
  binding an allele the network paper never tested would look
  spatial-specific for a purely methodological reason. Widen only if the
  overlap is abandoned. Flagged, not silently substituted.

PEPTIDES. Real protein context from Ensembl GRCh38 release 115
  pep.all.fa, not poly-alanine padding. Lookup chain is ENST, gene symbol,
  alias, ENSG, isoform scan, then a +/-30 offset scan for signal-peptide
  numbering shifts. The network paper reached 98.6-98.7% mapping this way.

THRESHOLDS. Binder IC50 < 500 nM, strong binder < 50 nM, differential
  neoantigen mutant < 500 nM with wild-type > 500 nM. Lengths 8 to 11.

-------------------------------------------------------------------------
WHAT DIFFERS HERE, AND WHY
-------------------------------------------------------------------------
GROUPS. The network paper split cells into SBS2_HIGH / CNV_HIGH / NORMAL
  from per-cell signature refitting. That is unavailable here: only 6 of
  435 mutation-carrying beads have 2 or more mutations, so per-bead
  refitting is not possible. The axis is therefore CELL TYPE, with
  epithelium as the tumor compartment and the immune and stromal
  compartments as the background. This is a different contrast and must
  be described as such, not as a replication of the group comparison.

SUBTRACTION. The network paper subtracts NORMAL-group variants as
  germline. Here the analogue is non-epithelial subtraction. SComatic
  already assigns each PASS variant to exactly one cell type (the 2,431
  variants sum exactly across the nine types), so this is close to a
  no-op; it is run and reported anyway so the count is auditable.

TWO VARIANT SETS, DELIBERATELY.
  - CATALOG (Step05b, 2,431 PASS variants, cell-type resolved) is what
    gets annotated and scored. It is the fuller set.
  - BEAD MAP (Step05c, 442 bead-variant pairs across 435 beads) says
    which individual beads carry which call. Only calls clearing 3 ALT
    reads and 5x depth in a single bead appear here.
  Every output carries `n_beads`, so a neoantigen with n_beads = 0 is
  real at the cell-type level but has no single-bead support. Do not
  quote bead counts as prevalence: they are a floor set by depth.

SPECTRUM CAVEAT. The Step05b spectrum is 32% T>C with A>G and T>C nearly
  balanced, consistent with unannotated ADAR editing surviving SComatic's
  filter. RNA editing produces real transcript-level changes, so some of
  these peptides may genuinely be presented, but they are not somatic DNA
  mutations. `is_tcw` and `sub_class` are emitted per neoantigen so an
  APOBEC-attributable subset can be separated from the rest.

Env: NEOANTIGEN (SnpEff, mhcflurry, pysam, pandas)
Author: Jake Lehle, Texas Biomedical Research Institute
Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
"""

import os
import re
import sys
import time
import subprocess
import collections

import numpy as np
import pandas as pd
import pysam

# =========================================================================
# CONFIGURATION
# =========================================================================

PROOT   = "/master/jlehle/WORKING/slide-TCR-seq-working"
SC_DIR  = f"{PROOT}/data/outputs/05_mutations/SComatic"
OUT     = f"{PROOT}/data/outputs/07_neoantigen"

CATALOG  = f"{SC_DIR}/FilteredVariants/pooled.calling.filtered.tsv"
BEAD_MAP = f"{SC_DIR}/SingleCell/FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv"
H5AD     = f"{PROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"
GENOME   = f"{PROOT}/data/inputs/ref/GRCh38/GRCh38.primary_assembly.genome.fa"

# Ensembl r115 proteome, shared with the network paper.
PROTEOME = "/master/jlehle/WORKING/2026_NMF_PAPER/data/reference/Homo_sapiens.GRCh38.pep.all.fa"

TUMOR_CELL_TYPE = "epithelial"
SNPEFF_GENOME   = "GRCh38.p14"
SNPEFF_XMX      = "-Xmx200g"

PEPTIDE_LENGTHS      = [8, 9, 10, 11]
MHC_BIND_THRESH      = 500
STRONG_BIND_THRESH   = 50
MAX_POSITION_OFFSET  = 30

# Network paper panel. See the HLA PANEL note in the docstring.
HLA_PANEL = [
    "HLA-A0201", "HLA-A0101", "HLA-A0301", "HLA-A2402",
    "HLA-B0702", "HLA-B0801", "HLA-B4402", "HLA-B3501",
    "HLA-C0701", "HLA-C0401",
]

PROTEIN_ALTERING = [
    "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
    "start_lost", "inframe_insertion", "inframe_deletion",
    "disruptive_inframe_insertion", "disruptive_inframe_deletion",
]

GENE_ALIASES = {"C4orf3": "C4orf33", "TMEM199": "VMA12", "SLC9A3R1": "NHERF1"}

AA_3TO1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
           "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
           "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
           "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
VALID_AA = set(AA_3TO1.values())
STANDARD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]

_report = []


def log(m=""):
    print(m, flush=True)
    _report.append(m)


def sep(t=""):
    log("")
    log("=" * 78)
    if t:
        log(f"  {t}")
        log("=" * 78)


def read_scomatic(path):
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


def pick(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None


# =========================================================================
# COORDINATES
# =========================================================================

def detect_offset(df, fa, chrom_c, pos_c, ref_c):
    """
    Is SComatic's Start 1-based (offset 0) or 0-based (offset 1)?
    The network paper established 1-based. Verified, not assumed.
    """
    hit = {0: 0, 1: 0}
    n = 0
    for _, r in df.head(400).iterrows():
        ref = str(r[ref_c])
        if len(ref) != 1:
            continue
        try:
            p = int(r[pos_c])
        except (TypeError, ValueError):
            continue
        n += 1
        for off in (0, 1):
            try:
                if fa.fetch(str(r[chrom_c]), p - 1 - off, p - off).upper() == ref:
                    hit[off] += 1
            except Exception:
                pass
    best = max(hit, key=hit.get)
    log(f"  coordinate probe on {n} variants: "
        f"1-based {hit[0]}, 0-based {hit[1]} -> offset {best}")
    if n == 0 or hit[best] / n < 0.95:
        log("  FATAL: neither convention matches REF on >95% of variants.")
        log("  Check that the BAM, the FASTA and the variant table agree.")
        sys.exit(1)
    if best != 0:
        log("  WARNING: this run says 0-based, but the network paper")
        log("  established 1-based. Reconcile before trusting the output.")
    return best


def write_vcf(df, path, chrom_c, pos_c, ref_c, alt_c, offset, sample="SPATIAL"):
    d = df[df[chrom_c].isin(STANDARD_CHROMS)].copy()
    d["_pos"] = pd.to_numeric(d[pos_c], errors="coerce") - offset
    d = d.dropna(subset=["_pos"])
    d = d.drop_duplicates(subset=[chrom_c, "_pos", ref_c, alt_c])
    order = {c: i for i, c in enumerate(STANDARD_CHROMS)}
    d = d.sort_values([chrom_c, "_pos"], key=lambda s: s.map(order)
                      if s.name == chrom_c else s)
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n##source=SComatic_spatial\n")
        f.write("##reference=GRCh38\n")
        for c in STANDARD_CHROMS:
            f.write(f"##contig=<ID={c}>\n")
        f.write('##INFO=<ID=CT,Number=1,Type=String,Description="Cell type">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for _, r in d.iterrows():
            ct = str(r.get("Cell_types", ".")).replace(";", ",") or "."
            f.write(f"{r[chrom_c]}\t{int(r['_pos'])}\t.\t{r[ref_c]}\t"
                    f"{r[alt_c]}\t.\tPASS\tCT={ct}\tGT\t0/1\n")
    return len(d)


# =========================================================================
# SnpEff
# =========================================================================

def run_snpeff(vcf_in, vcf_out):
    if os.path.exists(vcf_out) and os.path.getsize(vcf_out) > 100:
        log(f"    {os.path.basename(vcf_out)} exists, reusing")
        return True
    cmd = ["snpEff", SNPEFF_XMX, "ann", "-noStats", "-no-downstream",
           "-no-upstream", "-no-intergenic", "-canon", SNPEFF_GENOME, vcf_in]
    log(f"    $ {' '.join(cmd[:5])} ... {os.path.basename(vcf_in)}")
    try:
        with open(vcf_out, "w") as f:
            r = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE,
                               text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        log("    SnpEff timed out after 3600s")
        return False
    if r.returncode != 0:
        log(f"    SnpEff failed: {(r.stderr or '')[:500]}")
        return False
    return True


def parse_snpeff(vcf):
    """Parse the first (canonical) ANN entry per record."""
    out, warn = [], 0
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        info = f[7]
        if "WARNING_REF_DOES_NOT_MATCH_GENOME" in info:
            warn += 1
        ann = next((p[4:] for p in info.split(";") if p.startswith("ANN=")), "")
        if not ann:
            continue
        a = ann.split(",")[0].split("|")
        if len(a) < 11:
            continue
        ct = next((p[3:] for p in info.split(";") if p.startswith("CT=")), "")
        out.append({"chrom": f[0], "pos": f[1], "ref": f[3], "alt": f[4],
                    "effect": a[1], "impact": a[2], "gene": a[3],
                    "gene_id": a[4], "transcript_id": a[6],
                    "hgvs_c": a[9], "hgvs_p": a[10], "cell_type": ct})
    return pd.DataFrame(out), warn


# =========================================================================
# PROTEOME
# =========================================================================

def load_proteome(path):
    if not os.path.exists(path):
        gz = path + ".gz"
        if os.path.exists(gz):
            import gzip
            fh, path = gzip.open(gz, "rt"), gz
        else:
            log(f"  FATAL: proteome not found at {path}")
            log("  Download Ensembl release-115 pep.all.fa; see the network")
            log("  paper's RUN_NEOANTIGEN_PIEPLINE.sh pre-flight check.")
            sys.exit(1)
    else:
        fh = open(path)

    by_enst, by_symbol, by_ensg, isoforms = {}, {}, {}, collections.defaultdict(list)
    name = seq = None
    meta = {}

    def commit():
        if not name or not seq:
            return
        e = dict(meta); e["seq"] = seq; e["length"] = len(seq); e["ensp"] = name
        if e.get("enst"):
            by_enst[e["enst"]] = e
        sym = e.get("symbol")
        if sym:
            isoforms[sym].append(e)
            cur = by_symbol.get(sym)
            if (cur is None
                    or (e["is_canonical"] and not cur["is_canonical"])
                    or (e["is_canonical"] == cur["is_canonical"]
                        and e["length"] > cur["length"])):
                by_symbol[sym] = e
        g = e.get("ensg")
        if g:
            cur = by_ensg.get(g)
            if (cur is None
                    or (e["is_canonical"] and not cur["is_canonical"])
                    or (e["is_canonical"] == cur["is_canonical"]
                        and e["length"] > cur["length"])):
                by_ensg[g] = e

    for line in fh:
        line = line.rstrip("\n")
        if line.startswith(">"):
            commit()
            name = line[1:].split()[0].split(".")[0]
            seq = ""
            enst = re.search(r"transcript:(\S+)", line)
            ensg = re.search(r"gene:(\S+)", line)
            sym  = re.search(r"gene_symbol:(\S+)", line)
            meta = {"enst": enst.group(1).split(".")[0] if enst else None,
                    "ensg": ensg.group(1).split(".")[0] if ensg else None,
                    "symbol": sym.group(1) if sym else None,
                    "is_canonical": "Ensembl_canonical" in line}
        else:
            seq += line
    commit()
    fh.close()
    log(f"  proteome: {len(by_enst):,} transcripts, {len(by_symbol):,} symbols, "
        f"{len(by_ensg):,} genes")
    for g in ("TP53", "ANXA1", "KRT6B", "COX4I1", "SPRR1A", "B2M"):
        e = by_symbol.get(g)
        log(f"    {g:<8} {e['length']:>5} aa ({e['enst']})" if e
            else f"    {g:<8} NOT FOUND")
    return by_enst, by_symbol, by_ensg, isoforms


def parse_hgvs_p(h):
    if not h or pd.isna(h):
        return None
    m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", str(h).strip())
    if not m:
        return None
    wt, pos, mt = AA_3TO1.get(m.group(1)), int(m.group(2)), AA_3TO1.get(m.group(3))
    return (wt, pos, mt) if wt and mt else None


def peptides_from_seq(seq, wt, pos, mt, lengths, status):
    mut = seq[:pos - 1] + mt + seq[pos:]
    mp, wp, meta = [], [], []
    for L in lengths:
        for s in range(max(0, pos - L), min(pos, len(seq) - L + 1)):
            w, m = seq[s:s + L], mut[s:s + L]
            if "*" in w or "*" in m:
                continue
            if not (set(w) <= VALID_AA and set(m) <= VALID_AA):
                continue
            mp.append(m); wp.append(w)
            meta.append({"peptide_length": L,
                         "mut_position_in_peptide": (pos - 1) - s,
                         "protein_start": s + 1, "protein_end": s + L})
    status["n_peptides"] = len(mp)
    status["status"] = "SUCCESS" if mp else "NO_VALID_PEPTIDES"
    return mp, wp, meta, status


def generate_peptides(sym, ensg, enst, wt, pos, mt, lengths,
                      by_enst, by_symbol, by_ensg, isoforms):
    """Six-layer lookup: ENST, symbol, alias, ENSG, isoform scan, +/-30 offset."""
    st = {"gene_symbol": sym, "gene_id": ensg, "transcript_id": enst,
          "status": None, "detail": "", "lookup_method": None,
          "protein_length": None, "aa_match": None, "n_peptides": 0}
    entry, method = None, None
    for cand, meth in ((by_enst.get((enst or "").split(".")[0]), "enst"),
                       (by_symbol.get(sym), "symbol"),
                       (by_symbol.get(GENE_ALIASES.get(sym, "")), "alias"),
                       (by_ensg.get((ensg or "").split(".")[0]), "ensg")):
        if cand:
            entry, method = cand, meth
            break
    if entry is None:
        st["status"] = "GENE_NOT_FOUND"
        st["detail"] = f"{sym} / {ensg} / {enst} absent from proteome"
        return [], [], [], st

    seq = entry["seq"]
    st["lookup_method"], st["protein_length"] = method, len(seq)
    if 1 <= pos <= len(seq) and seq[pos - 1] == wt:
        st["aa_match"] = True
        return peptides_from_seq(seq, wt, pos, mt, lengths, st)

    tried = 0
    for iso in isoforms.get(sym, []):
        tried += 1
        s = iso["seq"]
        if 1 <= pos <= len(s) and s[pos - 1] == wt:
            st.update({"lookup_method": f"isoform_scan(was:{method})",
                       "aa_match": True, "protein_length": len(s),
                       "detail": f"matched isoform {iso['enst']}"})
            return peptides_from_seq(s, wt, pos, mt, lengths, st)
    st["n_isoforms_tried"] = tried

    for iso in ([entry] + isoforms.get(sym, [])):
        s = iso["seq"]
        for off in range(1, MAX_POSITION_OFFSET + 1):
            for tp in (pos + off, pos - off):
                if 1 <= tp <= len(s) and s[tp - 1] == wt:
                    st.update({"lookup_method": f"offset_scan({tp-pos:+d},was:{method})",
                               "aa_match": True, "protein_length": len(s),
                               "detail": f"offset {tp-pos:+d}: {pos}->{tp} "
                                         f"in {iso['enst']}"})
                    return peptides_from_seq(s, wt, tp, mt, lengths, st)

    if pos < 1 or pos > len(seq):
        st["status"] = "POSITION_OUT_OF_BOUNDS"
        st["detail"] = f"pos {pos} outside {len(seq)} aa, no isoform or offset"
    else:
        st.update({"aa_match": False, "status": "AA_MISMATCH",
                   "detail": f"expected {wt} at {pos}, found {seq[pos-1]} "
                             f"(tried {tried} isoforms, offset +/-{MAX_POSITION_OFFSET})"})
    return [], [], [], st


# =========================================================================
# MAIN
# =========================================================================

def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)

    sep("STEP 06: NEOANTIGEN PREDICTION (spatial)")
    log(f"  catalog  : {CATALOG}")
    log(f"  bead map : {BEAD_MAP}")
    log(f"  tumor    : {TUMOR_CELL_TYPE}")
    log(f"  HLA      : {len(HLA_PANEL)} alleles (network paper panel, ~80%)")

    # --- variants -------------------------------------------------------
    sep("STEP 1: load variants")
    cat = read_scomatic(CATALOG)
    if cat is None:
        sys.exit(f"FATAL: could not parse {CATALOG}")
    chrom_c = pick(cat, "#CHROM", "CHROM")
    pos_c   = pick(cat, "Start", "POS")
    ref_c   = pick(cat, "REF")
    alt_c   = pick(cat, "ALT")
    ct_c    = pick(cat, "Cell_types", "Cell_type")
    if not all([chrom_c, pos_c, ref_c, alt_c]):
        sys.exit(f"FATAL: unexpected columns in catalog: {list(cat.columns)}")
    if ct_c and ct_c != "Cell_types":
        cat = cat.rename(columns={ct_c: "Cell_types"})
    log(f"  catalog: {len(cat):,} PASS variants")
    if "Cell_types" in cat:
        for k, v in cat["Cell_types"].value_counts().items():
            log(f"    {k}: {v:,}")

    beads = pd.read_csv(BEAD_MAP, sep="\t")
    log(f"  bead map: {len(beads):,} pairs, {beads['CB'].nunique():,} beads")

    fa = pysam.FastaFile(GENOME)
    offset = detect_offset(cat, fa, chrom_c, pos_c, ref_c)

    # --- VCFs -----------------------------------------------------------
    sep("STEP 2: write VCFs")
    is_tum = cat["Cell_types"].astype(str).str.strip() == TUMOR_CELL_TYPE
    n_t = write_vcf(cat[is_tum], f"{OUT}/spatial_{TUMOR_CELL_TYPE}.vcf",
                    chrom_c, pos_c, ref_c, alt_c, offset, TUMOR_CELL_TYPE.upper())
    n_o = write_vcf(cat[~is_tum], f"{OUT}/spatial_other.vcf",
                    chrom_c, pos_c, ref_c, alt_c, offset, "OTHER")
    log(f"  {TUMOR_CELL_TYPE}: {n_t:,} unique variants")
    log(f"  non-{TUMOR_CELL_TYPE}: {n_o:,} unique variants")

    # --- SnpEff ---------------------------------------------------------
    sep("STEP 3: SnpEff annotation")
    ann = {}
    for tag in (TUMOR_CELL_TYPE, "other"):
        vin, vout = f"{OUT}/spatial_{tag}.vcf", f"{OUT}/{tag}.snpeff.vcf"
        log(f"  {tag}:")
        if not run_snpeff(vin, vout):
            sys.exit(f"FATAL: SnpEff failed for {tag}")
        df, warn = parse_snpeff(vout)
        log(f"    parsed {len(df):,} annotated variants")
        if warn:
            log(f"    WARNING_REF_DOES_NOT_MATCH_GENOME on {warn:,} records")
            if warn > 0.05 * max(len(df), 1):
                log("    FATAL: >5% reference mismatch. This is the signature")
                log("    of a coordinate shift. Do not proceed.")
                sys.exit(1)
        ann[tag] = df

    tum = ann[TUMOR_CELL_TYPE]
    pa = tum[tum["effect"].isin(PROTEIN_ALTERING)].copy()
    log(f"\n  protein-altering in {TUMOR_CELL_TYPE}: {len(pa):,}")
    for k, v in pa["effect"].value_counts().items():
        log(f"    {k}: {v:,}")

    # --- subtraction ----------------------------------------------------
    sep("STEP 4: non-epithelial subtraction")
    other_keys = set(zip(ann["other"]["chrom"], ann["other"]["pos"],
                         ann["other"]["alt"]))
    keys = list(zip(pa["chrom"], pa["pos"], pa["alt"]))
    keep = [k not in other_keys for k in keys]
    log(f"  removed {len(pa) - sum(keep):,} variants also seen outside "
        f"{TUMOR_CELL_TYPE}")
    pa = pa[keep].copy()
    log(f"  somatic protein-altering: {len(pa):,}")
    pa.to_csv(f"{OUT}/{TUMOR_CELL_TYPE}.somatic_protein_altering.tsv",
              sep="\t", index=False)

    mis = pa[pa["effect"].str.contains("missense", na=False)].copy()
    log(f"  missense (MHCflurry-eligible): {len(mis):,}")
    if len(mis) == 0:
        log("\n  No missense variants. Nothing to predict. Stopping.")
        sys.exit(0)

    # --- peptides -------------------------------------------------------
    sep("STEP 5: peptide generation from the reference proteome")
    by_enst, by_symbol, by_ensg, isoforms = load_proteome(PROTEOME)

    mut_p, wt_p, metas, statuses = [], [], [], []
    diag = collections.Counter()
    for _, v in mis.iterrows():
        p = parse_hgvs_p(v["hgvs_p"])
        if p is None:
            diag["HGVS_PARSE_FAILED"] += 1
            statuses.append({"gene_symbol": v["gene"], "hgvs_p": v["hgvs_p"],
                             "status": "HGVS_PARSE_FAILED", "n_peptides": 0})
            continue
        wt, pos, mt = p
        mp, wp, md, st = generate_peptides(
            v["gene"], v["gene_id"], v["transcript_id"], wt, pos, mt,
            PEPTIDE_LENGTHS, by_enst, by_symbol, by_ensg, isoforms)
        st["chrom"], st["pos"], st["hgvs_p"] = v["chrom"], v["pos"], v["hgvs_p"]
        statuses.append(st)
        diag[st["status"]] += 1
        for m in md:
            m.update({"gene": v["gene"], "gene_id": v["gene_id"],
                      "chrom": v["chrom"], "pos": v["pos"],
                      "ref": v["ref"], "alt": v["alt"],
                      "location": f"{v['chrom']}:{v['pos']}",
                      "hgvs_p": v["hgvs_p"], "wt_aa": wt, "mut_aa": mt,
                      "mut_pos_protein": pos, "effect": v["effect"]})
        mut_p.extend(mp); wt_p.extend(wp); metas.extend(md)

    pd.DataFrame(statuses).to_csv(f"{OUT}/proteome_mapping_diagnostics.tsv",
                                  sep="\t", index=False)
    log("  mapping outcomes:")
    for k, n in diag.most_common():
        log(f"    {k}: {n:,} ({100*n/max(len(mis),1):.1f}%)")
    log(f"  peptides generated: {len(mut_p):,}")
    if not mut_p:
        sys.exit("FATAL: no peptides generated")

    # --- MHCflurry ------------------------------------------------------
    sep("STEP 6: MHCflurry binding prediction")
    try:
        from mhcflurry import Class1AffinityPredictor
        predictor = Class1AffinityPredictor.load()
    except Exception as e:
        log(f"  FATAL: MHCflurry failed to load: {e}")
        log("  pip install mhcflurry && mhcflurry-downloads fetch")
        sys.exit(1)

    alleles = []
    for a in HLA_PANEL:
        try:
            predictor.predict_to_dataframe(peptides=["GILGFVFTL"], alleles=[a])
            alleles.append(a)
        except Exception:
            log(f"  WARNING: {a} unsupported, skipping")
    log(f"  alleles: {len(alleles)}/{len(HLA_PANEL)}")
    if not alleles:
        sys.exit("FATAL: no usable alleles")

    meta_df = pd.DataFrame(metas)
    rows = []
    for a in alleles:
        log(f"  scoring {a} ...")
        mr = predictor.predict(peptides=mut_p, allele=a)
        wr = predictor.predict(peptides=wt_p, allele=a)
        d = meta_df.copy()
        d["allele"] = a
        d["mut_peptide"], d["wt_peptide"] = mut_p, wt_p
        d["mut_ic50"], d["wt_ic50"] = mr, wr
        rows.append(d)
    res = pd.concat(rows, ignore_index=True)
    res["delta_ic50"] = res["wt_ic50"] - res["mut_ic50"]
    res["is_binder"] = res["mut_ic50"] < MHC_BIND_THRESH
    res["is_strong"] = res["mut_ic50"] < STRONG_BIND_THRESH
    res["is_differential"] = res["is_binder"] & (res["wt_ic50"] > MHC_BIND_THRESH)
    res.to_csv(f"{OUT}/{TUMOR_CELL_TYPE}_all_peptide_results.tsv",
               sep="\t", index=False)

    neo = res[res["is_binder"]].copy()
    log(f"\n  peptide-allele pairs : {len(res):,}")
    log(f"  binders (<{MHC_BIND_THRESH} nM) : {len(neo):,}")
    log(f"  strong (<{STRONG_BIND_THRESH} nM) : {int(res['is_strong'].sum()):,}")
    log(f"  differential          : {int(res['is_differential'].sum()):,}")
    log(f"  unique neoantigen mutations (gene, hgvs_p): "
        f"{neo.groupby(['gene','hgvs_p']).ngroups:,}")

    # --- bead mapping ---------------------------------------------------
    sep("STEP 7: map neoantigens back to beads")
    bc = pick(beads, "#CHROM", "CHROM")
    bp = pick(beads, "Start", "POS")
    bead_idx = collections.defaultdict(list)
    for _, r in beads.iterrows():
        bead_idx[(str(r[bc]), str(int(r[bp]) - offset))].append(r["CB"])

    per_mut = []
    for (g, h), grp in neo.groupby(["gene", "hgvs_p"]):
        best = grp.loc[grp["delta_ic50"].idxmax()]
        cbs = sorted(set(bead_idx.get((str(best["chrom"]),
                                       str(int(best["pos"]) - offset)), [])))
        ref_b, alt_b = str(best["ref"]).upper(), str(best["alt"]).upper()
        tri = ""
        try:
            p = int(best["pos"]) - offset
            tri = fa.fetch(str(best["chrom"]), p - 2, p + 1).upper()
        except Exception:
            pass
        sub = f"{ref_b}>{alt_b}"
        if ref_b in ("A", "G"):
            comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
            sub = f"{comp[ref_b]}>{comp[alt_b]}"
            tri = "".join(comp[b] for b in reversed(tri)) if len(tri) == 3 else ""
        is_tcw = (len(tri) == 3 and tri[0] == "T" and tri[1] == "C"
                  and tri[2] in ("A", "T"))
        per_mut.append({
            "gene": g, "hgvs_p": h,
            "location": best["location"], "ref": best["ref"], "alt": best["alt"],
            "sub_class": sub, "tri_context": tri,
            "is_tcw": is_tcw, "is_tcw_ct": bool(is_tcw and sub == "C>T"),
            "best_allele": best["allele"],
            "mut_peptide": best["mut_peptide"], "wt_peptide": best["wt_peptide"],
            "mut_ic50": float(best["mut_ic50"]), "wt_ic50": float(best["wt_ic50"]),
            "delta_ic50": float(best["delta_ic50"]),
            "mut_position_in_peptide": int(best["mut_position_in_peptide"]),
            "is_strong": bool(best["mut_ic50"] < STRONG_BIND_THRESH),
            "is_differential": bool(best["mut_ic50"] < MHC_BIND_THRESH
                                    and best["wt_ic50"] > MHC_BIND_THRESH),
            "n_alleles_bound": int(grp["allele"].nunique()),
            "n_peptides": int(len(grp)),
            "n_beads": len(cbs), "beads": ",".join(cbs),
        })
    pm = pd.DataFrame(per_mut).sort_values(
        ["n_beads", "delta_ic50"], ascending=[False, False])
    pm.to_csv(f"{OUT}/{TUMOR_CELL_TYPE}_neoantigens_per_mutation.tsv",
              sep="\t", index=False)
    neo.to_csv(f"{OUT}/{TUMOR_CELL_TYPE}_neoantigens.tsv", sep="\t", index=False)

    log(f"  neoantigen mutations      : {len(pm):,}")
    log(f"  with single-bead support  : {int((pm['n_beads'] > 0).sum()):,}")
    log(f"  TCW context               : {int(pm['is_tcw'].sum()):,}")
    log(f"  clean TCW C>T (SBS2-like) : {int(pm['is_tcw_ct'].sum()):,}")

    log("\n  substitution class of neoantigen mutations:")
    for k, n in pm["sub_class"].value_counts().items():
        log(f"    {k}: {n:,} ({100*n/len(pm):.1f}%)")
    log("  (a large T>C share reflects the ADAR component noted in the header)")

    if (pm["n_beads"] > 0).any():
        log("\n  top 15 by bead support, then binding gain:")
        for _, r in pm[pm["n_beads"] > 0].head(15).iterrows():
            log(f"    {r['gene']:<12} {r['hgvs_p']:<16} "
                f"{r['n_beads']:>3} beads  "
                f"wt {r['wt_ic50']:>9.1f} -> mut {r['mut_ic50']:>7.1f} nM"
                f"{'  TCW' if r['is_tcw'] else ''}")

    # --- per-bead table -------------------------------------------------
    rows = []
    for _, r in pm.iterrows():
        for cb in (r["beads"].split(",") if r["beads"] else []):
            rows.append({"CB": cb, "gene": r["gene"], "hgvs_p": r["hgvs_p"],
                         "location": r["location"],
                         "mut_peptide": r["mut_peptide"],
                         "best_allele": r["best_allele"],
                         "mut_ic50": r["mut_ic50"],
                         "delta_ic50": r["delta_ic50"],
                         "is_strong": r["is_strong"],
                         "is_tcw_ct": r["is_tcw_ct"]})
    pb = pd.DataFrame(rows)
    if len(pb):
        try:
            import scanpy as sc
            ad = sc.read_h5ad(H5AD)
            for c in ("puck_id", "unified_annotation", "x_coord", "y_coord"):
                if c in ad.obs:
                    pb[c] = pb["CB"].map(ad.obs[c].astype(str) if
                                         ad.obs[c].dtype.name in ("category", "object")
                                         else ad.obs[c])
        except Exception as e:
            log(f"  WARN: could not attach spatial metadata ({e})")
        pb.to_csv(f"{OUT}/neoantigens_per_bead.tsv", sep="\t", index=False)
        log(f"\n  per-bead neoantigen table: {len(pb):,} rows, "
            f"{pb['CB'].nunique():,} beads")
        if "puck_id" in pb:
            for k, n in pb["puck_id"].value_counts().items():
                log(f"    {k}: {n:,} bead-neoantigen pairs")

    with open(f"{OUT}/step06_report.txt", "w") as f:
        f.write("\n".join(_report))

    sep(f"STEP 06 COMPLETE in {(time.time()-t0)/60:.1f} min")
    log(f"  per mutation : {OUT}/{TUMOR_CELL_TYPE}_neoantigens_per_mutation.tsv")
    log(f"  per bead     : {OUT}/neoantigens_per_bead.tsv")
    log(f"  all peptides : {OUT}/{TUMOR_CELL_TYPE}_all_peptide_results.tsv")
    log("")
    log("Next: Step07 overlap against the network paper by (gene, hgvs_p)")
    log("and by peptide, then figures.")


if __name__ == "__main__":
    main()
