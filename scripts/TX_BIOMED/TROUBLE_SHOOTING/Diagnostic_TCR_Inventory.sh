#!/bin/bash
# =============================================================================
# Diagnostic_TCR_Inventory.sh
# -----------------------------------------------------------------------------
# READ-ONLY survey of the Slide-TCR-seq input tree to scope the clonotype layer.
# Writes NOTHING into the data directories. All output goes to stdout; redirect
# to a log and paste it back:
#
#   source ~/anaconda3/bin/activate sc_pre   # need conda samtools (system one broken)
#   bash Diagnostic_TCR_Inventory.sh 2>&1 | tee TCR_inventory.log
#
# It answers three questions:
#   (1) Is there a pre-made TCR / clonotype product anywhere on disk?
#   (2) Which flowcell (H52J2DMXY vs HLGH2BGXK) carries the TCR-enriched library?
#       -> a TCR library piles up at TRA (chr14) and TRB (chr7); transcriptome
#          is spread genome-wide.
#   (3) Do the BAM reads carry recoverable bead-barcode + UMI tags (XC/XB/XM)?
#
# NOTE: the BAM survey uses an UNBIASED random subsample (samtools -s), which
# reads through each BAM once. Expect a couple of minutes per large BAM. Lower
# SAMPLE_FRAC if it drags; raise it if a flowcell's BAMs are small.
# =============================================================================
set -uo pipefail

# ---- EDIT IF NEEDED ----------------------------------------------------------
INPUT_ROOT="/master/jlehle/WORKING/slide-TCR-seq-working/data/inputs"
FASTQ_DIR="${INPUT_ROOT}/fastq"
SAMPLE_FRAC="0.005"   # fraction of reads to randomly sample per BAM
# -----------------------------------------------------------------------------

PUCKS=(Puck_211214_29 Puck_211214_37 Puck_211214_40)
RUNDIRS=(H52J2DMXY HLGH2BGXK)
LANES=(L001 L002 L003 L004)

echo "samtools: $(command -v samtools)"
samtools --version 2>/dev/null | head -1
echo "INPUT_ROOT: ${INPUT_ROOT}"
echo

# -----------------------------------------------------------------------------
# 1. Any pre-made TCR / clonotype product anywhere?
# -----------------------------------------------------------------------------
echo "=============================================================="
echo " 1. Search for TCR / clonotype / VDJ-named files (depth <= 8)"
echo "=============================================================="
find "${FASTQ_DIR}" -maxdepth 8 -type f \
  \( -iname '*tcr*'  -o -iname '*clonotype*' -o -iname '*cdr3*' -o -iname '*vdj*' \
  -o -iname '*trac*' -o -iname '*trbc*'      -o -iname '*mixcr*' -o -iname '*trust4*' \
  -o -iname '*repertoire*' \) 2>/dev/null | sort
echo "  ( no lines above = no pre-made TCR product on disk )"
echo

# -----------------------------------------------------------------------------
# 2. What is inside the per-lane alignment/ subdirs (the depth-6 blind spot)?
# -----------------------------------------------------------------------------
echo "=============================================================="
echo " 2. Contents of one alignment/ dir per flowcell"
echo "=============================================================="
for rd in "${RUNDIRS[@]}"; do
  d=$(find "${FASTQ_DIR}/2022-01-28_${PUCKS[0]}/${rd}" -maxdepth 2 -type d -name alignment 2>/dev/null | head -1)
  echo "--- ${rd}: ${d:-<none>}"
  [ -n "${d}" ] && ls -lhR "${d}" 2>/dev/null | head -50
  echo
done

# -----------------------------------------------------------------------------
# 3. Read structure from the one raw run folder we have (HLGH2BGXK)
# -----------------------------------------------------------------------------
echo "=============================================================="
echo " 3. RunInfo.xml read structure"
echo "=============================================================="
RUNINFO=$(find "${FASTQ_DIR}" -name RunInfo.xml 2>/dev/null | head -1)
echo "RunInfo: ${RUNINFO:-<none>}"
[ -n "${RUNINFO}" ] && grep -Eo '<Read [^>]*/>' "${RUNINFO}" | sed 's/^/  /'
echo

# -----------------------------------------------------------------------------
# 4. BAM content survey
# -----------------------------------------------------------------------------
survey_bam () {
  local bam="$1" label="$2"
  if [ ! -f "${bam}" ]; then echo "  MISSING: ${bam}"; return; fi
  echo "----- ${label}"
  echo "  path : ${bam}  ($(du -h "${bam}" 2>/dev/null | cut -f1))"
  if samtools quickcheck "${bam}" 2>/dev/null; then echo "  check: OK"; else echo "  check: FAIL (truncated/corrupt)"; return; fi

  local tmp; tmp=$(mktemp)
  samtools view -s "${SAMPLE_FRAC}" "${bam}" 2>/dev/null > "${tmp}"
  local n; n=$(wc -l < "${tmp}")
  echo "  sampled reads (frac ${SAMPLE_FRAC}): ${n}"
  if [ "${n}" -eq 0 ]; then echo "  (empty sample; raise SAMPLE_FRAC)"; rm -f "${tmp}"; return; fi

  echo "  top R2 read lengths (col10):"
  awk '{print length($10)}' "${tmp}" | sort -n | uniq -c | sort -rn | head -4 | sed 's/^/    /'

  echo "  reads carrying each barcode/UMI tag:"
  for tag in XC XB XM CB UB; do
    printf "    %s: %s\n" "${tag}" "$(grep -c "${tag}:Z:" "${tmp}" 2>/dev/null)"
  done

  echo "  top reference contigs in sample:"
  awk '{print $3}' "${tmp}" | sort | uniq -c | sort -rn | head -8 | sed 's/^/    /'

  echo "  fraction landing in TCR loci:"
  awk 'BEGIN{tra=0;trb=0;trg=0;tot=0}
       {tot++; c=$3; sub(/^chr/,"",c); p=$4+0;
        if(c=="14" && p>=21500000 && p<=22600000) tra++;
        else if(c=="7" && p>=142200000 && p<=142900000) trb++;
        else if(c=="7" && p>=38200000  && p<=38400000)  trg++;}
       END{if(tot>0) printf "    TRA(chr14)=%.3f%%  TRB(chr7)=%.3f%%  TRG(chr7)=%.3f%%  (n=%d)\n",
            100*tra/tot,100*trb/tot,100*trg/tot,tot}' "${tmp}"
  rm -f "${tmp}"
}

echo "=============================================================="
echo " 4. final.bam survey, per flowcell x puck x lane"
echo "    (high TRA/TRB fraction => this is the TCR-enriched flowcell)"
echo "=============================================================="
for puck in "${PUCKS[@]}"; do
  for rd in "${RUNDIRS[@]}"; do
    for lane in "${LANES[@]}"; do
      bam="${FASTQ_DIR}/2022-01-28_${puck}/${rd}/${lane}/${puck}.final.bam"
      [ -f "${bam}" ] && survey_bam "${bam}" "${puck} / ${rd} / ${lane} / final"
    done
  done
done
echo

# -----------------------------------------------------------------------------
# 5. matched.bam (transcriptome) for barcode-format parity reference
# -----------------------------------------------------------------------------
echo "=============================================================="
echo " 5. matched.bam (transcriptome) tag/format reference"
echo "=============================================================="
for puck in "${PUCKS[@]}"; do
  survey_bam "${FASTQ_DIR}/2022-01-28_${puck}/${puck}.matched.bam" "${puck} / matched"
done
echo

echo "=============================================================="
echo " HOW TO READ THIS"
echo "  - Sec 1 empty + Sec 2 only .bam/.bai  -> clonotypes are a from-BAM build."
echo "  - The flowcell with TRA+TRB >> a few %  is the TCR library."
echo "  - Whichever BAM shows XB (or XC) + XM tags is one we can put a bead barcode"
echo "    and UMI on; that decides whether TRUST4 (-b BAM --barcode --UMI) runs"
echo "    straight off it or we extract reads first."
echo "=============================================================="
echo "DONE"
