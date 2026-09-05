#!/bin/bash
#SBATCH -J TCR_FULLSWEEP
#SBATCH -o /master/jlehle/WORKING/LOGS/TCR_fullsweep.%j.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/TCR_fullsweep.%j.e.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=END,FAIL
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH --mem=400G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
# (%j in the log names = job id, so an overnight rerun never clobbers a prior log)
# =============================================================================
# Diagnostic_TCR_FullSweep.sh
# -----------------------------------------------------------------------------
# EXHAUSTIVE, READ-ONLY sweep of the Slide-TCR-seq inputs tree. Goal: be
# absolutely certain there is (or is not) any rhTCRseq / VDJ / CDR3 / long-read
# (Nanopore) / MiSeq data hiding anywhere before asking Sophia for more.
#
# Writes NOTHING into the data dirs. All output -> stdout. Run as:
#   source ~/anaconda3/bin/activate sc_pre
#   bash Diagnostic_TCR_FullSweep.sh 2>&1 | tee TCR_fullsweep.log
#
# Phases run cheap -> expensive. The metadata/keyword/feature phases finish fast
# and usually settle the question; the per-BAM read-length scan (Phase 5) reads
# every BAM end-to-end (tens of minutes) and is the definitive long-read test.
# Safe to Ctrl-C after Phase 4 if you only want the metadata verdict.
# =============================================================================
set -uo pipefail

# batch jobs do not inherit the interactive shell env: activate sc_pre here so
# the conda samtools (system samtools is broken) is on PATH.
source ~/anaconda3/bin/activate
conda activate sc_pre

THREADS="${SLURM_CPUS_PER_TASK:-8}"   # for samtools BGZF decompression

# ---- EDIT IF NEEDED ----------------------------------------------------------
INPUT_ROOT="/master/jlehle/WORKING/slide-TCR-seq-working/data/inputs"
FASTQ_DIR="${INPUT_ROOT}/fastq"
# -----------------------------------------------------------------------------

# keywords that would betray a TCR / long-read / targeted library or product
KW='rhtcr|[^a-z]tcr[^a-z]|tcr-|tcrseq|vdj|cdr3|clonotype|repertoire|nanopore|minion|promethion|oxford|pacbio|miseq|mixcr|trust4|trav|trbv|traj|trbj|igblast'

echo "samtools: $(command -v samtools)"; samtools --version 2>/dev/null | head -1
echo "pdftotext: $(command -v pdftotext || echo 'NOT INSTALLED (will strings PDFs)')"
echo "INPUT_ROOT: ${INPUT_ROOT}"
echo "date: $(date)"
echo

# =============================================================================
echo "##############################################################"
echo "# PHASE 1  File census (every file accounted for, by type)"
echo "##############################################################"
echo "--- total files / dirs:"
find "${INPUT_ROOT}" -type f | wc -l | sed 's/^/  files: /'
find "${INPUT_ROOT}" -type d | wc -l | sed 's/^/  dirs : /'

echo "--- count by extension:"
find "${INPUT_ROOT}" -type f | sed -E 's/.*\.([A-Za-z0-9]+)$/\1/' | sort | uniq -c | sort -rn | sed 's/^/  /'

echo "--- 'file' type on anything NOT a known transcriptome artifact"
echo "    (looking for fastq / fast5 / pod5 / unexpected formats):"
find "${INPUT_ROOT}" -type f \
  ! -iname '*.bam' ! -iname '*.bai' ! -iname '*.txt' ! -iname '*.txt.gz' \
  ! -iname '*.tsv' ! -iname '*.tsv.gz' ! -iname '*.mtx.gz' ! -iname '*.pdf' \
  ! -iname '*.xml' ! -iname '*.bin' ! -iname '*.cfg' ! -iname '*.zip' \
  ! -iname '*.gtf' ! -iname '*.gtf.gz' ! -iname '*.fa' ! -iname '*.fa.gz' \
  ! -iname '*.tab' ! -iname '*.out' ! -iname '*.pickle' ! -iname '*.summary.txt' \
  ! -iname '*.cfg' ! -iname '*.xml' 2>/dev/null \
  | while read -r f; do printf "  %s : %s\n" "$f" "$(file -b "$f")"; done
echo "  (no lines above = every file is a known type)"

echo "--- explicit hunt for read-data file formats anywhere:"
find "${INPUT_ROOT}" -type f \
  \( -iname '*.fastq*' -o -iname '*.fq*' -o -iname '*.fast5' -o -iname '*.pod5' \
  -o -iname '*.blow5' -o -iname '*sequencing_summary*' -o -iname '*fastq_pass*' \) 2>/dev/null \
  | sed 's/^/  /'
echo "  (no lines = no raw fastq / no Nanopore signal files present)"
echo

# =============================================================================
echo "##############################################################"
echo "# PHASE 2  Keyword sweep (metadata, logs, pickles, PDFs, zips)"
echo "##############################################################"
echo "--- small text/metadata files containing TCR/long-read keywords:"
find "${INPUT_ROOT}" -type f \
  \( -iname '*.txt' -o -iname '*.summary.txt' -o -iname '*.xml' -o -iname '*.tab' \
  -o -iname '*.out' -o -iname '*.cfg' -o -iname '*.csv' -o -iname 'reference_info*' \
  -o -iname 'star_version*' \) \
  ! -iname '*digital_expression*' ! -iname '*barcode_distribution*' \
  ! -iname '*numReads*' 2>/dev/null \
  | while read -r f; do
      m=$(grep -iEo "${KW}" "$f" 2>/dev/null | sort -u | tr '\n' ' ')
      [ -n "$m" ] && printf "  [HIT] %s : %s\n" "$f" "$m"
    done
echo "  (no [HIT] lines = no keywords in any metadata file)"

echo "--- pickles (strings):"
find "${INPUT_ROOT}" -type f -iname '*.pickle' 2>/dev/null | while read -r f; do
  m=$(strings -n 5 "$f" 2>/dev/null | grep -iEo "${KW}" | sort -u | tr '\n' ' ')
  [ -n "$m" ] && printf "  [HIT] %s : %s\n" "$f" "$m"
done
echo "  (no [HIT] = pickles clean)"

echo "--- PDFs:"
find "${INPUT_ROOT}" -type f -iname '*.pdf' 2>/dev/null | while read -r f; do
  if command -v pdftotext >/dev/null 2>&1; then txt=$(pdftotext "$f" - 2>/dev/null); else txt=$(strings "$f" 2>/dev/null); fi
  m=$(printf '%s' "$txt" | grep -iEo "${KW}" | sort -u | tr '\n' ' ')
  [ -n "$m" ] && printf "  [HIT] %s : %s\n" "$f" "$m"
done
echo "  (no [HIT] = PDFs mention nothing TCR/long-read)"

echo "--- inside Logs.zip (listing + keyword grep):"
find "${INPUT_ROOT}" -type f -iname '*.zip' 2>/dev/null | while read -r z; do
  echo "  zip: $z"
  unzip -l "$z" 2>/dev/null | sed 's/^/    /'
  m=$(unzip -p "$z" 2>/dev/null | grep -iEo "${KW}" | sort -u | tr '\n' ' ')
  [ -n "$m" ] && printf "    [HIT] keywords inside zip: %s\n" "$m"
done
echo

# =============================================================================
echo "##############################################################"
echo "# PHASE 3  TCR feature control (V/J should be absent, C present)"
echo "##############################################################"
find "${FASTQ_DIR}" -type f -iname '*matched.digital_expression_features.tsv.gz' 2>/dev/null | while read -r feat; do
  echo "--- $feat"
  echo "  constant-region genes present (expected YES):"
  zcat "$feat" 2>/dev/null | grep -iE 'TRAC|TRBC|TRDC|TRGC' | sed 's/^/    /' | head
  nC=$(zcat "$feat" 2>/dev/null | grep -icE 'TRAC|TRBC|TRDC|TRGC')
  nV=$(zcat "$feat" 2>/dev/null | grep -icE 'TR[ABGD]V[0-9]')
  nJ=$(zcat "$feat" 2>/dev/null | grep -icE 'TR[ABGD]J[0-9]')
  printf "  counts -> constant(C)=%s  variable(V)=%s  joining(J)=%s\n" "$nC" "$nV" "$nJ"
  echo "  (paper predicts: C present, V/J near zero. Lots of V/J would be a surprise worth chasing.)"
done
echo

# =============================================================================
echo "##############################################################"
echo "# PHASE 4  Run-folder metadata (library names, chemistry, sheet)"
echo "##############################################################"
echo "--- RunInfo.xml read structure:"
find "${INPUT_ROOT}" -name RunInfo.xml 2>/dev/null | while read -r r; do
  echo "  $r"; grep -Eo '<Read [^>]*/>' "$r" | sed 's/^/    /'
done
echo "--- RunParameters.xml (library / experiment / chemistry fields):"
find "${INPUT_ROOT}" -iname 'RunParameters.xml' 2>/dev/null | while read -r r; do
  echo "  $r"
  grep -iE 'library|experiment|chemistry|<Read|ExperimentName|FlowCell|Application' "$r" | sed 's/^/    /' | head -30
done
echo "--- any SampleSheet anywhere:"
find "${INPUT_ROOT}" -iname '*samplesheet*' 2>/dev/null | sed 's/^/  /'
echo "  (no SampleSheet lines = none delivered)"
echo

# =============================================================================
echo "##############################################################"
echo "# PHASE 5  Exhaustive per-BAM read-length scan  (LONG-READ CATCHER)"
echo "#   one full pass per BAM: read count, max len, tail bins, tags"
echo "#   ANY reads >70bp => a non-transcriptome (MiSeq/Nanopore) population"
echo "##############################################################"
find "${FASTQ_DIR}" -type f -iname '*.bam' 2>/dev/null | sort | while read -r bam; do
  printf -- "----- %s  (%s)\n" "$bam" "$(du -h "$bam" 2>/dev/null | cut -f1)"
  if ! samtools quickcheck "$bam" 2>/dev/null; then echo "  check: FAIL"; continue; fi
  samtools view -@ "${THREADS}" "$bam" 2>/dev/null | awk '
    { n++
      L=length($10); if(L>maxL)maxL=L
      if(L>70)b70++; if(L>100)b100++; if(L>150)b150++; if(L>300)b300++; if(L>1000)b1000++
      if($0 ~ /XC:Z:/) xc++; if($0 ~ /XB:Z:/) xb++; if($0 ~ /XM:Z:/) xm++ }
    END {
      printf "  reads=%d  maxlen=%d\n", n, maxL
      printf "  len>70=%d  >100=%d  >150=%d  >300=%d  >1000=%d\n", b70+0,b100+0,b150+0,b300+0,b1000+0
      printf "  tags: XC=%d  XB=%d  XM=%d\n", xc+0,xb+0,xm+0 }'
done
echo

echo "##############################################################"
echo "# HOW TO READ THIS"
echo "#  - Phase 1: any fastq/fast5/pod5/unknown file = uninventoried read data."
echo "#  - Phase 2: any [HIT] = a TCR/long-read trace left in a log/pdf/zip."
echo "#  - Phase 3: V/J genes absent + C present == matches the paper's chemistry"
echo "#    (fragmentation keeps only the constant region; CDR3 lives in the"
echo "#    targeted fraction, which is what we'd be missing)."
echo "#  - Phase 4: a TCR library named in RunParameters/SampleSheet would mean"
echo "#    it was sequenced on THIS run after all."
echo "#  - Phase 5: maxlen ~60 and all >70 bins == 0 across EVERY bam => there is"
echo "#    no MiSeq/Nanopore VDJ population anywhere in this data. That is the"
echo "#    'absolutely sure' result that justifies going back to Sophia."
echo "##############################################################"
echo "DONE  $(date)"
