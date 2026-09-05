#!/bin/bash
#SBATCH -J MANIFEST_AUDIT
#SBATCH -o /master/jlehle/WORKING/LOGS/audit.%j.o.log
#SBATCH -e /master/jlehle/WORKING/LOGS/audit.%j.e.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=END,FAIL
#SBATCH -p normal
#SBATCH -t 00:30:00
#SBATCH --mem=16G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
# =============================================================================
# Diagnostic_Manifest_Audit_v2.sh
# -----------------------------------------------------------------------------
# Verifies every line of the inventory manifest against the DELIVERED tree
# (data/inputs). v2 fixes the v1 problems:
#   - v1 crawled $HOME + other projects -> 12h hang + false FASTQ/SampleSheet
#     FAILs from unrelated data. v2 scopes everything to data/inputs.
#   - v1 caught the wrong dir as the run folder (name glob). v2 finds it by
#     RunInfo.xml.
#   - v1 ran `samtools view -c` on a corrupt 8G unmapped.bam. v2 checks the
#     28-byte BGZF EOF marker instead (instant).
# No full-file reads in the default path -> finishes in minutes.
# =============================================================================
set -uo pipefail
source ~/anaconda3/bin/activate
conda activate sc_pre

# ---- roots / paths -----------------------------------------------------------
PROJECT="/master/jlehle/WORKING/slide-TCR-seq-working"
INPUTS="${PROJECT}/data/inputs"          # the delivered tree = the audit universe
FASTQ="${INPUTS}/fastq"
H5AD="${PROJECT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"
STEP01="$(find ${PROJECT}/scripts -name 'Step01_Alignment.sh' 2>/dev/null | head -1)"
PUCKS=(Puck_211214_29 Puck_211214_37 Puck_211214_40)
VERIFY_READ_COUNTS="${VERIFY_READ_COUNTS:-0}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ---- claimed values (from the inventory doc) --------------------------------
CLAIM_N_PUCKS=3
CLAIM_FINAL_H52=2; CLAIM_FINAL_HLGH=4
CLAIM_UNMAPPED_TOTAL=18
CLAIM_FASTQ=0; CLAIM_SAMPLESHEET=0; CLAIM_H52_RAW=0
CLAIM_CB_LEN=14; CLAIM_UMI_LEN=9
CLAIM_NOBS=99341; CLAIM_TCELL=5184
declare -A CLAIM_MATCHED_READS=( [Puck_211214_29]=121668905 [Puck_211214_37]=94572075 [Puck_211214_40]=141079737 )

PASS=0; FAIL=0
verdict(){ local v; if [ "$2" == "$3" ]; then v="PASS"; ((PASS++)); else v="FAIL"; ((FAIL++)); fi
  printf "  [%s] %-46s claim=%-10s observed=%-10s\n" "$v" "$1" "$2" "$3"; }
note(){ printf "  [....] %s\n" "$*"; }

echo "date: $(date)"; echo "samtools: $(command -v samtools)"

# ---- ONE filesystem walk, cached; everything greps this ---------------------
LIST="$(mktemp)"; trap 'rm -f "$LIST"' EXIT
t0=$SECONDS
find "${INPUTS}" -type f > "${LIST}" 2>/dev/null
echo "cached file list: $(wc -l < "${LIST}") files in $((SECONDS-t0))s"; echo

# =============================================================================
echo "############### A. STRUCTURE ###############"
n=$(for p in "${PUCKS[@]}"; do [ -d "${FASTQ}/2022-01-28_${p}" ] && echo x; done | wc -l)
verdict "puck directories present" "${CLAIM_N_PUCKS}" "${n}"
note "inputs total files: $(wc -l < "${LIST}")"
for sub in fastq panglaodb ref; do
  note "$(printf '%8s  %s' "$(grep -c "/inputs/${sub}/" "${LIST}")" "inputs/${sub}/")"
done
echo

# =============================================================================
echo "############### B. RAW DATA MATRIX (delivered tree only) ###############"
# run folder = the dir that actually contains a RunInfo.xml (not a name glob)
mapfile -t RUNINFOS < <(grep '/RunInfo.xml$' "${LIST}")
note "RunInfo.xml files (= sequencing run folders): ${#RUNINFOS[@]}"
h52raw=0
for ri in "${RUNINFOS[@]}"; do
  rf=$(dirname "$ri")
  fc=$(grep -Eo 'Flowcell="[^"]+"|<Flowcell>[^<]+' "$ri" 2>/dev/null | head -1 | sed -E 's/.*[">]//')
  note "  run folder: ${rf}   flowcell=${fc:-?}"
  echo "$fc" | grep -qi 'H52J2DMXY' && ((h52raw++))
  nbcl=$(grep -icE "^${rf}/.*/BaseCalls/.*\.(bcl|cbcl)" "${LIST}")
  note "    BaseCalls bcl/cbcl files: ${nbcl}   RTAComplete=$( [ -f "${rf}/RTAComplete.txt" ] && echo Y || echo N )"
  [ -f "$ri" ] && { grep -Eo '<Read [^>]*/>' "$ri" | sed 's/^/      read: /'; }
done
verdict "H52J2DMXY RAW run folders (expect none)" "${CLAIM_H52_RAW}" "${h52raw}"
verdict "FASTQ in delivered tree (expect none)"   "${CLAIM_FASTQ}"   "$(grep -icE '\.(fastq|fq)(\.gz)?$' "${LIST}")"
verdict "SampleSheet in delivered tree (expect 0)" "${CLAIM_SAMPLESHEET}" "$(grep -ic 'samplesheet' "${LIST}")"
verdict "Nanopore signal files (expect none)" "0" "$(grep -icE '\.(fast5|pod5)$|sequencing_summary' "${LIST}")"
echo

# =============================================================================
echo "############### C. ALIGNED TRANSCRIPTOME BAMS (per puck) ###############"
for p in "${PUCKS[@]}"; do
  echo "--- ${p}"
  verdict "  H52J2DMXY final.bam lanes" "${CLAIM_FINAL_H52}"  "$(grep -c "2022-01-28_${p}/H52J2DMXY/.*\.final\.bam$" "${LIST}")"
  verdict "  HLGH2BGXK final.bam lanes" "${CLAIM_FINAL_HLGH}" "$(grep -c "2022-01-28_${p}/HLGH2BGXK/.*\.final\.bam$" "${LIST}")"
  m="${FASTQ}/2022-01-28_${p}/${p}.matched.bam"; a="${FASTQ}/2022-01-28_${p}/${p}.all_illumina.bam"
  note "  matched.bam:      $( [ -f "$m" ] && echo "present $(du -h "$m"|cut -f1) quickcheck=$(samtools quickcheck "$m" 2>/dev/null && echo OK || echo FAIL)" || echo MISSING )"
  note "  all_illumina.bam: $( [ -f "$a" ] && echo "present $(du -h "$a"|cut -f1)" || echo MISSING )"
done
# barcode/UMI structure from a small sample (head stops the stream immediately)
mb="${FASTQ}/2022-01-28_${PUCKS[0]}/${PUCKS[0]}.matched.bam"
if [ -f "$mb" ]; then
  echo "--- barcode/UMI (sampled 5000 reads from ${PUCKS[0]} matched.bam)"
  smp="$(samtools view "$mb" 2>/dev/null | head -5000)"
  for tag in XC XB XM; do note "  reads with ${tag}: $(printf '%s\n' "$smp" | grep -c "${tag}:Z:") / 5000"; done
  cbl=$(printf '%s\n' "$smp" | grep -oE 'XC:Z:[ACGTN]+' | head -1 | sed 's/XC:Z://' | tr -d '\n' | wc -c)
  uml=$(printf '%s\n' "$smp" | grep -oE 'XM:Z:[ACGTN]+' | head -1 | sed 's/XM:Z://' | tr -d '\n' | wc -c)
  verdict "  cell-barcode (XC) length" "${CLAIM_CB_LEN}"  "${cbl}"
  verdict "  UMI (XM) length"          "${CLAIM_UMI_LEN}" "${uml}"
fi
echo

# =============================================================================
echo "############### D. unmapped.bam corruption claim ###############"
mapfile -t UNMAPPED < <(grep '\.unmapped\.bam$' "${LIST}")
verdict "unmapped.bam files total" "${CLAIM_UNMAPPED_TOTAL}" "${#UNMAPPED[@]}"
nfail=0; for u in "${UNMAPPED[@]}"; do samtools quickcheck "$u" 2>/dev/null || ((nfail++)); done
verdict "unmapped.bam FAILING quickcheck" "${CLAIM_UNMAPPED_TOTAL}" "${nfail}"
if [ "${#UNMAPPED[@]}" -gt 0 ]; then
  u="${UNMAPPED[0]}"
  BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
  got=$(tail -c 28 "$u" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')
  if [ "$got" == "$BGZF_EOF" ]; then note "one unmapped has intact BGZF EOF -> quickcheck fussy, not truncated"
  else note "one unmapped MISSING BGZF EOF marker -> genuinely truncated transfer"; fi
fi
echo

# =============================================================================
echo "############### E. DGE / SPATIAL MAPS / REF ###############"
for p in "${PUCKS[@]}"; do
  d="${FASTQ}/2022-01-28_${p}"; echo "--- ${p}"
  mtx="${d}/${p}.matched.digital_expression_matrix.mtx.gz"
  bcs="${d}/${p}.matched.digital_expression_barcodes.tsv.gz"
  fts="${d}/${p}.matched.digital_expression_features.tsv.gz"
  [ -f "$mtx" ] && note "  mtx dims (genes beads nnz): $(zcat "$mtx" 2>/dev/null | grep -m1 -vE '^%')" || note "  mtx MISSING"
  [ -f "$bcs" ] && note "  barcodes (beads): $(zcat "$bcs" 2>/dev/null | wc -l)"
  [ -f "$fts" ] && note "  features (genes): $(zcat "$fts" 2>/dev/null | wc -l)"
  bm="${d}/barcode_matching"
  note "  spatial: BeadBarcodes=$( [ -f "${bm}/BeadBarcodes.txt" ] && echo Y||echo N) BeadLocations=$( [ -f "${bm}/BeadLocations.txt" ] && echo Y||echo N) matching.gz=$(grep -c "2022-01-28_${p}/barcode_matching/.*barcode_matching.txt.gz$" "${LIST}") xy.gz=$(grep -c "2022-01-28_${p}/barcode_matching/.*barcode_xy.txt.gz$" "${LIST}")"
done
note "ref fasta: $(grep -m1 '/ref/GRCh38/.*\.fa$' "${LIST}" || echo MISSING)"
note "ref gtf:   $(grep -m1 '/ref/GRCh38/.*\.gtf$' "${LIST}" || echo MISSING)"
note "STAR idx:  SA=$(grep -c '/star/SA$' "${LIST}") SAindex=$(grep -c '/star/SAindex$' "${LIST}") Genome=$(grep -c '/star/Genome$' "${LIST}")"
note "PanglaoDB: $(grep -m1 '/panglaodb/PanglaoDB_markers_' "${LIST}" || echo MISSING)"
echo

# =============================================================================
echo "############### F. ANNOTATED OBJECT ###############"
python - "$H5AD" "$CLAIM_NOBS" "$CLAIM_TCELL" << 'PY'
import sys, os
h5, cnobs, ctcell = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
try:
    import anndata as ad
except Exception as e:
    print("  [....] anndata import failed:", e); sys.exit(0)
if not os.path.exists(h5):
    print("  [FAIL] annotated object MISSING:", h5); sys.exit(0)
a = ad.read_h5ad(h5, backed='r')
nobs, nvar = a.shape
print("  path:", h5)
print("  [%s] n_obs=%d (claim %d)" % ("PASS" if nobs==cnobs else "FAIL", nobs, cnobs))
print("  n_var=%d" % nvar)
col = "unified_annotation"
if col in a.obs:
    vc = a.obs[col].value_counts()
    for k, v in vc.items(): print("      %-14s %d" % (k, int(v)))
    t = int(vc.get("T_cell", 0))
    print("  [%s] T_cell=%d (claim %d)" % ("PASS" if t==ctcell else "FAIL", t, ctcell))
else:
    print("  [FAIL] no", col, "in obs. cols:", list(a.obs.columns)[:25])
PY
echo

# =============================================================================
echo "############### G. REPO-DERIVED VALUES (Step01) ###############"
if [ -n "${STEP01}" ]; then
  note "Step01: ${STEP01}"
  grep -E 'Puck_211214_(29|37|40),' "${STEP01}" 2>/dev/null | sed 's/^/        /'
  grep -iE 'CB_LEN|UMI_LEN' "${STEP01}" 2>/dev/null | grep -vE '^\s*#' | head -4 | sed 's/^/        /'
else note "Step01 not found under ${PROJECT}/scripts"; fi
echo

# =============================================================================
if [ "${VERIFY_READ_COUNTS}" == "1" ]; then
echo "############### H. matched.bam READ COUNTS (slow, full pass, threaded) ###############"
for p in "${PUCKS[@]}"; do
  m="${FASTQ}/2022-01-28_${p}/${p}.matched.bam"; [ -f "$m" ] || { note "${p}: MISSING"; continue; }
  verdict "${p} matched.bam reads" "${CLAIM_MATCHED_READS[$p]}" "$(samtools view -c -@ "${THREADS}" "$m" 2>/dev/null)"
done
echo
else echo "(read-count verification skipped; set VERIFY_READ_COUNTS=1 for a threaded full-pass count)"; echo; fi

echo "############### SUMMARY ###############"
echo "  PASS=${PASS}  FAIL=${FAIL}   elapsed=$((SECONDS))s"
echo "  Any FAIL = manifest wrong on that line; paste the log and I'll correct the doc."
echo "DONE  $(date)"
