#!/bin/bash

#SBATCH -J SPATIAL_05b_SCOMATIC
#SBATCH -o /master/jlehle/WORKING/LOGS/Step05b_SComatic.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step05b_SComatic.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 3-00:00:00
#SBATCH -p normal
#SBATCH --mem=900G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 80

#===============================================================================
# STEP 05b: SComatic CELL-TYPE-LEVEL VARIANT CALLING
#
# Wrapper for Step05b_SComatic_CellType_Variants.py. Runs the full SComatic
# chain on the pooled retagged BAM from Step05a:
#
#   SplitBam -> BaseCellCounter -> MergeCounts -> Calling step1 -> step2
#            -> BED filter -> callable sites -> trinucleotide background
#
# Every phase is checkpointed under SComatic/checkpoints/, so a resubmit
# resumes rather than restarting. BaseCellCounter is checkpointed PER CELL
# TYPE, which matters because epithelial is 60% of the work.
#
# Input  (from Step05a): pooled.matched.retagged.bam, meta_unified_annotation.tsv
# Output (05_mutations/SComatic/):
#   SplitBam/pooled.{cell_type}.bam
#   BaseCellCounts/*.tsv
#   MergedCounts/pooled.BaseCellCounts.AllCellTypes.tsv
#   VariantCalling/pooled.calling.step1.tsv, .step2.tsv
#   FilteredVariants/pooled.calling.filtered.tsv     <- the variant set
#   CellTypeCallableSites/, TrinucleotideBackground/
#   Step05b_summary.txt
#
# Runtime: dominated by BaseCellCounter on 136M epithelial reads.
#          Expect 6-24 h. Walltime set to 3 days.
#
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# Author:  Jake Lehle, Texas Biomedical Research Institute
# Server:  Zeus / Titan (Texas Biomed HPC)
#===============================================================================

set -o pipefail

source ~/anaconda3/bin/activate
conda activate slide-TCR-seq

#===============================================================================
# CONFIGURATION
#===============================================================================

PROJECT_ROOT="/master/jlehle/WORKING/slide-TCR-seq-working"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED"
PY_SCRIPT="${SCRIPT_DIR}/Step05b_SComatic_CellType_Variants.py"

OUTDIR="${PROJECT_ROOT}/data/outputs/05_mutations"
SC_DIR="${OUTDIR}/SComatic"

POOLED_BAM="${OUTDIR}/pooled.matched.retagged.bam"
META_FILE="${OUTDIR}/meta_unified_annotation.tsv"
GENOME_FA="${PROJECT_ROOT}/data/inputs/ref/GRCh38/GRCh38.primary_assembly.genome.fa"

SCOMATIC="/master/jlehle/WORKING/SComatic"
PON_FILE="${SCOMATIC}/PoNs/PoN.scRNAseq.hg38.tsv"
EDIT_SITES="${SCOMATIC}/RNAediting/AllEditingSites.hg38.txt"
BED_FILE="${SCOMATIC}/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed"

MIN_FREE_GB=200

#===============================================================================
# LOGGING
#===============================================================================

mkdir -p "${SC_DIR}"
LOG_FILE="${SC_DIR}/Step05b_$(date +%Y%m%d_%H%M%S).log"

log()         { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$1] ${@:2}" | tee -a "${LOG_FILE}"; }
log_info()    { log "INFO" "$@"; }
log_warn()    { log "WARN" "$@"; }
log_error()   { log "ERROR" "$@"; }
log_success() { log "SUCCESS" "$@"; }

section() {
    echo "" | tee -a "${LOG_FILE}"
    echo "========================================" | tee -a "${LOG_FILE}"
    echo "$1" | tee -a "${LOG_FILE}"
    echo "========================================" | tee -a "${LOG_FILE}"
}

#===============================================================================
# VALIDATION
#===============================================================================

validate() {
    section "VALIDATING INPUTS"
    local errors=0

    # --- tools ---
    command -v samtools &>/dev/null \
        && log_info "samtools: $(samtools --version | head -1)" \
        || { log_error "samtools not found"; ((errors++)); }
    command -v bedtools &>/dev/null \
        && log_info "bedtools: $(bedtools --version)" \
        || { log_error "bedtools not found (needed for phase 5)"; ((errors++)); }
    python -c "import pysam, pandas, scipy, numpy" 2>/dev/null \
        && log_info "python imports: OK" \
        || { log_error "python imports failed"; ((errors++)); }

    [ -f "${PY_SCRIPT}" ] && log_info "worker script: OK" \
        || { log_error "not found: ${PY_SCRIPT}"; ((errors++)); }

    # --- Step05a outputs ---
    if [ -s "${POOLED_BAM}" ]; then
        log_info "pooled BAM: OK ($(du -h ${POOLED_BAM} | cut -f1))"
        samtools quickcheck "${POOLED_BAM}" \
            && log_info "  quickcheck: OK" \
            || { log_error "  quickcheck FAILED"; ((errors++)); }
        [ -f "${POOLED_BAM}.bai" ] && log_info "  index: OK" \
            || log_warn "  no .bai; SComatic may need one"
        local first_sq=$(samtools view -H "${POOLED_BAM}" | grep -m1 '^@SQ' | sed 's/.*SN:\([^\t]*\).*/\1/')
        if [[ "${first_sq}" == chr* ]]; then
            log_info "  contigs: chr-prefixed (${first_sq})"
        else
            log_error "  contigs are '${first_sq}', expected chr-prefixed. Rerun Step05a."
            ((errors++))
        fi
    else
        log_error "not found: ${POOLED_BAM}. Run Step05a first."; ((errors++))
    fi

    if [ -s "${META_FILE}" ]; then
        local rows=$(( $(wc -l < ${META_FILE}) - 1 ))
        local hdr=$(head -1 "${META_FILE}")
        log_info "meta: OK (${rows} beads)"
        [ "${hdr}" == "$(printf 'Index\tCell_type')" ] \
            && log_info "  header: Index/Cell_type OK" \
            || { log_error "  header is '${hdr}', expected Index<TAB>Cell_type"; ((errors++)); }
        log_info "  cell types: $(tail -n +2 ${META_FILE} | cut -f2 | sort -u | tr '\n' ' ')"
    else
        log_error "not found: ${META_FILE}"; ((errors++))
    fi

    # --- reference and aux ---
    [ -f "${GENOME_FA}" ]     && log_info "genome FASTA: OK" || { log_error "not found: ${GENOME_FA}"; ((errors++)); }
    [ -f "${GENOME_FA}.fai" ] && log_info "FASTA index: OK"  || { log_error "not found: ${GENOME_FA}.fai (BaseCellCounter needs it)"; ((errors++)); }
    for f in "${PON_FILE}" "${EDIT_SITES}" "${BED_FILE}"; do
        [ -f "${f}" ] && log_info "aux OK: $(basename ${f})" \
            || { log_error "aux missing: ${f}"; ((errors++)); }
    done
    for f in SplitBam/SplitBamCellTypes.py \
             BaseCellCounter/BaseCellCounter.py \
             MergeCounts/MergeBaseCellCounts.py \
             BaseCellCalling/BaseCellCalling.step1.py \
             BaseCellCalling/BaseCellCalling.step2.py \
             GetCallableSites/GetAllCallableSites.py \
             TrinucleotideBackground/TrinucleotideContextBackground.py; do
        [ -f "${SCOMATIC}/scripts/${f}" ] || { log_error "SComatic script missing: ${f}"; ((errors++)); }
    done
    [ ${errors} -eq 0 ] && log_info "all SComatic scripts present"

    # --- disk ---
    local free_gb=$(df -BG "${OUTDIR}" | tail -1 | awk '{gsub("G","",$4); print $4}')
    [ "${free_gb}" -ge "${MIN_FREE_GB}" ] \
        && log_info "free space: ${free_gb}G (need ~${MIN_FREE_GB}G)" \
        || { log_error "free space ${free_gb}G below ${MIN_FREE_GB}G"; ((errors++)); }

    # --- resume state ---
    if [ -d "${SC_DIR}/checkpoints" ]; then
        local n=$(ls "${SC_DIR}/checkpoints"/*.done 2>/dev/null | wc -l)
        [ "${n}" -gt 0 ] && log_info "resuming: ${n} checkpoint(s) already set"
    fi

    [ ${errors} -gt 0 ] && { log_error "validation failed with ${errors} error(s)"; return 1; }
    log_success "all inputs validated"
    return 0
}

#===============================================================================
# MAIN
#===============================================================================

if [ "${1:-}" == "--help" ] || [ "${1:-}" == "-h" ]; then
    cat << EOF
Step 05b: SComatic cell-type-level variant calling

Usage: sbatch $(basename $0)
       $(basename $0) --validate-only   # run checks, do not call variants
       $(basename $0) --clean           # clear checkpoints, start fresh
       $(basename $0) --help

Runs SplitBam through BED filtering on the pooled retagged BAM from
Step05a, using meta_unified_annotation.tsv (9 cell types, 99,341 beads).

Every phase is checkpointed; BaseCellCounter is checkpointed per cell
type. A resubmit after a timeout resumes where it stopped.

Parameters match the network paper so the resulting variant set is
directly comparable to its Tier 1/2/3 neoantigen lists. Edit the
CONFIGURATION block of the .py to change them.
EOF
    exit 0
fi

section "STEP 05b: SComatic CELL-TYPE VARIANT CALLING"
log_info "Started at $(date) on $(hostname)"
log_info "Cores: ${SLURM_CPUS_PER_TASK:-80}   Mem: ${SLURM_MEM_PER_NODE:-900G}"

if [ "${1:-}" == "--clean" ]; then
    log_warn "clearing all checkpoints"
    rm -f "${SC_DIR}/checkpoints"/*.done
fi

validate || exit 1

if [ "${1:-}" == "--validate-only" ]; then
    log_success "validate-only requested; stopping here"
    exit 0
fi

section "RUNNING SComatic"
log_info "BaseCellCounter on 136M epithelial reads is the long pole; expect 6-24 h"

python "${PY_SCRIPT}" 2>&1 | tee -a "${LOG_FILE}"
RC=${PIPESTATUS[0]}

if [ ${RC} -ne 0 ]; then
    log_error "Step05b failed with exit code ${RC}"
    log_error "Checkpoints are preserved; resubmit to resume from the failed phase."
    exit 1
fi

#===============================================================================
# POST-RUN VERIFICATION
#===============================================================================

section "VERIFYING OUTPUT"

FILTERED="${SC_DIR}/FilteredVariants/pooled.calling.filtered.tsv"
for f in "${SC_DIR}/VariantCalling/pooled.calling.step1.tsv" \
         "${SC_DIR}/VariantCalling/pooled.calling.step2.tsv" \
         "${FILTERED}"; do
    [ -s "${f}" ] && log_info "OK: $(basename ${f}) ($(du -h ${f} | cut -f1))" \
        || { log_error "missing or empty: ${f}"; exit 1; }
done

log_info "PASS variants: $(grep -vc '^#' ${FILTERED})"

section "SUMMARY"
[ -f "${SC_DIR}/Step05b_summary.txt" ] && cat "${SC_DIR}/Step05b_summary.txt" | tee -a "${LOG_FILE}"

section "STEP 05b COMPLETE"
log_success "Finished at $(date)"
log_info ""
log_info "Next: SnpEff annotation and neoantigen prediction (Step06), then"
log_info "overlap against the network paper Tier 1/2/3 lists by gene + AA change."

exit 0
