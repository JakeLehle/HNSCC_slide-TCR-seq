#!/bin/bash

#SBATCH -J SPATIAL_05a_RETAG
#SBATCH -o /master/jlehle/WORKING/LOGS/Step05a_Retag.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step05a_Retag.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 1-00:00:00
#SBATCH -p normal
#SBATCH --mem=200G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16

#===============================================================================
# STEP 05a: RETAG matched.bam FOR SComatic
#
# Wrapper for Step05a_Retag_MatchedBAM_For_SComatic.py. Validates every input
# and every SComatic auxiliary file BEFORE launching the streaming pass, so a
# missing PoN is caught in seconds rather than after a 2-hour retag.
#
# Produces:
#   05_mutations/pooled.matched.retagged.bam(.bai)
#   05_mutations/meta_unified_annotation.tsv        <- primary
#   05_mutations/meta_consensus_annotation.tsv      <- alternative
#   05_mutations/celltype_read_counts.tsv           <- READ THIS BEFORE 05b
#   05_mutations/retag_stats.tsv
#   05_mutations/Step05a_summary.txt
#
# Runtime: ~1-3 h (streaming pass over ~357M reads, 3 processes)
# Disk:    ~25-35 GB peak (3 temp BAMs + merged output)
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
SCRIPT_DIR="${PROJECT_ROOT}/scripts/TX_BIOMED"
PY_SCRIPT="${SCRIPT_DIR}/Step05a_Retag_MatchedBAM_For_SComatic.py"

INPUTS="${PROJECT_ROOT}/data/inputs/fastq"
H5AD="${PROJECT_ROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"
OUTDIR="${PROJECT_ROOT}/data/outputs/05_mutations"
GENOME_FA="${PROJECT_ROOT}/data/inputs/ref/GRCh38/GRCh38.primary_assembly.genome.fa"

# --- SComatic install and auxiliary files (shared with the network paper) ---
SCOMATIC="/master/jlehle/WORKING/SComatic"
PON_FILE="${SCOMATIC}/PoNs/PoN.scRNAseq.hg38.tsv"
EDIT_SITES="${SCOMATIC}/RNAediting/AllEditingSites.hg38.txt"
BED_FILE="${SCOMATIC}/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed"

SAMPLES=("Puck_211214_29" "Puck_211214_37" "Puck_211214_40")
SAMPLE_DIRS=("2022-01-28_Puck_211214_29" "2022-01-28_Puck_211214_37" "2022-01-28_Puck_211214_40")

MIN_FREE_GB=60

#===============================================================================
# LOGGING
#===============================================================================

mkdir -p "${OUTDIR}"
LOG_FILE="${OUTDIR}/Step05a_$(date +%Y%m%d_%H%M%S).log"

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
    python -c "import pysam, scanpy, pandas" 2>/dev/null \
        && log_info "python: pysam $(python -c 'import pysam;print(pysam.__version__)'), scanpy OK" \
        || { log_error "python imports failed (pysam / scanpy / pandas)"; ((errors++)); }

    # --- our files ---
    [ -f "${PY_SCRIPT}" ] && log_info "worker script: OK" \
        || { log_error "not found: ${PY_SCRIPT}"; ((errors++)); }
    [ -f "${H5AD}" ] \
        && log_info "annotated h5ad: OK ($(du -h ${H5AD} | cut -f1))" \
        || { log_error "not found: ${H5AD}"; ((errors++)); }
    [ -f "${GENOME_FA}" ] && log_info "genome FASTA: OK" \
        || { log_error "not found: ${GENOME_FA}"; ((errors++)); }
    [ -f "${GENOME_FA}.fai" ] && log_info "genome FASTA index: OK" \
        || { log_warn "no .fai for ${GENOME_FA}; BaseCellCounter will need one"; }

    # --- Sophia's BAMs ---
    for i in "${!SAMPLES[@]}"; do
        local mb="${INPUTS}/${SAMPLE_DIRS[$i]}/${SAMPLES[$i]}.matched.bam"
        if [ -f "${mb}" ]; then
            log_info "matched.bam ${SAMPLES[$i]}: OK ($(du -h ${mb} | cut -f1))"
            samtools quickcheck "${mb}" 2>/dev/null \
                || { log_error "  quickcheck FAILED for ${mb}"; ((errors++)); }
        else
            log_error "not found: ${mb}"; ((errors++))
        fi
    done

    # --- SComatic install and aux files, checked NOW not after the retag ---
    [ -d "${SCOMATIC}/scripts" ] && log_info "SComatic scripts: OK" \
        || { log_error "not found: ${SCOMATIC}/scripts"; ((errors++)); }
    for f in "${SCOMATIC}/scripts/SplitBam/SplitBamCellTypes.py" \
             "${SCOMATIC}/scripts/BaseCellCounter/BaseCellCounter.py" \
             "${SCOMATIC}/scripts/MergeCounts/MergeBaseCellCounts.py" \
             "${SCOMATIC}/scripts/BaseCellCalling/BaseCellCalling.step1.py" \
             "${SCOMATIC}/scripts/BaseCellCalling/BaseCellCalling.step2.py"; do
        [ -f "${f}" ] || { log_error "SComatic script missing: ${f}"; ((errors++)); }
    done
    for f in "${PON_FILE}" "${EDIT_SITES}" "${BED_FILE}"; do
        [ -f "${f}" ] && log_info "aux OK: $(basename ${f})" \
            || { log_error "aux missing: ${f}"; ((errors++)); }
    done

    # Aux files must be chr-prefixed to match the relabelled BAM.
    if [ -f "${BED_FILE}" ]; then
        local first=$(head -1 "${BED_FILE}" | cut -f1)
        if [[ "${first}" == chr* ]]; then
            log_info "BED contig style: chr-prefixed (matches relabelled BAM)"
        else
            log_error "BED first contig is '${first}', expected chr-prefixed"
            ((errors++))
        fi
    fi

    # --- disk ---
    local free_gb=$(df -BG "${OUTDIR}" | tail -1 | awk '{gsub("G","",$4); print $4}')
    if [ "${free_gb}" -ge "${MIN_FREE_GB}" ]; then
        log_info "free space: ${free_gb}G (need ~${MIN_FREE_GB}G)"
    else
        log_error "free space ${free_gb}G below ${MIN_FREE_GB}G"; ((errors++))
    fi

    [ ${errors} -gt 0 ] && { log_error "validation failed with ${errors} error(s)"; return 1; }
    log_success "all inputs validated"
    return 0
}

#===============================================================================
# MAIN
#===============================================================================

section "STEP 05a: RETAG matched.bam FOR SComatic"
log_info "Started at $(date) on $(hostname)"
log_info "Project root: ${PROJECT_ROOT}"
log_info "Output:       ${OUTDIR}"
log_info "Threads:      ${SLURM_CPUS_PER_TASK:-16}"

if [ "${1:-}" == "--help" ] || [ "${1:-}" == "-h" ]; then
    cat << EOF
Step 05a: Retag Sophia Liu's matched.bam files for SComatic

Usage: sbatch $(basename $0)
       $(basename $0) --validate-only   # run checks, do not retag
       $(basename $0) --help

Builds one pooled BAM whose CB tag equals the h5ad obs_name, so the
SComatic meta Index matches the BAM with no lookup table. Applies three
transforms: contig relabel (Ensembl -> GENCODE), nM copied from NM, and
CB set from XB plus the puck id.

Edit the CONFIGURATION block of the .py for filter thresholds,
annotation column, or PER_PUCK mode.
EOF
    exit 0
fi

validate || exit 1

if [ "${1:-}" == "--validate-only" ]; then
    log_success "validate-only requested; stopping here"
    exit 0
fi

section "RUNNING RETAG"
log_info "Streaming ~357M reads across 3 pucks; expect 1-3 hours"

python "${PY_SCRIPT}" 2>&1 | tee -a "${LOG_FILE}"
RC=${PIPESTATUS[0]}

if [ ${RC} -ne 0 ]; then
    log_error "Step05a failed with exit code ${RC}"
    exit 1
fi

#===============================================================================
# POST-RUN VERIFICATION
#===============================================================================

section "VERIFYING OUTPUT"

POOLED="${OUTDIR}/pooled.matched.retagged.bam"
for f in "${POOLED}" "${POOLED}.bai" \
         "${OUTDIR}/meta_unified_annotation.tsv" \
         "${OUTDIR}/celltype_read_counts.tsv"; do
    [ -s "${f}" ] && log_info "OK: $(basename ${f}) ($(du -h ${f} | cut -f1))" \
        || { log_error "missing or empty: ${f}"; exit 1; }
done

samtools quickcheck "${POOLED}" \
    && log_info "quickcheck: OK" \
    || { log_error "quickcheck FAILED on ${POOLED}"; exit 1; }

log_info "Pooled BAM reads: $(samtools view -c -@ ${SLURM_CPUS_PER_TASK:-16} ${POOLED})"
log_info "Meta file rows:   $(( $(wc -l < ${OUTDIR}/meta_unified_annotation.tsv) - 1 ))"

section "READS PER CELL TYPE (primary column)"
awk -F'\t' 'NR==1 || $1=="unified_annotation"' \
    "${OUTDIR}/celltype_read_counts.tsv" | column -t | tee -a "${LOG_FILE}"

section "STEP 05a COMPLETE"
log_success "Finished at $(date)"
log_info ""
log_info "Review celltype_read_counts.tsv before Step05b. It shows whether each"
log_info "cell type carries the depth to clear min_cov 5 / min_cells 5."
log_info ""
log_info "Next: SplitBamCellTypes.py, then BaseCellCounter (Step05b)."

exit 0
