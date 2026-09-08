#!/bin/bash

#SBATCH -J SPATIAL_05c_SCGENO
#SBATCH -o /master/jlehle/WORKING/LOGS/Step05c_SCGenotype.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step05c_SCGenotype.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -p normal
#SBATCH --mem=900G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 80

#===============================================================================
# STEP 05c: PER-BEAD GENOTYPING
#
# Wrapper for Step05c_SingleCellGenotype.py. Assigns the Step05b variant set
# to individual beads, then annotates each call with trinucleotide context
# taken FROM THE GENOME (SComatic's context fields are genome-correct in only
# 23 of 355 cases per the network paper).
#
# Produces the table signature_analysis.py consumes, and the per-bead
# mutation map the neoantigen and spatial colocalization work needs.
#
# Input  (from Step05b): SplitBam/pooled.{cell_type}.bam
#                        FilteredVariants/pooled.calling.filtered.tsv
#                        meta_unified_annotation.tsv
# Output (05_mutations/SComatic/SingleCell/):
#   FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv
#   mutations_per_bead.tsv
#   Step05c_feasibility.txt          <- READ THIS, it decides the next step
#   CombinedCallableSites/complete_callable_sites.tsv
#
# Runtime: phases 1-3 fast (2,431 sites). Phase 4 (SitesPerCell) is the long
#          pole because it iterates the 3.8 GB step1 file once per cell type;
#          set RUN_SITES_PER_CELL=False in the .py to skip it.
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
PY_SCRIPT="${SCRIPT_DIR}/Step05c_SingleCellGenotype.py"

OUTDIR="${PROJECT_ROOT}/data/outputs/05_mutations"
SC_DIR="${OUTDIR}/SComatic"
OUT_SC="${SC_DIR}/SingleCell"

SPLIT_DIR="${SC_DIR}/SplitBam"
FILTERED="${SC_DIR}/FilteredVariants/pooled.calling.filtered.tsv"
META_FILE="${OUTDIR}/meta_unified_annotation.tsv"
GENOME_FA="${PROJECT_ROOT}/data/inputs/ref/GRCh38/GRCh38.primary_assembly.genome.fa"

SCOMATIC="/master/jlehle/WORKING/SComatic"

#===============================================================================
# LOGGING
#===============================================================================

mkdir -p "${OUT_SC}"
LOG_FILE="${OUT_SC}/Step05c_$(date +%Y%m%d_%H%M%S).log"

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

    python -c "import pysam, pandas, numpy" 2>/dev/null \
        && log_info "python imports: OK" \
        || { log_error "python imports failed"; ((errors++)); }
    [ -f "${PY_SCRIPT}" ] && log_info "worker script: OK" \
        || { log_error "not found: ${PY_SCRIPT}"; ((errors++)); }

    # --- Step05b outputs ---
    local n_bam=$(ls ${SPLIT_DIR}/pooled.*.bam 2>/dev/null | wc -l)
    if [ "${n_bam}" -gt 0 ]; then
        log_info "split BAMs: ${n_bam} found"
        for b in ${SPLIT_DIR}/pooled.*.bam; do
            [ -f "${b}.bai" ] || log_warn "  no index for $(basename ${b})"
        done
    else
        log_error "no split BAMs in ${SPLIT_DIR}. Run Step05b first."; ((errors++))
    fi

    if [ -s "${FILTERED}" ]; then
        log_info "variant file: OK ($(grep -vc '^#' ${FILTERED}) PASS variants)"
    else
        log_error "not found or empty: ${FILTERED}"; ((errors++))
    fi

    [ -s "${META_FILE}" ] \
        && log_info "meta: OK ($(( $(wc -l < ${META_FILE}) - 1 )) beads)" \
        || { log_error "not found: ${META_FILE}"; ((errors++)); }

    [ -f "${GENOME_FA}" ]     && log_info "genome FASTA: OK" || { log_error "not found: ${GENOME_FA}"; ((errors++)); }
    [ -f "${GENOME_FA}.fai" ] && log_info "FASTA index: OK"  || { log_error "not found: ${GENOME_FA}.fai (needed for context lookup)"; ((errors++)); }

    for f in SingleCellGenotype/SingleCellGenotype.py \
             SitesPerCell/SitesPerCell.py; do
        [ -f "${SCOMATIC}/scripts/${f}" ] && log_info "SComatic: $(basename ${f}) OK" \
            || { log_error "SComatic script missing: ${f}"; ((errors++)); }
    done

    local free_gb=$(df -BG "${OUTDIR}" | tail -1 | awk '{gsub("G","",$4); print $4}')
    log_info "free space: ${free_gb}G"

    if [ -d "${OUT_SC}/checkpoints" ]; then
        local n=$(ls "${OUT_SC}/checkpoints"/*.done 2>/dev/null | wc -l)
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
Step 05c: Per-bead genotyping

Usage: sbatch $(basename $0)
       $(basename $0) --validate-only
       $(basename $0) --clean          # clear checkpoints
       $(basename $0) --help

Runs SingleCellGenotype.py on each cell-type BAM against the 2,431 PASS
variants from Step05b, filters per-bead calls, and adds trinucleotide
context from GRCh38 (not from SComatic's context fields).

Two things this does NOT do, deliberately:
  - It does not append a suffix to CB. CB is already the h5ad obs_name.
  - It does not trust SComatic's Up_context / Down_context fields.

Read Step05c_feasibility.txt afterwards. It reports how many beads carry
enough mutations for a per-bead signature fit to mean anything.
EOF
    exit 0
fi

section "STEP 05c: PER-BEAD GENOTYPING"
log_info "Started at $(date) on $(hostname)"
log_info "Cores: ${SLURM_CPUS_PER_TASK:-80}"

if [ "${1:-}" == "--clean" ]; then
    log_warn "clearing all checkpoints"
    rm -f "${OUT_SC}/checkpoints"/*.done
fi

validate || exit 1

if [ "${1:-}" == "--validate-only" ]; then
    log_success "validate-only requested; stopping here"
    exit 0
fi

section "RUNNING"

python "${PY_SCRIPT}" 2>&1 | tee -a "${LOG_FILE}"
RC=${PIPESTATUS[0]}

if [ ${RC} -ne 0 ]; then
    log_error "Step05c failed with exit code ${RC}"
    log_error "Checkpoints preserved; resubmit to resume."
    exit 1
fi

#===============================================================================
# POST-RUN
#===============================================================================

section "VERIFYING OUTPUT"

COMBINED="${OUT_SC}/FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv"
for f in "${COMBINED}" "${OUT_SC}/mutations_per_bead.tsv" \
         "${OUT_SC}/Step05c_feasibility.txt"; do
    [ -s "${f}" ] && log_info "OK: $(basename ${f}) ($(du -h ${f} | cut -f1))" \
        || { log_error "missing or empty: ${f}"; exit 1; }
done

# The columns signature_analysis.py hard-requires.
HDR=$(head -1 "${COMBINED}")
for col in REF ALT_expected REF_TRI ALT_TRI CB; do
    echo "${HDR}" | tr '\t' '\n' | grep -qx "${col}" \
        && log_info "column ${col}: present" \
        || { log_error "column ${col} MISSING from ${COMBINED}"; exit 1; }
done

section "FEASIBILITY READOUT"
cat "${OUT_SC}/Step05c_feasibility.txt" | tee -a "${LOG_FILE}"

section "STEP 05c COMPLETE"
log_success "Finished at $(date)"
log_info ""
log_info "Pick --mutation-threshold from the table above, then run"
log_info "signature_analysis.py with --hnscc-only --use-scree."

exit 0
