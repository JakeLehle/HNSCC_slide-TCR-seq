#!/bin/bash

#SBATCH -J SPATIAL_05d_SIGS
#SBATCH -o /master/jlehle/WORKING/LOGS/Step05d_Signatures.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step05d_Signatures.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 0-04:00:00
#SBATCH -p normal
#SBATCH --mem=200G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16

#===============================================================================
# STEP 05d: SEMI-SUPERVISED COSMIC SIGNATURE REFITTING
#
# Wrapper for Step05d_Signature_Refitting.py. Fits COSMIC signatures per bead
# using the network paper's approach: a fixed core set (SBS2, SBS13, SBS5)
# plus HNSCC-associated candidates admitted by scree-plot elbow detection.
#
# Fast: ~430 beads x 96 contexts. Minutes, not hours.
#
# Input  (from Step05c):
#   SingleCell/FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv
#   SingleCell/CombinedCallableSites/complete_callable_sites.tsv   (optional)
# Output (06_signatures/):
#   signature_weights_per_cell.txt           raw NNLS, signatures x beads
#   signature_weights_per_cell_relative.txt  normalized per bead
#   signature_weights_per_celltype.txt       pseudobulk companion
#   mutation_matrix_96contexts.txt
#   per_bead_fit_quality.tsv
#   Step05d_scree.tsv, Step05d_summary.txt
#   all_pucks_annotated_signatures.h5ad
#   figures/  (PDF + PNG at 300 DPI)
#
# NOTE: signature_weights_per_cell.txt is written signatures-as-rows to match
# the network paper convention. It requires .T on load.
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
PY_SCRIPT="${SCRIPT_DIR}/Step05d_Signature_Refitting.py"

OUTDIR="${PROJECT_ROOT}/data/outputs/05_mutations"
OUT_SC="${OUTDIR}/SComatic/SingleCell"
SIGDIR="${PROJECT_ROOT}/data/outputs/06_signatures"

MUTATIONS="${OUT_SC}/FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv"
CALLABLE="${OUT_SC}/CombinedCallableSites/complete_callable_sites.tsv"
H5AD="${PROJECT_ROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"

#===============================================================================
# LOGGING
#===============================================================================

mkdir -p "${SIGDIR}"
LOG_FILE="${SIGDIR}/Step05d_$(date +%Y%m%d_%H%M%S).log"

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
# COSMIC LOCATOR
#
# COSMIC_FILE in the .py is a best guess. This surfaces what is actually on
# disk so the path can be corrected before the run rather than after.
#===============================================================================

COSMIC_FILE="/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/COSMIC_v3.4_SBS_GRCh38.txt"

find_cosmic() {
    section "CHECKING COSMIC SBS MATRIX"
    if [ -r "${COSMIC_FILE}" ]; then
        log_info "COSMIC: ${COSMIC_FILE} ($(du -h ${COSMIC_FILE} | cut -f1))"
        log_info "  rows: $(( $(wc -l < ${COSMIC_FILE}) - 1 )) (expect 96)"
        log_info "  signatures: $(( $(head -1 ${COSMIC_FILE} | tr '\t' '\n' | wc -l) - 1 ))"
        log_info "  first context: $(sed -n 2p ${COSMIC_FILE} | cut -f1) (expect A[C>A]A)"
        # v3.4 splits SBS40; the .py expands the requested SBS40 accordingly.
        local sbs40=$(head -1 "${COSMIC_FILE}" | tr '\t' '\n' | grep -c '^SBS40')
        log_info "  SBS40 variants present: ${sbs40}"
        for s in SBS1 SBS2 SBS4 SBS5 SBS7a SBS7b SBS13 SBS16 SBS17a SBS17b \
                 SBS18 SBS29 SBS39 SBS44; do
            head -1 "${COSMIC_FILE}" | tr '\t' '\n' | grep -qx "${s}" \
                || log_warn "  requested signature ${s} NOT in this file"
        done
        return 0
    fi

    log_warn "not readable: ${COSMIC_FILE}"
    log_warn "searching for alternatives ..."
    for pat in \
        "/master/jlehle/WORKING/SC/fastq/Head_and_neck_cancer/*SBS*" \
        "/master/jlehle/WORKING/SC/ref/COSMIC/*SBS*" \
        "/master/jlehle/WORKING/*/COSMIC*SBS*" \
        "/master/jlehle/WORKING/*/*/COSMIC*SBS*"; do
        for f in ${pat}; do
            [ -f "${f}" ] && log_info "candidate: ${f}"
        done
    done
    log_warn "Set COSMIC_FILE at the top of the .py before running."
    log_warn "Expected: 96 trinucleotide contexts as rows, SBS columns."
}

#===============================================================================
# VALIDATION
#===============================================================================

validate() {
    section "VALIDATING INPUTS"
    local errors=0

    python -c "import scanpy, pandas, numpy, scipy, matplotlib" 2>/dev/null \
        && log_info "python imports: OK" \
        || { log_error "python imports failed"; ((errors++)); }
    [ -f "${PY_SCRIPT}" ] && log_info "worker script: OK" \
        || { log_error "not found: ${PY_SCRIPT}"; ((errors++)); }

    if [ -s "${MUTATIONS}" ]; then
        local n=$(( $(wc -l < ${MUTATIONS}) - 1 ))
        log_info "mutations: OK (${n} bead-variant rows)"
        local hdr=$(head -1 "${MUTATIONS}")
        for col in CB REF_TRI ALT_TRI; do
            echo "${hdr}" | tr '\t' '\n' | grep -qx "${col}" \
                && log_info "  column ${col}: present" \
                || { log_error "  column ${col} MISSING"; ((errors++)); }
        done
        local beads=$(tail -n +2 "${MUTATIONS}" | cut -f$(echo "${hdr}" | tr '\t' '\n' | grep -nx CB | cut -d: -f1) | sort -u | wc -l)
        log_info "  distinct beads: ${beads}"
    else
        log_error "not found or empty: ${MUTATIONS}. Run Step05c first."; ((errors++))
    fi

    [ -s "${CALLABLE}" ] \
        && log_info "callable sites: OK" \
        || log_warn "callable sites absent; the mask step will be skipped"

    [ -f "${H5AD}" ] && log_info "annotated h5ad: OK" \
        || { log_error "not found: ${H5AD}"; ((errors++)); }

    [ ${errors} -gt 0 ] && { log_error "validation failed with ${errors} error(s)"; return 1; }
    log_success "all inputs validated"
    return 0
}

#===============================================================================
# MAIN
#===============================================================================

if [ "${1:-}" == "--help" ] || [ "${1:-}" == "-h" ]; then
    cat << EOF
Step 05d: Semi-supervised COSMIC signature refitting

Usage: sbatch $(basename $0)
       $(basename $0) --validate-only   # checks + COSMIC locator
       $(basename $0) --help

Core signatures  : SBS2, SBS13, SBS5
HNSCC candidates : SBS1, SBS4, SBS7a, SBS7b, SBS16, SBS17a, SBS17b,
                   SBS18, SBS29, SBS39, SBS40, SBS44
Selection        : scree elbow (second derivative, L-method cross-check)

Edit the CONFIGURATION block of the .py for COSMIC_FILE, MUT_THRESHOLD,
or ELBOW_METHOD.

The per-bead fit is reported stratified by mutation count. Beads with a
single mutation are degenerate by construction; the per-cell-type
pseudobulk companion is the stable comparison.
EOF
    exit 0
fi

section "STEP 05d: SIGNATURE REFITTING"
log_info "Started at $(date) on $(hostname)"

find_cosmic
validate || exit 1

if [ "${1:-}" == "--validate-only" ]; then
    log_success "validate-only requested; stopping here"
    exit 0
fi

section "RUNNING"

python "${PY_SCRIPT}" 2>&1 | tee -a "${LOG_FILE}"
RC=${PIPESTATUS[0]}

if [ ${RC} -ne 0 ]; then
    log_error "Step05d failed with exit code ${RC}"
    exit 1
fi

section "VERIFYING OUTPUT"
for f in "${SIGDIR}/signature_weights_per_cell.txt" \
         "${SIGDIR}/signature_weights_per_cell_relative.txt" \
         "${SIGDIR}/Step05d_summary.txt" \
         "${SIGDIR}/all_pucks_annotated_signatures.h5ad"; do
    [ -s "${f}" ] && log_info "OK: $(basename ${f}) ($(du -h ${f} | cut -f1))" \
        || { log_error "missing or empty: ${f}"; exit 1; }
done
log_info "figures: $(ls ${SIGDIR}/figures/*.pdf 2>/dev/null | wc -l) PDF"

section "SUMMARY"
cat "${SIGDIR}/Step05d_summary.txt" | tee -a "${LOG_FILE}"

section "STEP 05d COMPLETE"
log_success "Finished at $(date)"
log_info ""
log_info "Next: SnpEff annotation of the epithelial variants, then MHCflurry"
log_info "against MHC-I supertype representatives, then overlap with the"
log_info "network paper Tier 1/2/3 lists by gene + AA substitution."

exit 0
