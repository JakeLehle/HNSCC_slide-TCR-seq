#!/bin/bash

#SBATCH -J SPATIAL_07_OVERLAP
#SBATCH -o /master/jlehle/WORKING/LOGS/Step07_Overlap.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step07_Overlap.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 0-04:00:00
#SBATCH -p normal
#SBATCH --mem=200G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16

#===============================================================================
# STEP 07: OVERLAP WITH THE NETWORK PAPER + SPATIAL FIGURES
#
# Wrapper for Step07_Overlap_And_Spatial_Figures.py.
#
# Compares the 411 spatial neoantigen mutations against the network paper's
# 775 at four levels (locus, mutation, peptide, gene), then draws the spatial
# figures. Fast: minutes.
#
# Input:
#   07_neoantigen/epithelial_neoantigens_per_mutation.tsv
#   07_neoantigen/neoantigens_per_bead.tsv
#   2026_NMF_PAPER/data/FIG_7/06_prevalence_ranking/neoantigen_prevalence_ranking_full.tsv
#
# Output (08_overlap_figures/):
#   shared_neoantigen_mutations.tsv, shared_gene.tsv, shared_peptide.tsv
#   gene_overlap_stats.tsv
#   neoantigen_beads_by_puck_gene.tsv
#   step07_report.txt
#   figures/  Fig_A .. Fig_E, PDF + PNG at 300 DPI, individual panels
#
# The network paper tree must be readable from this node. If it is not,
# the script falls back to the per-group binder sets and says so.
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
PY_SCRIPT="${SCRIPT_DIR}/Step07_Overlap_And_Spatial_Figures.py"

NEO="${PROJECT_ROOT}/data/outputs/07_neoantigen"
OUT="${PROJECT_ROOT}/data/outputs/08_overlap_figures"
H5AD="${PROJECT_ROOT}/data/outputs/04_annotation/all_pucks_annotated_unified.h5ad"

NMF="/master/jlehle/WORKING/2026_NMF_PAPER/data/FIG_7"
NMF_RANKING="${NMF}/06_prevalence_ranking/neoantigen_prevalence_ranking_full.tsv"

#===============================================================================
# LOGGING
#===============================================================================

mkdir -p "${OUT}"
LOG_FILE="${OUT}/Step07_$(date +%Y%m%d_%H%M%S).log"

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

    python -c "import scanpy, pandas, numpy, scipy, matplotlib" 2>/dev/null \
        && log_info "python imports: OK" \
        || { log_error "python imports failed"; ((errors++)); }
    [ -f "${PY_SCRIPT}" ] && log_info "worker script: OK" \
        || { log_error "not found: ${PY_SCRIPT}"; ((errors++)); }

    # --- Step06 outputs ---
    for f in "${NEO}/epithelial_neoantigens_per_mutation.tsv" \
             "${NEO}/neoantigens_per_bead.tsv"; do
        [ -s "${f}" ] \
            && log_info "OK: $(basename ${f}) ($(( $(wc -l < ${f}) - 1 )) rows)" \
            || { log_error "not found or empty: ${f}. Run Step06."; ((errors++)); }
    done
    [ -f "${H5AD}" ] && log_info "annotated h5ad: OK" \
        || { log_error "not found: ${H5AD}"; ((errors++)); }

    # --- network paper, readable from this node? ---
    if [ -r "${NMF_RANKING}" ]; then
        log_info "network ranking: OK ($(( $(wc -l < ${NMF_RANKING}) - 1 )) mutations, expect 775)"
        log_info "  tier counts:"
        awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="tier") t=i; next}
                    t{c[$t]++} END{for(k in c) printf "    %s: %d\n", k, c[k]}' \
            "${NMF_RANKING}" | tee -a "${LOG_FILE}"
    else
        log_warn "network ranking not readable: ${NMF_RANKING}"
        log_warn "  the script will fall back to the per-group binder sets"
        for g in SBS2_HIGH CNV_HIGH; do
            f="${NMF}/03_mhc_binding/${g}_neoantigens.tsv"
            [ -r "${f}" ] && log_info "  fallback OK: ${g}" \
                || log_warn "  fallback MISSING: ${f}"
        done
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
Step 07: Overlap with the network paper + spatial figures

Usage: sbatch $(basename $0)
       $(basename $0) --validate-only
       $(basename $0) --help

Overlap is computed at four levels: locus, mutation (gene + hgvs_p),
peptide, and gene. Only the gene level carries a hypergeometric test;
the others have no clean null and are reported as counts and names.

Figures (individual PDF + PNG, 300 DPI, for Illustrator):
  Fig_A  cell types per puck                       [v2: 2x point size]
  Fig_B  mutation-carrying and neoantigen beads in space
  Fig_C  compartment size and detection rate by cell type
  Fig_D  top candidates by bead support, TCW flagged
  Fig_E  overlap with the network paper
  Fig_F  NEW: neoantigen-expressing beads only, large red points,
         hotspot rings (solid = beads sharing one mutation)
  Fig_G  NEW: multi-bead neoantigens, carriers of one mutation joined
  Fig_H  NEW: permutation test, are neoantigen beads clustered

New tables: neoantigen_clonality.tsv, neoantigen_hotspots.tsv,
spatial_clustering_stats.tsv
EOF
    exit 0
fi

section "STEP 07: OVERLAP AND SPATIAL FIGURES"
log_info "Started at $(date) on $(hostname)"

validate || exit 1

if [ "${1:-}" == "--validate-only" ]; then
    log_success "validate-only requested; stopping here"
    exit 0
fi

section "RUNNING"
python "${PY_SCRIPT}" 2>&1 | tee -a "${LOG_FILE}"
RC=${PIPESTATUS[0]}

if [ ${RC} -ne 0 ]; then
    log_error "Step07 failed with exit code ${RC}"
    exit 1
fi

section "VERIFYING OUTPUT"
n_fig=$(ls "${OUT}/figures"/*.pdf 2>/dev/null | wc -l)
log_info "figures written: ${n_fig} PDF"
[ "${n_fig}" -ge 7 ] || { log_error "expected at least 7 figures, got ${n_fig}"; exit 1; }
for f in "${OUT}/figures"/*.pdf; do
    log_info "  $(basename ${f}) ($(du -h ${f} | cut -f1))"
done
[ -s "${OUT}/step07_report.txt" ] && log_info "report: OK" \
    || { log_error "report missing"; exit 1; }

for f in neoantigen_clonality.tsv neoantigen_hotspots.tsv \
         spatial_clustering_stats.tsv; do
    [ -s "${OUT}/${f}" ] \
        && log_info "  ${f}: $(( $(wc -l < ${OUT}/${f}) - 1 )) rows" \
        || log_warn "  ${f}: not written (nothing met the criteria)"
done

section "SUMMARY"
sed -n '/STEP 3: overlap/,/STEP 5: figures/p' "${OUT}/step07_report.txt" \
    | tee -a "${LOG_FILE}"

section "STEP 07 COMPLETE"
log_success "Finished at $(date)"

exit 0
