#!/bin/bash

#SBATCH -J SPATIAL_06_NEOANTIGEN
#SBATCH -o /master/jlehle/WORKING/LOGS/Step06_Neoantigen.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step06_Neoantigen.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 0-12:00:00
#SBATCH -p normal
#SBATCH --mem=300G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16

#===============================================================================
# STEP 06: NEOANTIGEN PREDICTION (spatial)
#
# Wrapper for Step06_Neoantigen_Prediction.py. Runs in the NEOANTIGEN env,
# which carries SnpEff and MHCflurry. The network paper split this across
# three scripts with an env switch; the spatial variant set is small enough
# that one script in one env is simpler and has fewer seams.
#
# Input:
#   05_mutations/SComatic/FilteredVariants/pooled.calling.filtered.tsv
#   05_mutations/SComatic/SingleCell/FilteredSingleCellAlleles/
#       all_cell.single_cell_genotype.filtered.tsv
# Output (07_neoantigen/):
#   epithelial_neoantigens_per_mutation.tsv   <- the ranked catalog
#   neoantigens_per_bead.tsv                  <- which bead carries what
#   epithelial_neoantigens.tsv                <- binder rows, network-paper shape
#   epithelial_all_peptide_results.tsv        <- full audit
#   proteome_mapping_diagnostics.tsv
#   step06_report.txt
#
# TWO THINGS THAT MUST NOT DRIFT
#   1. SComatic Start is 1-based. VCF POS = Start, no +1. The script probes
#      the genome and aborts if the convention disagrees, and aborts again if
#      SnpEff emits WARNING_REF_DOES_NOT_MATCH_GENOME on more than 5% of
#      records, which is what a coordinate shift looks like.
#   2. The HLA panel is the network paper's 10 alleles. Widening it breaks
#      the overlap comparison in Step07.
#
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# Author:  Jake Lehle, Texas Biomedical Research Institute
# Server:  Zeus / Titan (Texas Biomed HPC)
#===============================================================================

set -o pipefail

source ~/anaconda3/bin/activate
conda activate NEOANTIGEN

#===============================================================================
# CONFIGURATION
#===============================================================================

PROJECT_ROOT="/master/jlehle/WORKING/slide-TCR-seq-working"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED"
PY_SCRIPT="${SCRIPT_DIR}/Step06_Neoantigen_Prediction.py"

SC_DIR="${PROJECT_ROOT}/data/outputs/05_mutations/SComatic"
OUT="${PROJECT_ROOT}/data/outputs/07_neoantigen"

CATALOG="${SC_DIR}/FilteredVariants/pooled.calling.filtered.tsv"
BEAD_MAP="${SC_DIR}/SingleCell/FilteredSingleCellAlleles/all_cell.single_cell_genotype.filtered.tsv"
GENOME="${PROJECT_ROOT}/data/inputs/ref/GRCh38/GRCh38.primary_assembly.genome.fa"
PROTEOME="/master/jlehle/WORKING/2026_NMF_PAPER/data/reference/Homo_sapiens.GRCh38.pep.all.fa"

#===============================================================================
# LOGGING
#===============================================================================

mkdir -p "${OUT}"
LOG_FILE="${OUT}/Step06_$(date +%Y%m%d_%H%M%S).log"

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

    # --- env ---
    log_info "conda env: ${CONDA_DEFAULT_ENV:-unknown}"
    command -v snpEff &>/dev/null \
        && log_info "SnpEff: $(snpEff -version 2>&1 | head -1)" \
        || { log_error "snpEff not on PATH (need the NEOANTIGEN env)"; ((errors++)); }
    python -c "import mhcflurry" 2>/dev/null \
        && log_info "mhcflurry: $(python -c 'import mhcflurry;print(mhcflurry.__version__)')" \
        || { log_error "mhcflurry import failed"; ((errors++)); }
    python -c "import pysam, pandas, numpy" 2>/dev/null \
        && log_info "pysam / pandas / numpy: OK" \
        || { log_error "python imports failed"; ((errors++)); }

    # MHCflurry model weights are a separate download from the package.
    python - <<'PY' 2>/dev/null || { log_error "MHCflurry models not downloaded: run 'mhcflurry-downloads fetch'"; ((errors++)); }
from mhcflurry import Class1AffinityPredictor
p = Class1AffinityPredictor.load()
p.predict(peptides=["GILGFVFTL"], allele="HLA-A0201")
PY
    [ ${errors} -eq 0 ] && log_info "MHCflurry models: OK"

    # --- SnpEff database ---
    if snpEff databases 2>/dev/null | grep -q '^GRCh38.p14'; then
        log_info "SnpEff database GRCh38.p14: available"
    else
        log_warn "GRCh38.p14 not listed by 'snpEff databases'"
        log_warn "  it may still be installed locally; SnpEff will download if not"
    fi

    # --- inputs ---
    [ -f "${PY_SCRIPT}" ] && log_info "worker script: OK" \
        || { log_error "not found: ${PY_SCRIPT}"; ((errors++)); }
    [ -s "${CATALOG}" ] \
        && log_info "variant catalog: OK ($(grep -vc '^#' ${CATALOG}) PASS variants)" \
        || { log_error "not found: ${CATALOG}. Run Step05b."; ((errors++)); }
    [ -s "${BEAD_MAP}" ] \
        && log_info "bead map: OK ($(( $(wc -l < ${BEAD_MAP}) - 1 )) pairs)" \
        || { log_error "not found: ${BEAD_MAP}. Run Step05c."; ((errors++)); }
    [ -f "${GENOME}" ] && log_info "genome FASTA: OK" \
        || { log_error "not found: ${GENOME}"; ((errors++)); }
    [ -f "${GENOME}.fai" ] && log_info "FASTA index: OK" \
        || { log_error "not found: ${GENOME}.fai"; ((errors++)); }

    if [ -f "${PROTEOME}" ]; then
        log_info "proteome: OK ($(grep -c '^>' ${PROTEOME}) isoforms, expect ~245,535)"
    elif [ -f "${PROTEOME}.gz" ]; then
        log_info "proteome (gzipped): OK"
    else
        log_error "proteome not found: ${PROTEOME}"
        log_error "  wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"
        ((errors++))
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
Step 06: Neoantigen prediction (spatial)

Usage: sbatch $(basename $0)
       $(basename $0) --validate-only
       $(basename $0) --help

Annotates the epithelial PASS variants with SnpEff, subtracts anything also
called outside the epithelium, builds 8-11mer mutant and wild-type peptides
from real Ensembl r115 protein context, scores them with MHCflurry against
the network paper's 10-allele panel, and maps every neoantigen back to the
beads that carry it.

Edit the CONFIGURATION block of the .py for the tumor cell type, thresholds,
or the HLA panel. Note that changing the panel breaks the Step07 overlap.
EOF
    exit 0
fi

section "STEP 06: NEOANTIGEN PREDICTION"
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
    log_error "Step06 failed with exit code ${RC}"
    exit 1
fi

section "VERIFYING OUTPUT"
for f in "${OUT}/epithelial_neoantigens_per_mutation.tsv" \
         "${OUT}/epithelial_all_peptide_results.tsv" \
         "${OUT}/proteome_mapping_diagnostics.tsv" \
         "${OUT}/step06_report.txt"; do
    [ -s "${f}" ] && log_info "OK: $(basename ${f}) ($(du -h ${f} | cut -f1))" \
        || { log_error "missing or empty: ${f}"; exit 1; }
done
[ -s "${OUT}/neoantigens_per_bead.tsv" ] \
    && log_info "OK: neoantigens_per_bead.tsv ($(( $(wc -l < ${OUT}/neoantigens_per_bead.tsv) - 1 )) rows)" \
    || log_warn "neoantigens_per_bead.tsv absent: no neoantigen had single-bead support"

section "STEP 06 COMPLETE"
log_success "Finished at $(date)"
log_info ""
log_info "Next: Step07 overlap against the network paper by (gene, hgvs_p) and"
log_info "by peptide, then the spatial and per-bead figures."

exit 0
