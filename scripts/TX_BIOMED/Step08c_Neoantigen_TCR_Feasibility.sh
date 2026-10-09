#!/usr/bin/env bash
#SBATCH --job-name=neo_tcr_feas
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/neo_tcr_feas_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/neo_tcr_feas_%j.err
# =============================================================================
# Step08c_Neoantigen_TCR_Feasibility.sh
# =============================================================================
# READ-ONLY. Decides whether the neoantigen-anchored TCR analysis is runnable,
# and builds the case/control sets it would use. Runs no hypothesis test.
#
#   sbatch Step08c_Neoantigen_TCR_Feasibility.sh
#   sbatch Step08c_Neoantigen_TCR_Feasibility.sh --chain TRB
#   sbatch Step08c_Neoantigen_TCR_Feasibility.sh --radius 100
#   bash   Step08c_Neoantigen_TCR_Feasibility.sh --validate-only
#
# Run it twice, once with no --chain and once with --chain TRB, so the cost of
# the TRB-primary definition is visible. Outputs are suffixed by chain.
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
CANONICAL_DIR="$PROOT/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED"
PYNAME=Step08c_Neoantigen_TCR_Feasibility.py
ENV_NAME=slide-TCR-seq

PY_OVERRIDE=""
PASS_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --script) PY_OVERRIDE="$2"; shift 2 ;;
        -h|--help) sed -n '10,26p' "$0"; exit 0 ;;
        *) PASS_ARGS+=("$1"); shift ;;
    esac
done

echo "=============================================================="
echo "Step08c_Neoantigen_TCR_Feasibility"
echo "  started: $(date)"
echo "  host:    $(hostname)"
echo "  job:     ${SLURM_JOB_ID:-interactive}"
echo "  submit:  ${SLURM_SUBMIT_DIR:-n/a}"
echo "  args:    ${PASS_ARGS[*]:-none}"
echo "=============================================================="

PY=""
CANDIDATES=()
[[ -n "$PY_OVERRIDE" ]] && CANDIDATES+=("$PY_OVERRIDE")
[[ -n "${SLURM_SUBMIT_DIR:-}" ]] && CANDIDATES+=("$SLURM_SUBMIT_DIR/$PYNAME")
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" || SELF_DIR=""
[[ -n "$SELF_DIR" ]] && CANDIDATES+=("$SELF_DIR/$PYNAME")
CANDIDATES+=("$CANONICAL_DIR/$PYNAME")
CANDIDATES+=("$CANONICAL_DIR/TROUBLE_SHOOTING/$PYNAME")

echo "resolving worker:"
for c in "${CANDIDATES[@]}"; do
    if [[ -f "$c" ]]; then echo "  FOUND   $c"; PY="$c"; break
    else echo "  no      $c"; fi
done
[[ -n "$PY" ]] || { echo; echo "FAILED: could not locate $PYNAME"; \
    echo "  sbatch $0 --script /full/path/to/$PYNAME"; exit 1; }

source ~/anaconda3/bin/activate "$ENV_NAME" || { echo "FAILED to activate $ENV_NAME"; exit 1; }
python -c "import anndata,pandas,numpy,scipy,matplotlib; \
print('anndata',anndata.__version__,'| scipy',scipy.__version__)" || exit 1

mkdir -p "$PROOT/data/outputs/11_neoantigen_tcr/figures"

echo
python "$PY" "${PASS_ARGS[@]}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date)  (exit $RC)"
ls -la "$PROOT/data/outputs/11_neoantigen_tcr/" 2>/dev/null | sed 's/^/    /'
echo "=============================================================="
exit $RC
