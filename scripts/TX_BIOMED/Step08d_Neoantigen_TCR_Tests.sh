#!/usr/bin/env bash
#SBATCH --job-name=neo_tcr_tests
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/neo_tcr_tests_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/neo_tcr_tests_%j.err
# =============================================================================
# Step08d_Neoantigen_TCR_Tests.sh
# =============================================================================
# The four neoantigen-anchored TCR tests. Every setting was fixed in Step08c
# BEFORE any outcome was computed.
#
#   sbatch Step08d_Neoantigen_TCR_Tests.sh                 # TRB, 75 and 100 um
#   sbatch Step08d_Neoantigen_TCR_Tests.sh --chain all     # sensitivity arm
#   sbatch Step08d_Neoantigen_TCR_Tests.sh --skip-C        # C is the slow test
#   bash   Step08d_Neoantigen_TCR_Tests.sh --validate-only
#
# REQUIRES Step08c outputs: matched_pairs.tsv, testable_cases.tsv,
# background_pool.tsv. Run 08c first if any are missing.
#
# Test C does N_ROT=1000 rotations x 3 lambdas x 3 pucks x 2 radii, each a
# least-squares fit over ~17,000 epithelial beads. Expect it to dominate the
# runtime; --skip-C drops it.
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
CANONICAL_DIR="$PROOT/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED"
PYNAME=Step08d_Neoantigen_TCR_Tests.py
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
echo "Step08d_Neoantigen_TCR_Tests"
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

FEAS="$PROOT/data/outputs/11_neoantigen_tcr"
for f in matched_pairs.tsv testable_cases.tsv background_pool.tsv; do
    if [[ -f "$FEAS/$f" ]]; then echo "  OK      $FEAS/$f"
    else echo "  MISSING $FEAS/$f  -> run Step08c first"; exit 1; fi
done

mkdir -p "$PROOT/data/outputs/12_neoantigen_tcr_tests/figures"

echo
python "$PY" "${PASS_ARGS[@]}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date)  (exit $RC)"
ls -la "$PROOT/data/outputs/12_neoantigen_tcr_tests/" 2>/dev/null | sed 's/^/    /'
echo "=============================================================="
exit $RC
