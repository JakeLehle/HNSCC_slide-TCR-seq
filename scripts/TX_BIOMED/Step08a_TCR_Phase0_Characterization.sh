#!/usr/bin/env bash
#SBATCH --job-name=tcr_phase0
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/tcr_phase0_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/tcr_phase0_%j.err
# =============================================================================
# Step08a_TCR_Phase0_Characterization.sh
# =============================================================================
# READ-ONLY. Characterizes the processed TCR tables before any analysis
# decision is made. Sets no thresholds and applies no filters.
#
# Gated on a Step05c preflight: every neoantigen bead ID must be a member of
# the annotated bead set. Those IDs are reconstructions from restore_cb() and
# that rerun was never confirmed, so a bad restore would silently dilute every
# downstream join rather than failing.
#
#   sbatch Step08a_TCR_Phase0_Characterization.sh
#   bash   Step08a_TCR_Phase0_Characterization.sh --validate-only
#   bash   Step08a_TCR_Phase0_Characterization.sh --help
#   sbatch Step08a_TCR_Phase0_Characterization.sh --script /path/to/worker.py
#
# Memory is set high because the annotated object is 99,341 x 19,822 and the
# three TCR tables together carry ~270,000 rows.
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
CANONICAL_DIR="$PROOT/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED"
PYNAME=Step08a_TCR_Phase0_Characterization.py
ENV_NAME=slide-TCR-seq

PY_OVERRIDE=""
PASS_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --script) PY_OVERRIDE="$2"; shift 2 ;;
        -h|--help) sed -n '10,30p' "$0"; exit 0 ;;
        *) PASS_ARGS+=("$1"); shift ;;
    esac
done

echo "=============================================================="
echo "Step08a_TCR_Phase0_Characterization"
echo "  started: $(date)"
echo "  host:    $(hostname)"
echo "  job:     ${SLURM_JOB_ID:-interactive}"
echo "  submit:  ${SLURM_SUBMIT_DIR:-n/a}"
echo "  args:    ${PASS_ARGS[*]:-none}"
echo "=============================================================="

# SLURM runs this from /var/spool/slurmd/job<N>/, so BASH_SOURCE is not enough
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
python - <<'PYCHK' || exit 1
import anndata, pandas, numpy, scipy, matplotlib
print("anndata", anndata.__version__, "| pandas", pandas.__version__,
      "| numpy", numpy.__version__, "| scipy", scipy.__version__,
      "| matplotlib", matplotlib.__version__)
try:
    import scirpy
    print("scirpy", scirpy.__version__, "(available for Phase 2)")
except ImportError:
    print("scirpy NOT installed. Needed for Phase 2, not for Phase 0.")
    print("  conda install -c conda-forge -c bioconda scirpy")
PYCHK

mkdir -p "$PROOT/data/outputs/09_tcr_phase0/figures"

echo
python "$PY" "${PASS_ARGS[@]}"
RC=$?
if (( RC != 0 )); then
    echo
    echo "worker exited $RC"
    echo "If the preflight failed, rerun Step05c -> Step06 -> Step07 first:"
    echo "  S=$PROOT/data/outputs/05_mutations/SComatic/SingleCell"
    echo "  rm -f \$S/checkpoints/phase2_filter.done \$S/checkpoints/phase3_feasibility.done"
    echo "  sbatch Step05c_SingleCellGenotype.sh"
    exit $RC
fi

echo
echo "=============================================================="
echo "  finished: $(date)"
ls -la "$PROOT/data/outputs/09_tcr_phase0/" | sed 's/^/    /'
echo "=============================================================="
