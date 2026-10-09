#!/usr/bin/env bash
#SBATCH --job-name=tcr_salvage
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/tcr_salvage_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/tcr_salvage_%j.err
# =============================================================================
# Step08b_TCR_Salvage_Assessment.sh
# =============================================================================
# READ-ONLY. Decides whether the processed TCR tables can support spatial
# inference, and at what cost.
#
#   sbatch Step08b_TCR_Salvage_Assessment.sh
#   sbatch Step08b_TCR_Salvage_Assessment.sh --skip-ont     # BEAT 8 is the slow one
#   sbatch Step08b_TCR_Salvage_Assessment.sh --radius 50    # locality radius
#   bash   Step08b_TCR_Salvage_Assessment.sh --validate-only
#
# BEAT 8 streams ~500k reads from each of three gzipped ONT files (9-14 GB
# each). Expect roughly 10-20 minutes for that beat alone; --skip-ont drops it.
#
# statsmodels is REQUIRED here, unlike Step08a where the GLM was optional:
# BEAT 2 and BEAT 5 are both model based.
#   conda install -c conda-forge statsmodels
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
CANONICAL_DIR="$PROOT/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED"
PYNAME=Step08b_TCR_Salvage_Assessment.py
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
echo "Step08b_TCR_Salvage_Assessment"
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
import sys
import anndata, pandas, numpy, scipy, matplotlib
print("anndata", anndata.__version__, "| pandas", pandas.__version__,
      "| numpy", numpy.__version__, "| scipy", scipy.__version__)
try:
    import statsmodels
    print("statsmodels", statsmodels.__version__)
except ImportError:
    print("FATAL: statsmodels is required for BEAT 2 and BEAT 5.")
    print("  conda install -c conda-forge statsmodels")
    sys.exit(1)
PYCHK

mkdir -p "$PROOT/data/outputs/10_tcr_salvage/figures"

echo
python "$PY" "${PASS_ARGS[@]}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date)  (exit $RC)"
ls -la "$PROOT/data/outputs/10_tcr_salvage/" 2>/dev/null | sed 's/^/    /'
echo "=============================================================="
exit $RC
