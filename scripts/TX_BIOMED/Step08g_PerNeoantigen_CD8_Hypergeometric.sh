#!/usr/bin/env bash
#SBATCH --job-name=cd8_peranchor
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/cd8_peranchor_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/cd8_peranchor_%j.err
# =============================================================================
# Step08g_PerNeoantigen_CD8_Hypergeometric.sh
# =============================================================================
# Per-anchor test on EVERY predicted neoantigen: is any TCR clone locally
# over-represented on CD8-dominated beads around this specific neoantigen?
#
#   sbatch Step08g_PerNeoantigen_CD8_Hypergeometric.sh
#   sbatch Step08g_PerNeoantigen_CD8_Hypergeometric.sh --diffusion-radius 30
#   sbatch Step08g_PerNeoantigen_CD8_Hypergeometric.sh --n-perm 500
#   sbatch Step08g_PerNeoantigen_CD8_Hypergeometric.sh --skip-spike
#   bash   Step08g_PerNeoantigen_CD8_Hypergeometric.sh --validate-only
#   bash   Step08g_PerNeoantigen_CD8_Hypergeometric.sh --help
#   sbatch Step08g_PerNeoantigen_CD8_Hypergeometric.sh --script /path/to/worker.py
#
# NO rarefaction floor, so no dependence on Step08c. It needs Step08b's MEASURED
# diffusion radius, resolved from 10_tcr_salvage/salvage_summary.tsv. If the
# schema cannot be resolved the worker STOPS and asks for --diffusion-radius
# rather than guessing. Same discipline as the floors, for the same reason.
#
# The spike-in now sweeps the clone's PUCK-WIDE focus count K as well as its
# local count k, so its runtime is 5x the first build's. --skip-spike drops it.
#
# RUNTIME is dominated by the rotation calibration (--n-perm rotations x anchors
# x two puck-wide arms x three pucks). --n-perm 50 for a fast look; the
# empirical p resolution is bounded by 1/(1+n_null), so a hit list that needs to
# survive BH wants the default or more.
#
# WORKER PATH RESOLUTION
#   Under sbatch, SLURM copies this script to /var/spool/slurmd/job<N>/ and runs
#   it from there, so ${BASH_SOURCE[0]} points at the spool copy and the Python
#   worker is not beside it. Resolution order:
#     1. --script PATH            explicit override
#     2. $SLURM_SUBMIT_DIR        where sbatch was invoked
#     3. dirname $BASH_SOURCE     correct for interactive `bash` runs
#     4. CANONICAL_DIR            the repo location, as a last resort
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
CANONICAL_DIR="$PROOT/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED"
PYNAME=Step08g_PerNeoantigen_CD8_Hypergeometric.py
ENV_NAME=slide-TCR-seq

PY_OVERRIDE=""
PASS_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --script) PY_OVERRIDE="$2"; shift 2 ;;
        -h|--help) sed -n '10,40p' "$0"; exit 0 ;;
        *) PASS_ARGS+=("$1"); shift ;;
    esac
done

echo "=============================================================="
echo "Step08g_PerNeoantigen_CD8_Hypergeometric"
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
python - <<'PYCHK' || exit 1
import sys, anndata, pandas, numpy, scipy, matplotlib
print("python", ".".join(str(x) for x in sys.version_info[:3]),
      "| anndata", anndata.__version__, "| pandas", pandas.__version__,
      "| numpy", numpy.__version__, "| scipy", scipy.__version__,
      "| matplotlib", matplotlib.__version__)
PYCHK

[[ -d "$PROOT" ]] || { echo "FAILED: $PROOT not found"; exit 1; }
mkdir -p "$PROOT/data/outputs/15_per_neoantigen_cd8/figures"

echo
python "$PY" "${PASS_ARGS[@]}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date)  (exit $RC)"
echo "  outputs:  $PROOT/data/outputs/15_per_neoantigen_cd8/"
ls -la "$PROOT/data/outputs/15_per_neoantigen_cd8/" 2>/dev/null | sed 's/^/    /'
echo "=============================================================="
exit $RC
