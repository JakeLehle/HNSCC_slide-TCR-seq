#!/usr/bin/env bash
#SBATCH --job-name=cd8_tests
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/cd8_tests_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/cd8_tests_%j.err
# =============================================================================
# Step08f_CD8_Neoantigen_Tests.sh
# =============================================================================
# The CD8-conditioned neoantigen tests. Designs 2 and 3 from Step08e; design 1
# is documented as not runnable rather than tested.
#
#   sbatch Step08f_CD8_Neoantigen_Tests.sh
#   sbatch Step08f_CD8_Neoantigen_Tests.sh --chain TRB --radius 100
#   sbatch Step08f_CD8_Neoantigen_Tests.sh --skip-design3
#   sbatch Step08f_CD8_Neoantigen_Tests.sh --skip-partC
#   bash   Step08f_CD8_Neoantigen_Tests.sh --validate-only
#   bash   Step08f_CD8_Neoantigen_Tests.sh --help
#   sbatch Step08f_CD8_Neoantigen_Tests.sh --script /path/to/worker.py
#
# ORDER OF OPERATIONS
#   Step08c --chain <X>  ->  Step08e  ->  Step08f --chain <X>
#   Floors are read from feasibility_summary_<X>.tsv and the worker hard-stops
#   if it is absent. Never hardcoded; that was the v2 bug. Step08e's
#   pair_attrition.tsv feeds PART F, which degrades gracefully if missing.
#
# PART D dominates the runtime: 400 rarefaction replicates per neighbourhood
# across ~103 pairs x 2 members. --skip-design3 cuts it to minutes.
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
PYNAME=Step08f_CD8_Neoantigen_Tests.py
ENV_NAME=slide-TCR-seq

PY_OVERRIDE=""
PASS_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --script) PY_OVERRIDE="$2"; shift 2 ;;
        -h|--help) sed -n '10,45p' "$0"; exit 0 ;;
        *) PASS_ARGS+=("$1"); shift ;;
    esac
done

echo "=============================================================="
echo "Step08f_CD8_Neoantigen_Tests"
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
mkdir -p "$PROOT/data/outputs/14_cd8_tests/figures"

echo
python "$PY" "${PASS_ARGS[@]}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date)  (exit $RC)"
echo "  outputs:  $PROOT/data/outputs/14_cd8_tests/"
ls -la "$PROOT/data/outputs/14_cd8_tests/" 2>/dev/null | sed 's/^/    /'
echo "=============================================================="
exit $RC
