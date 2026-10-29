#!/usr/bin/env bash
#SBATCH --job-name=cd8_phase0
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=08:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/cd8_phase0_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/cd8_phase0_%j.err
# =============================================================================
# Step08e_CD8_Phase0_Diagnostic.sh
# =============================================================================
# Wrapper for Step08e_CD8_Phase0_Diagnostic.py. Runs NO hypothesis test. Sizes
# whether a CD8-restricted version of the Step08d tests can run, and on what
# terms. Its one write outside its own output directory is a backup of the
# Step08c design artifacts, which are not chain-suffixed and are overwritten by
# every 08c run.
#
#   sbatch Step08e_CD8_Phase0_Diagnostic.sh
#   sbatch Step08e_CD8_Phase0_Diagnostic.sh --chain TRB
#   sbatch Step08e_CD8_Phase0_Diagnostic.sh --radii 30 50 75 100 150 200
#   sbatch Step08e_CD8_Phase0_Diagnostic.sh --skip-markers
#   bash   Step08e_CD8_Phase0_Diagnostic.sh --validate-only
#   bash   Step08e_CD8_Phase0_Diagnostic.sh --help
#   sbatch Step08e_CD8_Phase0_Diagnostic.sh --script /path/to/worker.py
#
# Memory is set at 160G because BEAT 3 slices gene columns out of the 99,341 x
# 19,822 annotated object. Pass --skip-markers and 64G is enough.
#
# Step08c --chain <X> MUST have run first: the rarefaction floors are read from
# feasibility_summary_<X>.tsv and the script hard-stops if it is absent. Floors
# are never hardcoded; that was the v2 bug.
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
PYNAME=Step08e_CD8_Phase0_Diagnostic.py
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
echo "Step08e_CD8_Phase0_Diagnostic"
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
import anndata, pandas, numpy, scipy, matplotlib
print("anndata", anndata.__version__, "| pandas", pandas.__version__,
      "| numpy", numpy.__version__, "| scipy", scipy.__version__,
      "| matplotlib", matplotlib.__version__)
PYCHK

[[ -d "$PROOT" ]] || { echo "FAILED: $PROOT not found"; exit 1; }
mkdir -p "$PROOT/data/outputs/13_cd8_phase0/figures"

echo
python "$PY" "${PASS_ARGS[@]}"
RC=$?

echo
echo "=============================================================="
echo "  finished: $(date)  (exit $RC)"
echo "  outputs:  $PROOT/data/outputs/13_cd8_phase0/"
ls -la "$PROOT/data/outputs/13_cd8_phase0/" 2>/dev/null | sed 's/^/    /'
echo "  backup:   $PROOT/data/outputs/11_neoantigen_tcr/BACKUP_pre_cd8_*/"
ls -d "$PROOT/data/outputs/11_neoantigen_tcr/BACKUP_pre_cd8_"* 2>/dev/null | sed 's/^/    /'
echo "=============================================================="
exit $RC
