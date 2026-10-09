#!/usr/bin/env bash
#SBATCH --job-name=bam_inventory
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/bam_inventory_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/bam_inventory_%j.err
# =============================================================================
# Diagnostic_Full_BAM_Inventory.sh
# =============================================================================
# READ-ONLY wrapper for Diagnostic_Full_BAM_Inventory.py.
#
#   sbatch Diagnostic_Full_BAM_Inventory.sh                 # sampled, ~10 min
#   sbatch Diagnostic_Full_BAM_Inventory.sh --counts        # + exact counts
#   bash   Diagnostic_Full_BAM_Inventory.sh --validate-only # paths only
#   bash   Diagnostic_Full_BAM_Inventory.sh --help
#   sbatch Diagnostic_Full_BAM_Inventory.sh --script /path/to/worker.py
#
# The --counts pass proves the uBAMs hold the complete library, by checking
# uBAM record count against 2x the cellular_tagging total. Worth the wall time
# once, before the archive is built.
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
CANONICAL_DIR="$PROOT/scripts/HNSCC_slide-TCR-seq/scripts/TX_BIOMED/TROUBLE_SHOOTING"
PYNAME=Diagnostic_Full_BAM_Inventory.py
ENV_NAME=slide-TCR-seq

# ---- parse our own args, pass the rest through -------------------------------
PY_OVERRIDE=""
PASS_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --script) PY_OVERRIDE="$2"; shift 2 ;;
        -h|--help) sed -n '10,34p' "$0"; exit 0 ;;
        *) PASS_ARGS+=("$1"); shift ;;
    esac
done

echo "=============================================================="
echo "Diagnostic_Full_BAM_Inventory"
echo "  started: $(date)"
echo "  host:    $(hostname)"
echo "  job:     ${SLURM_JOB_ID:-interactive}"
echo "  submit:  ${SLURM_SUBMIT_DIR:-n/a}"
echo "  args:    ${PASS_ARGS[*]:-none}"
echo "=============================================================="

# ---- resolve the worker ------------------------------------------------------
PY=""
CANDIDATES=()
[[ -n "$PY_OVERRIDE" ]] && CANDIDATES+=("$PY_OVERRIDE")
[[ -n "${SLURM_SUBMIT_DIR:-}" ]] && CANDIDATES+=("$SLURM_SUBMIT_DIR/$PYNAME")
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" || SELF_DIR=""
[[ -n "$SELF_DIR" ]] && CANDIDATES+=("$SELF_DIR/$PYNAME")
CANDIDATES+=("$CANONICAL_DIR/$PYNAME")

echo "resolving worker:"
for c in "${CANDIDATES[@]}"; do
    if [[ -f "$c" ]]; then
        echo "  FOUND   $c"
        PY="$c"
        break
    else
        echo "  no      $c"
    fi
done

if [[ -z "$PY" ]]; then
    echo
    echo "FAILED: could not locate $PYNAME"
    echo "  Pass it explicitly:  sbatch $0 --script /full/path/to/$PYNAME"
    exit 1
fi

# ---- environment -------------------------------------------------------------
source ~/anaconda3/bin/activate "$ENV_NAME" || { echo "FAILED to activate $ENV_NAME"; exit 1; }
python -c "import pysam; print('pysam', pysam.__version__)" || exit 1
samtools --version | head -1 || exit 1

[[ -d "$PROOT" ]] || { echo "FAILED: $PROOT not found"; exit 1; }
mkdir -p "$PROOT/data/outputs/00_inventory"

# ---- run ---------------------------------------------------------------------
echo
python "$PY" "${PASS_ARGS[@]}" || { echo "FAILED: python worker exited nonzero"; exit 1; }

echo
echo "=============================================================="
echo "  finished: $(date)"
echo "  outputs:  $PROOT/data/outputs/00_inventory/"
ls -la "$PROOT/data/outputs/00_inventory/" | sed 's/^/    /'
echo "=============================================================="
