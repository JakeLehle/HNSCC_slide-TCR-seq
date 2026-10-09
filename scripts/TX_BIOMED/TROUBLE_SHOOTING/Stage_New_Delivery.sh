#!/usr/bin/env bash
#SBATCH --job-name=stage_delivery
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/stage_delivery_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/stage_delivery_%j.err
# =============================================================================
# Stage_New_Delivery.sh
# =============================================================================
# Stages the 2026-09 Sophia delivery from BACKUP into data/inputs, so that
# data/inputs becomes the single complete source that the archive is built from.
#
# DESIGN RULES
#   1. BACKUP is never modified and never emptied. It is the only evidence of
#      what was delivered. Small files are COPIED, large files SYMLINKED.
#   2. Nothing is ever overwritten. A size mismatch is a stop, not a merge.
#   3. Phase 2 refuses to run until Diagnostic_Backup_vs_Inputs_Verify.sh has
#      produced its lists.
#
# USAGE
#   bash   Stage_New_Delivery.sh                         # dry run, prints plan
#   sbatch Stage_New_Delivery.sh --go                    # phase 1 (TCR, 36 GB cp)
#   sbatch Stage_New_Delivery.sh --go --phase2           # + symlink new paths
#   sbatch Stage_New_Delivery.sh --go --phase2 --phase3  # + samtools index
#
# Phase 1 copies ~36 GB of ONT FASTQ, which is why this runs on the scheduler.
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
BK="$PROOT/BACKUP"
IN="$PROOT/data/inputs/fastq"
TCR="$PROOT/data/inputs/tcr"
INV="$PROOT/data/outputs/00_inventory"
ENV_NAME=slide-TCR-seq

PUCKS=(29 37 40)
GO=0; PHASE2=0; PHASE3=0

for a in "$@"; do
    case "$a" in
        --go)     GO=1 ;;
        --phase2) PHASE2=1 ;;
        --phase3) PHASE3=1 ;;
        -h|--help) sed -n '9,32p' "$0"; exit 0 ;;
        *) echo "unknown arg: $a"; exit 1 ;;
    esac
done

run() {
    if (( GO )); then
        echo "  RUN: $*"
        "$@" || { echo "  FAILED: $*"; exit 1; }
    else
        echo "  DRY: $*"
    fi
}
banner() { echo; echo "===================================================================="; echo "$1"; echo "===================================================================="; }

echo "=============================================================="
echo "Stage_New_Delivery"
echo "  started: $(date)"
echo "  host:    $(hostname)"
echo "  job:     ${SLURM_JOB_ID:-interactive}"
echo "  args:    $*"
echo "=============================================================="
(( GO )) || echo "*** DRY RUN. Nothing will be changed. Add --go to act. ***"

# -----------------------------------------------------------------------------
banner "PHASE 0  protect the as-delivered archive"
# -----------------------------------------------------------------------------
# Stops an accidental mv out of BACKUP later in the project.
run chmod -R a-w "$BK"

# -----------------------------------------------------------------------------
banner "PHASE 1  TCR files  (no counterpart in data/inputs)"
# -----------------------------------------------------------------------------
run mkdir -p "$TCR/processed" "$TCR/ont"

declare -A CSV_MAP=(
  ["B59_29_hTCR_tcr.csv"]="Puck_211214_29"
  ["B59_37_hTCR_tcr.csv"]="Puck_211214_37"
  ["B59_40_hTCR_tcr.csv"]="Puck_211214_40"
)
declare -A ONT_MAP=(
  ["TCR_20220224_Puck_211214_29.gz"]="Puck_211214_29"
  ["TCR_20220228_Puck_211214_37.gz"]="Puck_211214_37"
  ["TCR_20220127_Puck_211214_40.gz"]="Puck_211214_40"
)

for f in "${!CSV_MAP[@]}"; do
    if [[ -e "$TCR/processed/$f" ]]; then
        echo "  SKIP (exists): $TCR/processed/$f"
    else
        run cp -p "$BK/$f" "$TCR/processed/$f"
    fi
done

# ~36 GB of ONT FASTQ. rsync so an interrupted job resumes instead of restarting.
for f in "${!ONT_MAP[@]}"; do
    if [[ -e "$TCR/ont/$f" ]] && [[ $(stat -c%s "$TCR/ont/$f") -eq $(stat -c%s "$BK/$f") ]]; then
        echo "  SKIP (complete): $TCR/ont/$f"
    else
        run rsync -a --partial --inplace --info=progress2 "$BK/$f" "$TCR/ont/$f"
    fi
done

# provenance map so no script has to parse a puck out of a filename
MAPFILE="$TCR/README_puck_mapping.tsv"
if (( GO )); then
    {
        printf "file\tpuck_id\tsubdir\trole\tplatform\tdelivered\n"
        for f in "${!CSV_MAP[@]}"; do
            printf "%s\t%s\tprocessed\tMiXCR_clonotype_table\tIllumina\t2026-09\n" "$f" "${CSV_MAP[$f]}"
        done
        for f in "${!ONT_MAP[@]}"; do
            printf "%s\t%s\tont\trhTCRseq_Fraction2_reads\tOxfordNanopore\t2026-09\n" "$f" "${ONT_MAP[$f]}"
        done
    } | sort -k2,2 -k3,3 > "$MAPFILE"
    echo "  wrote $MAPFILE"
else
    echo "  DRY: write $MAPFILE"
fi

# -----------------------------------------------------------------------------
banner "PHASE 2  per-puck paths that are genuinely new  (symlink, no copy)"
# -----------------------------------------------------------------------------
if (( ! PHASE2 )); then
    echo "  skipped (pass --phase2 to enable)"
else
    MISS=0
    for P in "${PUCKS[@]}"; do
        for f in only_in_backup_puck${P}.txt size_mismatch_puck${P}.txt; do
            [[ -f "$INV/$f" ]] || { echo "  MISSING $INV/$f"; MISS=1; }
        done
    done
    if (( MISS )); then
        echo; echo "  STOP. Run Diagnostic_Backup_vs_Inputs_Verify.sh first."; exit 1
    fi

    # a size difference is a decision, not an automatic action
    NMM=0
    for P in "${PUCKS[@]}"; do
        n=$(wc -l < "$INV/size_mismatch_puck${P}.txt"); NMM=$((NMM + n))
        (( n > 0 )) && echo "  puck $P has $n size mismatches"
    done
    if (( NMM > 0 )); then
        echo; echo "  STOP. $NMM files differ in size between BACKUP and data/inputs."
        echo "  Decide per file whether inputs was truncated or Sophia reprocessed."
        echo "  See $INV/size_mismatch_puck*.txt"; exit 1
    fi

    for P in "${PUCKS[@]}"; do
        D="2022-01-28_Puck_211214_${P}"; N=0
        while IFS= read -r rel; do
            [[ -n "$rel" ]] || continue
            DST="$IN/$D/$rel"
            [[ -e "$DST" ]] && continue
            run mkdir -p "$(dirname "$DST")"
            run ln -s "$BK/$D/$rel" "$DST"
            N=$((N+1))
        done < "$INV/only_in_backup_puck${P}.txt"
        echo "  puck $P: $N new paths staged as symlinks"
    done
fi

# -----------------------------------------------------------------------------
banner "PHASE 3  generate what no delivery supplies"
# -----------------------------------------------------------------------------
if (( ! PHASE3 )); then
    echo "  skipped (pass --phase3 to enable)"
else
    source ~/anaconda3/bin/activate "$ENV_NAME" || { echo "FAILED to activate $ENV_NAME"; exit 1; }
    samtools --version | head -1

    # only puck 40 shipped a matched.bam index
    for P in 29 37; do
        B="$IN/2022-01-28_Puck_211214_${P}/Puck_211214_${P}.matched.bam"
        if [[ -f "$B.bai" ]]; then
            echo "  SKIP (indexed): $B.bai"
        elif [[ -f "$B" ]]; then
            run samtools index -@ "${SLURM_CPUS_PER_TASK:-8}" "$B"
        else
            echo "  MISSING: $B"
        fi
    done

    if (( GO )); then
        find "$PROOT/data/inputs" -type f -printf "%P\t%s\n" 2>/dev/null \
            | sort > "$PROOT/data/inputs/input_inventory.txt"
        echo "  refreshed input_inventory.txt ($(wc -l < "$PROOT/data/inputs/input_inventory.txt") files)"
    else
        echo "  DRY: refresh input_inventory.txt"
    fi
fi

banner "DONE"
echo "  finished: $(date)"
if (( GO )); then
    echo "  data/inputs is now the single complete source."
    du -sh "$PROOT/data/inputs"/* 2>/dev/null | sed 's/^/    /'
    echo
    echo "  NEXT"
    echo "    sbatch Diagnostic_Full_BAM_Inventory.sh --counts"
    echo "    sbatch Build_Archive_Tarballs.sh --go"
fi
