#!/bin/bash
#SBATCH -J SPATIAL_04a_POPV
#SBATCH -o /work/sdz852/WORKING/LOGS/Step04a_popV.o.%j.log
#SBATCH -e /work/sdz852/WORKING/LOGS/Step04a_popV.e.%j.log
#SBATCH --mail-user=jake.lehle@utsa.edu
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -p compute2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 80
# ============================================================================
# Step 04a: popV ONLY, on CPU. First half of the decoupled Step04.
#
# Reads:  data/outputs/03_qc/all_pucks_merged_QC.h5ad
# Writes: data/outputs/04_annotation/all_pucks_popv_annotated.h5ad
#
# Uses the ORIGINAL 'spatial' env (untouched, known-good). No GPU: popV holds
# no card, so the idle-GPU reaper cannot kill it, and the slow CPU ensemble
# just runs to completion once and checkpoints. compute2 gives a long wall.
#
# popV downloads the Tabula Sapiens model (~2 GB) on first run. If compute2
# nodes have no outbound internet, pre-stage it from a login node first:
#   module load anaconda3; conda activate spatial
#   python -c "import popv; popv.hub.HubModel.pull_from_huggingface_hub('popV/tabula_sapiens_All_Cells', cache_dir='/work/sdz852/WORKING/slide-TCR-seq-working/data/outputs/04_annotation/popv_cache')"
#
# After this finishes, submit the GPU Step04 launcher: it finds the checkpoint,
# skips popV, and runs cell2location on the a100.
# ============================================================================
module load anaconda3
conda activate spatial

SCRIPT_DIR="/work/sdz852/WORKING/slide-TCR-seq-working/scripts/HNSCC_slide-TCR-seq/scripts/UTSA"
conda run --no-capture-output -n spatial python "${SCRIPT_DIR}/Step04a_popV_only.py"
