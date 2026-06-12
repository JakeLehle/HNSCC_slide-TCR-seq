#!/bin/bash
#SBATCH -J SPATIAL_02_ANNDATA
#SBATCH -o /work/sdz852/WORKING/LOGS/Step02_AnnData.o.%j.log
#SBATCH -e /work/sdz852/WORKING/LOGS/Step02_AnnData.e.%j.log
#SBATCH --mail-user=jake.lehle@utsa.edu
#SBATCH --mail-type=ALL
#SBATCH -t 04:00:00
#SBATCH -p compute1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
# ============================================================================
# Step 02: Load Sophia's processed Slide-seq data into AnnData (UTSA)
#
# Reads from:  data/inputs/fastq/2022-01-28_Puck_*/
# Writes to:   data/outputs/02_anndata/
#
# CPU-only step: pandas reads of the matched digital-expression text matrices
# plus barcode-to-coordinate matching. No GPU, so this goes to compute1.
#
# Memory note: the .matched.digital_expression.txt.gz files are read with
# pandas as DENSE gene x barcode matrices, and all three per-puck AnnData
# objects are held in memory at once before saving. Size for the dense
# footprint, not the sparse one. 128G is generous; tune down after the first
# run once the log shows actual peak usage. (TX Biomed's 800G was overkill.)
#
# Expected runtime: 10-30 minutes
# ============================================================================
module load anaconda3
conda activate spatial

SCRIPT_DIR="/work/sdz852/WORKING/slide-TCR-seq-working/scripts/HNSCC_slide-TCR-seq/scripts/UTSA"
conda run -n spatial python "${SCRIPT_DIR}/Step02_AnnData_Prep.py"

exit
