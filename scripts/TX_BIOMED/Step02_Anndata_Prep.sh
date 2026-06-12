#!/bin/bash

#SBATCH -J SPATIAL_02_ANNDATA
#SBATCH -o /master/jlehle/WORKING/LOGS/Step02_AnnData.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step02_AnnData.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 1-00:00:00
#SBATCH -p normal
#SBATCH --mem=800G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 80

# ============================================================================
# Step 02: Load Sophia's processed Slide-seq data into AnnData
#
# Reads from:  data/inputs/fastq/2022-01-28_Puck_*/
# Writes to:   data/outputs/02_anndata/
#
# Expected runtime: 10-30 minutes (loading DGE matrices + coordinate matching)
# ============================================================================

source ~/anaconda3/bin/activate
conda activate slide-TCR-seq

SCRIPT_DIR="/master/jlehle/WORKING/slide-TCR-seq-working/scripts"

python "${SCRIPT_DIR}/Step02_AnnData_Prep.py"
