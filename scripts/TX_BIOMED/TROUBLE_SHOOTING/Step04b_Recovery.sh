#!/bin/bash

#SBATCH -J SPATIAL_04b_RECOVERY
#SBATCH -o /master/jlehle/WORKING/LOGS/Step04b_Recovery.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step04b_Recovery.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 1-00:00:00
#SBATCH -p normal
#SBATCH --mem=900G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 80

# ============================================================================
# Step 04b: Recovery from popV crash
#
# popV annotation completed successfully (~5.5 hrs) but the h5ad write
# crashed due to nullable string dtype incompatibility. This script
# recovers using the saved predictions.csv and runs all post-processing.
#
# Expected runtime: ~30-60 minutes (PCA + UMAP + Leiden + plots, no popV)
# ============================================================================

source ~/anaconda3/bin/activate
conda activate sc_pre

SCRIPT_DIR="/master/jlehle/WORKING/slide-TCR-seq-working/scripts"

conda run -n sc_pre python "${SCRIPT_DIR}/Step04b_Recovery.py"
