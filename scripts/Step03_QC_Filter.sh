#!/bin/bash

#SBATCH -J SPATIAL_03_QC
#SBATCH -o /master/jlehle/WORKING/LOGS/Step03_QC.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step03_QC.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 1-00:00:00
#SBATCH -p normal
#SBATCH --mem=500G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16

# ============================================================================
# Step 03: QC filtering and puck merging
#
# Reads from:  data/outputs/02_anndata/
# Writes to:   data/outputs/03_qc/
#
# Needs ~500G because the dense h5ad files from Step02 are 12-14 GB each.
# After QC and sparse conversion, downstream files will be much smaller.
# Expected runtime: 30-60 minutes
# ============================================================================

source ~/anaconda3/bin/activate
conda activate slide-TCR-seq

SCRIPT_DIR="/master/jlehle/WORKING/slide-TCR-seq-working/scripts"

python "${SCRIPT_DIR}/Step03_QC_Filter.py"
