#!/bin/bash

#SBATCH -J SPATIAL_04_ANNOTATE
#SBATCH -o /master/jlehle/WORKING/LOGS/Step04_Annotate.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step04_Annotate.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00
#SBATCH -p normal
#SBATCH --mem=900G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 180

# ============================================================================
# Step 04: Cell type annotation (popV + Leiden clustering)
#
# Reads from:  data/outputs/03_qc/all_pucks_merged_QC.h5ad
# Writes to:   data/outputs/04_annotation/
#
# NOTE: Using NETWORK env (confirmed to have popV, scanpy, squidpy).
# If using a different env, update the conda activate line below.
#
# popV downloads the Tabula Sapiens model on first run (~2 GB).
# If compute nodes lack internet, pre-download on a login node first:
#   conda activate NETWORK
#   python -c "
#   import popv
#   popv.hub.HubModel.pull_from_huggingface_hub(
#       'popV/tabula_sapiens_All_Cells',
#       cache_dir='/master/jlehle/WORKING/slide-TCR-seq-working/data/outputs/04_annotation/popv_cache'
#   )"
#
# Expected runtime: 1-3 hours (popV is the bottleneck)
# ============================================================================

source ~/anaconda3/bin/activate
conda activate sc_pre

SCRIPT_DIR="/master/jlehle/WORKING/slide-TCR-seq-working/scripts"

python "${SCRIPT_DIR}/Step04_Cell_Type_Annotation.py"
