#!/bin/bash
#SBATCH -J SPATIAL_03_QC
#SBATCH -o /work/sdz852/WORKING/LOGS/Step03_QC.o.%j.log
#SBATCH -e /work/sdz852/WORKING/LOGS/Step03_QC.e.%j.log
#SBATCH --mail-user=jake.lehle@utsa.edu
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -p compute1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
# ============================================================================
# Step 03: QC filtering and puck merging (UTSA)
#
# Reads from:  data/outputs/02_anndata/
# Writes to:   data/outputs/03_qc/
#
# CPU-only (MAD outliers, sparse conversion, plotting, inner-join merge). No GPU.
#
# Memory: Step02 writes DENSE h5ads (~12-14 GB each on disk). All three are
# loaded at once, so the baseline is ~40 GB plus QC copies and the transient
# dense->sparse conversion. Estimated peak is ~60-80 GB; 128G gives headroom.
# TX Biomed's 500G was very conservative. If the log shows an OOM kill, switch
# "-p compute1" to "-p bigmem" and raise --mem; the real fix is upstream (see
# note below). The run is short (30-60 min), so an OOM retry is cheap.
#
# NOTE: the 500G ask exists only because Step02 stores X dense. If Step02 ever
# gets re-run to write sparse counts, this step's memory drops to a few GB. Not
# changing that now since Step02 already produced its outputs cleanly.
#
# Expected runtime: 30-60 minutes
# ============================================================================
module load anaconda3
conda activate spatial

SCRIPT_DIR="/work/sdz852/WORKING/slide-TCR-seq-working/scripts/HNSCC_slide-TCR-seq/scripts/UTSA"
conda run -n spatial python "${SCRIPT_DIR}/Step03_QC_Filter.py"

exit
