#!/bin/bash
#SBATCH -J SPATIAL_04_ANNOTATE
#SBATCH -o /work/sdz852/WORKING/LOGS/Step04_Annotate.o.%j.log
#SBATCH -e /work/sdz852/WORKING/LOGS/Step04_Annotate.e.%j.log
#SBATCH --mail-user=jake.lehle@utsa.edu
#SBATCH --mail-type=ALL
#SBATCH -t 3-00:00:00
#SBATCH -p gpu1a100
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
# ============================================================================
# Step 04: Cell type annotation, v5 (popV + cell2location + Pearson embedding)
#
# Reads:  data/outputs/03_qc/all_pucks_merged_QC.h5ad
#         data/inputs/reference/reference_for_cell2location.h5ad   (c2l ref)
# Writes: data/outputs/04_annotation/
#
# GPU step. cell2location's 30000-epoch spatial mapping is the whole reason
# this moved to UTSA; it runs on one a100 (gpu1a100, --gres=gpu:1). Wall time
# is capped at the GPU partition max of 3 days (was 7d on TX Biomed, far more
# than an a100 needs). --mem and -c are estimates for the a100 node; tune once
# the node specs are known. The script is device-aware (accelerator='gpu' if
# torch.cuda.is_available() else 'cpu'), so the only thing that matters is that
# the env's torch actually sees the GPU (see pre-flight step 1).
#
# -----------------------------------------------------------------------------
# PRE-FLIGHT (fresh start = no checkpoints, so the FULL heavy path runs:
# popV from scratch + reference regression + 30000-epoch map):
#
# 1) Confirm the env sees CUDA on a GPU node (this is what bit us on TX Biomed).
#      srun --partition=gpu1a100 --gres=gpu:1 -N 1 -n 1 -c 10 --time=00:30:00 --pty bash
#      module load anaconda3; conda activate spatial
#      python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
#    Want a real torch.version.cuda (e.g. 12.x) and True. If cuda is None or it
#    prints False, the env pulled a CPU-only torch wheel; fix before submitting.
#
# 2) popV pulls the Tabula Sapiens model (~2 GB) from HuggingFace on first run.
#    If a100 nodes have no outbound internet, pre-stage from a login node:
#      module load anaconda3; conda activate spatial
#      python -c "import popv; popv.hub.HubModel.pull_from_huggingface_hub('popV/tabula_sapiens_All_Cells', cache_dir='/work/sdz852/WORKING/slide-TCR-seq-working/data/outputs/04_annotation/popv_cache')"
#    Then uncomment the HF_HUB_OFFLINE line below so the GPU job never reaches
#    for the network.
#
# 3) Stage the c2l reference at the config path:
#      data/inputs/reference/reference_for_cell2location.h5ad
#
# 4) Confirm the spatial env has: popv, cell2location 0.1.5, scvi-tools
#    1.4.0.post1, bbknn, adjustText, squidpy. Missing bbknn is silently caught
#    and clustering runs UNCORRECTED (you lose batch correction), so check it.
#
# Expected runtime on a100: a few hours for the map + ~1-2h popV, well under 3d.
# ============================================================================
module load anaconda3
conda activate spatial

# export HF_HUB_OFFLINE=1   # uncomment AFTER pre-staging the popV model (step 2)

# NOTE: --no-capture-output is REQUIRED. Plain `conda run` buffers the child's
# stdout/stderr and only flushes on a clean exit, so a cancelled job (e.g. an
# idle-GPU reaper kill) discards the entire log and you get an empty file. With
# --no-capture-output the log streams live, which is what lets us see popV's
# per-member lightning device banners ("GPU available ... used: True/False").
SCRIPT_DIR="/work/sdz852/WORKING/slide-TCR-seq-working/scripts/HNSCC_slide-TCR-seq/scripts/UTSA"
conda run --no-capture-output -n spatial python "${SCRIPT_DIR}/Step04b_Cell_Type_Annotation.py"
