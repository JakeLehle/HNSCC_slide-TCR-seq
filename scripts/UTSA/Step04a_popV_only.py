#!/usr/bin/env python3
"""
Step04a_popV_only.py

Run ONLY the popV stage, on a CPU node, write the popV checkpoint, and exit.
This is the first half of the decoupled Step04: popV holds no GPU (so a GPU-idle
reaper cannot touch it), and the heavy cell2location half then runs as a separate
GPU job that loads this checkpoint and skips popV entirely.

It imports run_popv_annotation/safe_write straight from Step04b_Cell_Type_Annotation.py,
so the popV logic (popv_* whitelist, var_name restoration, query filtering) is a
single source of truth and the checkpoint is identical to the monolithic run.

Reads:   data/outputs/03_qc/all_pucks_merged_QC.h5ad
Writes:  data/outputs/04_annotation/all_pucks_popv_annotated.h5ad
Env:     spatial   (the original, untouched env -- NOT spatial_cuml)
Run on:  a CPU partition (compute1/compute2). No GPU needed.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import scanpy as sc
import popv

from spatial_config import (
    DIR_03_QC, DIR_04_ANNOTATION, ensure_dir, banner, log,
)
# Reuse the exact popV implementation from Step04b (the cell2location-onward
# half; renamed from Step04_Cell_Type_Annotation.py to Step04b_...py).
from Step04b_Cell_Type_Annotation import run_popv_annotation, safe_write


def main():
    banner("STEP 04a: popV ONLY (CPU node) -> checkpoint")

    # IMPORTANT: keep n_jobs = 1. popV reuses settings.n_jobs as the lightning
    # `devices` count for its scVI/scANVI members when cuml is off
    # (algorithms/_scvi.py: devices=[settings.device] if settings.cuml else
    # settings.n_jobs). So n_jobs > 1 launches multi-device DDP on CPU and
    # crashes scVI in the covariate encoder with
    # "expand(... {[1,1,512,281]}, size=[1,1,512])". n_jobs=1 -> single-device,
    # which is correct. The slow popV steps are the bbknn/harmony integrations,
    # which n_jobs does not accelerate anyway, so there is nothing to gain by
    # raising it. accelerator stays 'auto' (resolves to CPU here, no GPU).
    popv.settings.n_jobs = 1
    log(f"  popv.settings.n_jobs = {popv.settings.n_jobs}, "
        f"accelerator = {popv.settings.accelerator!r}")

    ann_dir = ensure_dir(DIR_04_ANNOTATION)
    cache_dir = os.path.join(ann_dir, "popv_cache")
    popv_path = os.path.join(ann_dir, "all_pucks_popv_annotated.h5ad")

    if os.path.exists(popv_path):
        log(f"  checkpoint already exists: {popv_path}")
        log(f"  nothing to do; the GPU job will load this. exiting.")
        return

    merged = os.path.join(DIR_03_QC, "all_pucks_merged_QC.h5ad")
    if not os.path.exists(merged):
        log(f"  FATAL: {merged} not found. Run Step03 first.")
        sys.exit(1)

    log(f"  loading {merged}")
    adata = sc.read_h5ad(merged)
    log(f"  loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    log(f"  pucks: {adata.obs['puck_id'].value_counts().to_dict()}")

    adata = run_popv_annotation(adata, cache_dir)
    safe_write(adata, popv_path)

    banner("STEP 04a COMPLETE")
    log(f"  checkpoint written: {popv_path}")
    log(f"  next: submit the GPU Step04 job; it will load this checkpoint,")
    log(f"        skip popV, and run cell2location on the a100.")


if __name__ == "__main__":
    main()
