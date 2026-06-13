#!/usr/bin/env python3
"""
Step04_popv_gpu_probe.py

Find out, on real data, WHETHER and WHEN popV uses the GPU, and which knobs the
installed popV exposes to control it. Run INTERACTIVELY on an a100:

    srun --partition=gpu1a100 --gres=gpu:1 -N 1 -n 1 -c 10 --time=01:00:00 --pty bash
    module load anaconda3
    conda activate spatial
    python Step04_popv_gpu_probe.py            # default 3000-bead subset

Optional:
    python Step04_popv_gpu_probe.py --n 5000 --h5ad /path/to/merged_QC.h5ad

What it does:
  [A] Introspect the installed popV API:
        - popv.__version__
        - attributes on popv.settings (look for accelerator / device / n_jobs)
        - signature of HubModel.annotate_data (look for accelerator=, methods=,
          devices=, etc.)  -> tells us if there is an explicit GPU knob and/or a
          way to drop the slow CPU-only members (OnClass, scanorama, ...).
  [B] Run popV on a SMALL real subset with a live GPU-utilization timeline
      printed every few seconds, so you can SEE the profile: which stretches
      pin the card (scVI/scANVI members) and which sit at ~0% (kNN/harmony/
      scanorama/RF/OnClass/celltypist). Long ~0% stretches are what an idle
      reaper trips on.

This does NOT change the pipeline. It is read-only except for the popV model
download into the cache dir (shared with the real run, so it is reused later).
"""

import os
import sys
import time
import argparse
import inspect
import threading
import subprocess


def query_gpu():
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used",
             "--format=csv,noheader,nounits"],
            stderr=subprocess.DEVNULL, timeout=10,
        ).decode().strip().splitlines()
        util, mem = out[0].split(",")
        return float(util), float(mem)
    except Exception:
        return None, None


class GpuTimeline:
    """Background nvidia-smi poller that ALSO prints a periodic util timeline."""

    def __init__(self, sample=1.0, print_every=5.0):
        self.sample = sample
        self.print_every = print_every
        self._stop = threading.Event()
        self._t = None
        self.utils = []
        self.mems = []

    def __enter__(self):
        self._t = threading.Thread(target=self._run, daemon=True)
        self._t.start()
        return self

    def _run(self):
        t0 = time.time()
        last_print = 0.0
        while not self._stop.is_set():
            u, m = query_gpu()
            if u is not None:
                self.utils.append(u)
                self.mems.append(m)
                elapsed = time.time() - t0
                if elapsed - last_print >= self.print_every:
                    print(f"      [t+{elapsed:6.0f}s] GPU util {u:5.0f}%   "
                          f"mem {m:7.0f} MiB", flush=True)
                    last_print = elapsed
            time.sleep(self.sample)

    def __exit__(self, *exc):
        self._stop.set()
        if self._t:
            self._t.join(timeout=3)

    def summary(self):
        if not self.utils:
            print("    no nvidia-smi samples collected")
            return
        n = len(self.utils)
        busy = sum(1 for u in self.utils if u >= 20)
        print(f"    util peak {max(self.utils):.0f}%, mean "
              f"{sum(self.utils)/n:.0f}%, "
              f"fraction of samples >=20% util: {100*busy/n:.0f}%   "
              f"mem peak {max(self.mems):.0f} MiB")
        print(f"    (a low busy-fraction over a long run = long CPU-only "
              f"stretches = what an idle-GPU reaper kills)")


def introspect_popv():
    print("\n[A] popV API introspection")
    import popv
    popv.settings.cuml = True
    popv.settings.n_jobs = 16
    print(f"    popv.__version__ = {getattr(popv, '__version__', 'unknown')}")

    # settings
    try:
        s = popv.settings
        attrs = [a for a in dir(s) if not a.startswith('_')]
        print(f"    popv.settings attrs: {attrs}")
        for a in attrs:
            if any(k in a.lower() for k in ("accel", "device", "gpu", "cuda", "jobs", "worker")):
                try:
                    print(f"      popv.settings.{a} = {getattr(s, a)!r}")
                except Exception:
                    pass
    except Exception as e:
        print(f"    popv.settings not available ({e})")

    # annotate_data signature
    try:
        sig = inspect.signature(popv.hub.HubModel.annotate_data)
        print(f"    HubModel.annotate_data{sig}")
        hits = [p for p in sig.parameters
                if any(k in p.lower() for k in
                       ("accel", "device", "gpu", "method", "n_jobs"))]
        print(f"    -> device/method-relevant params: {hits or 'NONE'}")
    except Exception as e:
        print(f"    could not inspect annotate_data ({e})")


def probe_popv(h5ad, cache_dir, n_sub, batch_key, gene_symbol_key):
    print(f"\n[B] popV on a {n_sub}-bead subset, with live GPU timeline")
    import numpy as np
    import scanpy as sc
    import popv

    print(f"    loading {h5ad}")
    adata = sc.read_h5ad(h5ad)
    print(f"    full: {adata.n_obs:,} x {adata.n_vars:,}")

    # stratified-ish subsample by puck if present
    if n_sub < adata.n_obs:
        rng = np.random.default_rng(0)
        if batch_key in adata.obs.columns:
            idx = []
            pucks = adata.obs[batch_key].astype(str).unique()
            per = max(1, n_sub // len(pucks))
            for p in pucks:
                pool = np.where(adata.obs[batch_key].astype(str).values == p)[0]
                idx.extend(rng.choice(pool, size=min(per, len(pool)), replace=False))
            idx = np.array(idx)
        else:
            idx = rng.choice(adata.n_obs, size=n_sub, replace=False)
        adata = adata[idx].copy()
    print(f"    subset: {adata.n_obs:,} x {adata.n_vars:,}")

    if gene_symbol_key not in adata.var.columns:
        adata.var[gene_symbol_key] = adata.var_names

    os.makedirs(cache_dir, exist_ok=True)
    print(f"    pulling model (cache: {cache_dir})")
    hmo = popv.hub.HubModel.pull_from_huggingface_hub(
        "popV/tabula_sapiens_All_Cells", cache_dir=cache_dir
    )

    print(f"    running annotate_data (watch the timeline; per-member lightning")
    print(f"    banners will print 'GPU available ... used: True/False')")
    with GpuTimeline(sample=1.0, print_every=5.0) as tl:
        t0 = time.time()
        annotated = hmo.annotate_data(
            adata,
            query_batch_key=batch_key,
            prediction_mode="inference",
            gene_symbols=gene_symbol_key,
        )
        dur = time.time() - t0
        tl.summary()
    print(f"    annotate_data wall time: {dur:.0f}s for {adata.n_obs:,} beads")
    if 'popv_prediction' in annotated.obs.columns:
        print(f"    popv_prediction types: "
              f"{annotated.obs['popv_prediction'].nunique()}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default=(
        "/work/sdz852/WORKING/slide-TCR-seq-working/data/outputs/03_qc/"
        "all_pucks_merged_QC.h5ad"))
    ap.add_argument("--cache", default=(
        "/work/sdz852/WORKING/slide-TCR-seq-working/data/outputs/04_annotation/"
        "popv_cache"))
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--batch-key", default="puck_id")
    ap.add_argument("--gene-symbol-key", default="feature_name")
    args = ap.parse_args()

    print("=" * 70)
    print("  popV GPU engagement probe")
    print(f"  host: {os.uname().nodename}")
    print("=" * 70)

    introspect_popv()

    if not os.path.exists(args.h5ad):
        print(f"\n[B] SKIP: merged QC not found at {args.h5ad}")
        print("    pass --h5ad to point at it; [A] above is still informative.")
        return
    probe_popv(args.h5ad, args.cache, args.n, args.batch_key, args.gene_symbol_key)

    print("\n" + "=" * 70)
    print("  Read [A] for the device/method knobs popV actually exposes.")
    print("  Read [B]'s busy-fraction + timeline: if the GPU is pinned the")
    print("  whole run, the full job is fine once logs stream. If there are")
    print("  long ~0% stretches, either restrict popV's methods to the fast/")
    print("  GPU members (if [A] shows a methods= knob) or run popV off-GPU.")
    print("=" * 70)


if __name__ == "__main__":
    main()
