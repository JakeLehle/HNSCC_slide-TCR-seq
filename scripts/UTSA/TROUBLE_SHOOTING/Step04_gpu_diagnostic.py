#!/usr/bin/env python3
"""
Step04_gpu_diagnostic.py

Standalone GPU-engagement diagnostic for the Step04 annotation pipeline. Run
INTERACTIVELY on an a100 node to confirm that each layer that is supposed to
use the GPU actually does, BEFORE committing a multi-day batch job.

    srun --partition=gpu1a100 --gres=gpu:1 -N 1 -n 1 -c 10 --time=00:30:00 --pty bash
    module load anaconda3
    conda activate spatial
    python Step04_gpu_diagnostic.py

Three escalating checks; the GPU-bound ones poll nvidia-smi in the background
so you can SEE utilization rise without opening a second terminal:

  [1] torch sees CUDA          version, cuda build, is_available, device name,
                               CUDA_VISIBLE_DEVICES
  [2] torch computes on GPU    a ~20s matmul burn; reports peak/mean util + mem
  [3] scvi trains on GPU       a tiny RegressionModel with accelerator='gpu'
                               (the EXACT path cell2location's reference
                               regression uses); also reports which device the
                               model parameters ended up on, and falls back to
                               use_gpu=True if the accelerator= kwarg is rejected

Why this exists: popV is a CPU-heavy ensemble and will NOT keep the GPU busy.
If the real job spends its first hours in popV while holding an idle a100, a
cluster GPU-idle reaper can cancel it before cell2location (the GPU-heavy
stage) starts. This script isolates "can the GPU be used at all" from that
scheduling problem, so we don't confuse the two.
"""

import os
import sys
import time
import threading
import subprocess


def query_gpu():
    """Return (util_percent, mem_used_MiB) for the first visible GPU, or (None, None)."""
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


class GpuMonitor:
    """Poll nvidia-smi in a background thread; report peak/mean util + peak mem."""

    def __init__(self, interval=0.5):
        self.interval = interval
        self._stop = threading.Event()
        self._t = None
        self.utils = []
        self.mems = []

    def __enter__(self):
        self._t = threading.Thread(target=self._run, daemon=True)
        self._t.start()
        return self

    def _run(self):
        while not self._stop.is_set():
            u, m = query_gpu()
            if u is not None:
                self.utils.append(u)
                self.mems.append(m)
            time.sleep(self.interval)

    def __exit__(self, *exc):
        self._stop.set()
        if self._t:
            self._t.join(timeout=2)

    def report(self, label):
        if not self.utils:
            print(f"    [{label}] nvidia-smi returned no samples "
                  f"(is it on PATH on this node?)")
            return
        print(f"    [{label}] GPU util: peak {max(self.utils):.0f}%, "
              f"mean {sum(self.utils) / len(self.utils):.0f}%   "
              f"mem: peak {max(self.mems):.0f} MiB")


def check_torch():
    print("\n[1] torch / CUDA visibility")
    try:
        import torch
    except Exception as e:
        print(f"    FAIL: could not import torch ({e})")
        return False
    print(f"    torch                {torch.__version__}")
    print(f"    torch.version.cuda   {torch.version.cuda}")
    print(f"    cuda.is_available()  {torch.cuda.is_available()}")
    print(f"    CUDA_VISIBLE_DEVICES {os.environ.get('CUDA_VISIBLE_DEVICES')!r}")
    if torch.cuda.is_available():
        print(f"    device               {torch.cuda.get_device_name(0)}")
        return True
    print("    -> torch cannot see a GPU. If nvidia-smi shows a card on this "
          "node,\n       the env likely has a CPU-only torch wheel "
          "(torch.version.cuda is None)\n       or CUDA_VISIBLE_DEVICES is "
          "empty. Fix before anything else.")
    return False


def check_matmul(seconds=20):
    print(f"\n[2] torch matmul burn on GPU (~{seconds}s; watch util climb)")
    import torch
    if not torch.cuda.is_available():
        print("    SKIP: no CUDA device")
        return False
    dev = torch.device("cuda")
    a = torch.randn((8192, 8192), device=dev)
    b = torch.randn((8192, 8192), device=dev)
    n = 0
    with GpuMonitor() as mon:
        t0 = time.time()
        while time.time() - t0 < seconds:
            c = a @ b
            torch.cuda.synchronize()
            n += 1
        mon.report("matmul")
    print(f"    completed {n} matmuls in ~{seconds}s "
          f"(peak util well below ~80% would be suspicious)")
    return True


def check_scvi(epochs=20):
    print(f"\n[3] scvi RegressionModel on GPU (accelerator='gpu', {epochs} epochs)")
    print("    this is the exact code path cell2location's reference regression uses")
    try:
        import numpy as np
        import pandas as pd
        import anndata as ad
        from cell2location.models import RegressionModel
    except Exception as e:
        print(f"    SKIP: could not import cell2location/scvi ({e})")
        return False
    import torch
    if not torch.cuda.is_available():
        print("    SKIP: no CUDA device")
        return False

    rng = np.random.default_rng(0)
    n_cells, n_genes = 2000, 500
    X = rng.poisson(0.5, size=(n_cells, n_genes)).astype("float32")
    obs = pd.DataFrame({
        "batch": pd.Categorical(rng.integers(0, 3, n_cells).astype(str)),
        "cell_type": pd.Categorical(rng.integers(0, 5, n_cells).astype(str)),
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{j}" for j in range(n_genes)])
    a = ad.AnnData(X=X, obs=obs, var=var)

    RegressionModel.setup_anndata(adata=a, batch_key="batch", labels_key="cell_type")
    mod = RegressionModel(a)

    used = None
    with GpuMonitor() as mon:
        try:
            mod.train(max_epochs=epochs, accelerator="gpu")
            used = "accelerator='gpu'"
        except TypeError as e:
            print(f"    accelerator= kwarg rejected ({e}); retrying use_gpu=True")
            mod = RegressionModel(a)
            mod.train(max_epochs=epochs, use_gpu=True)
            used = "use_gpu=True"
        mon.report("scvi")

    try:
        dev = next(mod.module.parameters()).device
        print(f"    model parameters live on: {dev}")
    except Exception:
        pass
    print(f"    scvi training completed via {used}")
    print(f"    -> in the pipeline, use the kwarg form that worked here.")
    return True


def main():
    print("=" * 70)
    print("  Step04 GPU engagement diagnostic")
    print(f"  host: {os.uname().nodename}")
    print("=" * 70)

    ok_torch = check_torch()
    if ok_torch:
        check_matmul()
        check_scvi()

    print("\n" + "=" * 70)
    print("  VERDICT")
    if not ok_torch:
        print("  torch cannot see the GPU. Nothing downstream will use it.")
        print("  Fix the env's CUDA torch first (see check [1] notes).")
    else:
        print("  If [2] and [3] both showed high utilization and finished, the")
        print("  GPU is usable from this env and the cell2location stage WILL")
        print("  drive it. A 2h idle-then-cancel is then a scheduling issue:")
        print("  popV (CPU-heavy) holds an idle GPU and gets reaped before")
        print("  cell2location starts. The fix is to stop holding a GPU during")
        print("  popV (decouple popV onto a CPU node), not to change torch.")
    print("=" * 70)


if __name__ == "__main__":
    main()
