#!/usr/bin/env python3
"""
Step04d_panglaodb_setup_and_enum.py  (download + enumerate, mostly read-only)

Fetches the PanglaoDB marker table into data/inputs/panglaodb/ and enumerates
the human (Hs) organ / cell-type universe so we can build the PanglaoDB -> 8
compartment cross-map. Reports, per PanglaoDB cell type, how many of its markers
are measurable in the spatial bead panel.

Env: sc_pre. Run on a theia LOGIN node (needs internet for the single download).
Hardcoded TX Biomed paths; does NOT import spatial_config.
"""

import os
import urllib.request
import numpy as np
import pandas as pd
import scanpy as sc

PROJECT_ROOT = "/master/jlehle/WORKING/slide-TCR-seq-working"
INPUTS = os.path.join(PROJECT_ROOT, "data/inputs/panglaodb")
PG_URL = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
PG_LOCAL = os.path.join(INPUTS, "PanglaoDB_markers_27_Mar_2020.tsv.gz")
ANNOTATED = os.path.join(PROJECT_ROOT, "data/outputs/04_annotation/all_pucks_annotated.h5ad")

os.makedirs(INPUTS, exist_ok=True)

if not os.path.exists(PG_LOCAL):
    print(f"Downloading PanglaoDB -> {PG_LOCAL}")
    urllib.request.urlretrieve(PG_URL, PG_LOCAL)
    print("  download complete")
else:
    print(f"[cache] {PG_LOCAL}")

pg = pd.read_table(PG_LOCAL)
print(f"\nraw rows: {len(pg):,}")
print(f"columns: {list(pg.columns)}")

GENE = 'official gene symbol'
pg_hs = pg[pg['species'].astype(str).str.contains('Hs', na=False)].copy()
print(f"\nhuman (Hs) rows: {len(pg_hs):,}  "
      f"cell types: {pg_hs['cell type'].nunique()}  organs: {pg_hs['organ'].nunique()}")

# spatial panel genes (backed read of var_names only, cheap)
var_names = set()
if os.path.exists(ANNOTATED):
    ad = sc.read_h5ad(ANNOTATED, backed='r')
    var_names = set(map(str, ad.var_names))
    print(f"spatial panel genes: {len(var_names):,}")
else:
    print(f"WARNING: spatial object missing at {ANNOTATED}; n_in_spatial will be 0")

print("\n=== Hs organs (by marker-row count) ===")
for organ, n in pg_hs['organ'].value_counts().items():
    print(f"  {organ}: {n}")

has_canon = 'canonical marker' in pg_hs.columns
rows = []
for ct, sub in pg_hs.groupby('cell type'):
    organs = ",".join(sorted(sub['organ'].dropna().astype(str).unique()))
    genes = sub[GENE].dropna().astype(str).unique()
    in_sp = [g for g in genes if g in var_names]
    n_canon_in_sp = 0
    if has_canon:
        canon = sub[pd.to_numeric(sub['canonical marker'], errors='coerce') == 1][GENE]
        canon = canon.dropna().astype(str).unique()
        n_canon_in_sp = sum(g in var_names for g in canon)
    rows.append((ct, organs, len(genes), len(in_sp), n_canon_in_sp))

df = (pd.DataFrame(rows, columns=['cell_type', 'organs', 'n_markers',
                                  'n_in_spatial', 'n_canonical_in_spatial'])
      .sort_values(['organs', 'cell_type']))

pd.set_option('display.max_rows', None)
pd.set_option('display.width', 200)
print("\n=== Hs cell types ===")
print(df.to_string(index=False))

out = os.path.join(INPUTS, "panglaodb_Hs_celltype_enumeration.tsv")
df.to_csv(out, sep='\t', index=False)
print(f"\nsaved: {out}")
print("\nPaste the organ list and cell-type table back so we can build the cross-map.")
