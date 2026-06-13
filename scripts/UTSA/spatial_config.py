#!/usr/bin/env python3
"""
spatial_config.py

Centralized configuration for the Slide-TCR-seq spatial analysis pipeline.
All Python scripts in this project import paths and parameters from here.

Mirrors the role of network_config.py / network_config_SC.py in the NMF paper.

Project: HPV16+ HNSCC Spatial Transcriptomics
Collaborator: Sophia Liu, Ragon Institute
Author: Jake Lehle, Texas Biomedical Research Institute

UTSA MIGRATION NOTE
-------------------
Ported from TX Biomed to UTSA. The only structural change is PROJECT_ROOT
(everything project-internal hangs off it). The NMF cross-reference block
below points at a SEPARATE project tree that is not part of this repo or the
backup pulled to UTSA; see the loud comment there before running Steps 04/08/09.
"""

import os
import sys
from datetime import datetime

# =============================================================================
# PROJECT PATHS
# =============================================================================

PROJECT_ROOT = "/work/sdz852/WORKING/slide-TCR-seq-working"

# --- Input directories (read-only source data) ---
INPUT_DIR       = os.path.join(PROJECT_ROOT, "data", "inputs")
FASTQ_DIR       = os.path.join(INPUT_DIR, "fastq")
REF_DIR         = os.path.join(INPUT_DIR, "ref", "GRCh38")
STAR_INDEX      = os.path.join(REF_DIR, "star")
GENOME_FA       = os.path.join(REF_DIR, "GRCh38.primary_assembly.genome.fa")
GTF_FILE        = os.path.join(REF_DIR, "gencode.v49.primary_assembly.annotation.gtf")
BCL_DIR         = os.path.join(FASTQ_DIR, "220116_NB501164_1345_AHLGH2BGXK")

# --- Output directories (one per pipeline step) ---
OUTPUT_DIR              = os.path.join(PROJECT_ROOT, "data", "outputs")
DIR_01_ALIGNMENT        = os.path.join(OUTPUT_DIR, "01_alignment")
DIR_02_ANNDATA          = os.path.join(OUTPUT_DIR, "02_anndata")
DIR_03_QC               = os.path.join(OUTPUT_DIR, "03_qc")
DIR_04_ANNOTATION       = os.path.join(OUTPUT_DIR, "04_annotation")
DIR_05_HPV              = os.path.join(OUTPUT_DIR, "05_hpv")
DIR_06_MUTATIONS        = os.path.join(OUTPUT_DIR, "06_mutations")
DIR_07_SIGNATURES       = os.path.join(OUTPUT_DIR, "07_signatures")
DIR_08_NEOANTIGENS      = os.path.join(OUTPUT_DIR, "08_neoantigens")
DIR_09_NMF_CROSSREF     = os.path.join(OUTPUT_DIR, "09_nmf_crossref")
DIR_10_TCR              = os.path.join(OUTPUT_DIR, "10_tcr")
DIR_11_COLOCALIZATION   = os.path.join(OUTPUT_DIR, "11_colocalization")
DIR_FIGURES             = os.path.join(OUTPUT_DIR, "figures")
DIR_TROUBLESHOOTING     = os.path.join(OUTPUT_DIR, "TROUBLESHOOTING")

# --- NMF paper cross-reference paths ---
# !!! CONFIRM BEFORE STEPS 04 / 08 / 09 !!!
# The 2026 NMF paper is a SEPARATE project tree. It is NOT part of this repo
# and was NOT in the backup pulled to UTSA. The files referenced below
# (adata_final, three_group_assignments, proteome) plus the cell2location
# reference used by Step04 must be staged onto UTSA before those steps run.
# Steps 02 and 03 do not touch any of these, so this can wait.
# The value below is a best-guess mirror path, not a verified location.
NMF_PAPER_ROOT      = "/work/sdz852/WORKING/2026_NMF_PAPER"   # <-- verify / stage files
NMF_NEOANTIGEN_SCRIPTS = os.path.join(NMF_PAPER_ROOT, "scripts", "NEOANTIGEN")
NMF_NEOANTIGEN_DATA    = os.path.join(NMF_PAPER_ROOT, "data", "FIG_7")
NMF_ADATA_FINAL        = os.path.join(NMF_PAPER_ROOT, "data", "FIG_3",
                                       "results_NMF_v0.1.1", "signatures",
                                       "adata_final.h5ad")
NMF_THREE_GROUPS       = os.path.join(NMF_PAPER_ROOT, "data", "FIG_4",
                                       "00_input", "three_group_assignments.tsv")
NMF_PROTEOME           = os.path.join(NMF_PAPER_ROOT, "data", "reference",
                                       "Homo_sapiens.GRCh38.pep.all.fa")

# --- cell2location single-cell reference (Step 04) ---
# The HNSCC single-cell object built from the NMF analysis
# (build_cell2location_reference.py): raw counts, obs['final_annotation']
# (12 types), batch column 'subject id'. Staged INTO this project so the UTSA
# pipeline is self-contained and does not depend on the NMF tree above. Step04
# imports this name; if it is ever missing, the script falls back to a stale
# path, so keep it defined here. Move the .h5ad to this location:
#   <PROJECT_ROOT>/data/inputs/annotation/reference_for_cell2location.h5ad
C2L_REFERENCE_PATH = os.path.join(INPUT_DIR, "annotation",
                                  "reference_for_cell2location.h5ad")

# =============================================================================
# PUCK DEFINITIONS
# =============================================================================

PUCKS = {
    "Puck_211214_29": {
        "sample_barcode": "AGATTTAA",
        "input_dir": "2022-01-28_Puck_211214_29",
        "description": "HPV16+ HNSCC section 1",
    },
    "Puck_211214_37": {
        "sample_barcode": "GGCGTCGA",
        "input_dir": "2022-01-28_Puck_211214_37",
        "description": "HPV16+ HNSCC section 2",
    },
    "Puck_211214_40": {
        "sample_barcode": "ATCACTCG",
        "input_dir": "2022-01-28_Puck_211214_40",
        "description": "HPV16+ HNSCC section 3 (best RNA recovery per Sophia)",
    },
}

PUCK_NAMES = list(PUCKS.keys())

# Derived paths for each puck's input data
def puck_input_dir(puck_name):
    """Full path to a puck's raw input directory."""
    return os.path.join(FASTQ_DIR, PUCKS[puck_name]["input_dir"])

def puck_barcode_matching(puck_name):
    """Path to the barcode_matching.txt.gz for a puck."""
    return os.path.join(puck_input_dir(puck_name), "barcode_matching",
                        f"{puck_name}_barcode_matching.txt.gz")

def puck_expression_file(puck_name):
    """Path to the matched digital expression file for a puck."""
    return os.path.join(puck_input_dir(puck_name),
                        f"{puck_name}.matched.digital_expression.txt.gz")

def puck_bam_file(puck_name):
    """Path to Sophia's matched BAM file for a puck."""
    return os.path.join(puck_input_dir(puck_name),
                        f"{puck_name}.matched.bam")

# =============================================================================
# ALIGNMENT PARAMETERS (Step 01)
# =============================================================================

CB_LEN    = 14       # Cell barcode length (Slide-seq, confirmed)
UMI_START = 15       # UMI start position in read
UMI_LEN   = 9        # UMI length (confirmed from Sophia's BAM XM tag)

# =============================================================================
# QC PARAMETERS (Step 03)
# =============================================================================

MIN_GENES         = 100      # Min genes per cell (lower threshold for spatial)
MIN_CELLS         = 10       # Min cells per gene
MT_THRESHOLD      = 25.0     # Max mitochondrial % (slightly higher for spatial)
N_HVG             = 2000     # Number of highly variable genes
LEIDEN_RESOLUTION = 1.0      # Clustering resolution
N_PCS             = 30       # PCA components for neighbors

# =============================================================================
# SCOMATIC PARAMETERS (Step 06)
# =============================================================================

SCOMATIC_MIN_VAF      = 0.05    # Minimum variant allele frequency
SCOMATIC_MIN_COVERAGE = 10      # Minimum read coverage
GNOMAD_MAX_AF         = 0.01    # Maximum gnomAD allele frequency (germline filter)

# =============================================================================
# NEOANTIGEN PARAMETERS (Step 08)
# =============================================================================

IC50_THRESHOLD    = 500      # nM, strong/intermediate binder cutoff
MIN_EXPRESSION    = 1.0      # TPM, minimum expression for neoantigen gene
ENSEMBL_RELEASE   = 115      # Match NMF paper proteome reference
OFFSET_SCAN_RANGE = 30       # SnpEff signal peptide offset (match NMF paper)

# =============================================================================
# FIGURE PARAMETERS
# =============================================================================

FONT_SIZE   = 30
RANDOM_SEED = 42
DPI         = 300

# Population colors (matching NMF paper color scheme)
COLOR_SBS2_HIGH = "#ed6a5a"   # coral
COLOR_CNV_HIGH  = "#F6D155"   # mustard
COLOR_NORMAL    = "#4682b4"   # steelblue
COLOR_OTHER     = "#e0e0e0"   # gray

# Puck colors (for multi-puck visualizations)
COLOR_PUCK_29 = "#7eb0d5"    # blue
COLOR_PUCK_37 = "#b2e061"    # green
COLOR_PUCK_40 = "#fd7f6f"    # salmon

# =============================================================================
# MATPLOTLIB DEFAULTS
# =============================================================================

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size':         FONT_SIZE,
    'axes.titlesize':    FONT_SIZE,
    'axes.labelsize':    FONT_SIZE,
    'xtick.labelsize':   FONT_SIZE - 4,
    'ytick.labelsize':   FONT_SIZE - 4,
    'legend.fontsize':   FONT_SIZE - 6,
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Arial', 'DejaVu Sans'],
    'savefig.dpi':       DPI,
    'savefig.bbox':      'tight',
})

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

REPORT_FILE = None   # Set by individual scripts if they want file logging

def banner(title, char="="):
    """Print a section banner."""
    log("")
    log(char * 80)
    log(f"  {title}")
    log(char * 80)

def log(msg=""):
    """Print a timestamped log message."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    if REPORT_FILE is not None:
        REPORT_FILE.write(line + "\n")

def ensure_dir(path):
    """Create directory if it does not exist. Returns the path."""
    os.makedirs(path, exist_ok=True)
    return path

def save_fig(fig, name, output_dir=None):
    """Save figure as PDF and PNG at 300 DPI. No figure titles over panels."""
    if output_dir is None:
        output_dir = DIR_FIGURES
    ensure_dir(output_dir)
    for ext in ['pdf', 'png']:
        fig.savefig(os.path.join(output_dir, f"{name}.{ext}"),
                    dpi=DPI, bbox_inches='tight')
    plt.close(fig)
    log(f"  Saved: {name}.pdf/.png -> {output_dir}")

# =============================================================================
# STARTUP SANITY CHECK
# =============================================================================
# Guard against importing a stale config from another machine. If PROJECT_ROOT
# does not exist on this host (e.g. a leftover /master/... path carried over to
# UTSA), fail loudly here with a clear message instead of dying deep inside
# os.makedirs several frames later with an opaque PermissionError.
if not os.path.isdir(PROJECT_ROOT):
    raise SystemExit(
        f"[spatial_config] PROJECT_ROOT does not exist on this host:\n"
        f"    {PROJECT_ROOT}\n"
        f"You are probably importing a stale config from another machine.\n"
        f"Confirm the UTSA copy (PROJECT_ROOT under /work/sdz852/...) is the\n"
        f"one beside the script on sys.path."
    )

print(f"[spatial_config] PROJECT_ROOT = {PROJECT_ROOT}", flush=True)
print(f"[spatial_config] OUTPUT_DIR   = {OUTPUT_DIR}", flush=True)
