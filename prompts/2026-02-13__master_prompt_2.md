# Slide-TCR-seq Spatial Transcriptomics Analysis
## Master Project Prompt - Comprehensive Record

**Project:** HPV+ HNSCC Tumor-Immune Microenvironment Characterization  
**PI:** Jake Lehle (Texas Biomedical Research Institute)  
**Collaborator:** Dr. Sophia Liu (Ragon Institute)  
**Timeline:** March - June 2026  
**Target:** K99/R00 Preliminary Data  

---

# Executive Summary

This project performs comprehensive spatial transcriptomics analysis of HPV+ head and neck 
squamous cell carcinoma (HNSCC) Slide-TCR-seq data to identify functional TCR-neoantigen 
pairs and characterize the tumor-immune microenvironment.

## Primary Research Questions
1. Do mutation-bearing spots co-localize with expanded T cell clones?
2. Are predicted neoantigens expressed in tumor regions where T cells infiltrate?
3. Is there evidence of T cell exclusion from neoantigen-rich areas?
4. Can we identify candidate TCR-neoantigen pairs for therapeutic validation?

## Central Hypothesis
APOBEC3B-driven mutagenesis creates neoantigens that should be recognized by 
tumor-infiltrating CD8+ T cells, yet these tumors maintain immune evasion. 
Understanding the spatial relationships between mutated cells, viral infection 
status, and T cell infiltration may reveal mechanisms of immune escape.

---

# Project Infrastructure

## Directory Structure
```
slide-TCR-seq-working/
├── data/
│   ├── inputs/
│   │   ├── fastq/
│   │   │   ├── 2022-01-28_Puck_211214_29/
│   │   │   ├── 2022-01-28_Puck_211214_37/
│   │   │   └── 2022-01-28_Puck_211214_40/
│   │   └── ref/GRCh38/
│   └── outputs/
│       ├── analysis/
│       │   ├── anndata/
│       │   ├── annotated/
│       │   ├── figures/
│       │   ├── fusions/
│       │   ├── integrated/
│       │   ├── mutations/
│       │   ├── neoantigens/
│       │   ├── processed/
│       │   ├── tcr/
│       │   ├── tls/
│       │   └── viral/
│       └── working/
├── scripts/
│   └── example_reference_scripts/
├── prompts/
├── reports/
├── envs/
└── logs/
```

## Starting Data (From Completed QC Pipeline)
| File | Location | Description |
|------|----------|-------------|
| all_pucks_processed.h5ad | data/outputs/analysis/processed/ | Merged, QC'd AnnData (~90K cells) |
| Puck_211214_29_processed.h5ad | data/outputs/analysis/processed/ | Individual puck (~30K cells) |
| Puck_211214_37_processed.h5ad | data/outputs/analysis/processed/ | Individual puck (~30K cells) |
| Puck_211214_40_processed.h5ad | data/outputs/analysis/processed/ | Individual puck (~30K cells) |
| Puck_*.matched.bam | data/inputs/fastq/2022-01-28_Puck_*/ | Aligned BAM files |
| GRCh38.primary_assembly.genome.fa | data/inputs/ref/GRCh38/ | Reference genome |
| gencode.v49.primary_assembly.annotation.gtf | data/inputs/ref/GRCh38/ | Gene annotations |

## Reference Scripts (Adapt for Spatial Data)
| Script | Purpose |
|--------|---------|
| SC_Cluster_Annotation.py | popV annotation + cluster assignment |
| SComatic_script.py | Somatic mutation calling |
| Kraken2_script.py | Viral sequence detection |
| 2025-12-20-COSMIC_Signatures.py | Mutational signature analysis |
| 2025-12-19_Virus_Porcessing.py | Viral read processing |

## Technical Constraints
- ~30K cells per puck after QC filtering
- Read depth may limit single-cell mutation calling - aggregate by cell type first
- Spatial coordinates from original Slide-seq bead positions
- BAM cell barcodes: XB tag (matched files) or XC tag (all_illumina)
- Reference genome: GRCh38 with GENCODE v49

---

# Phase Overview

The analysis is divided into 6 sequential phases, each designed to complete in 
6-12 hours with 5-12 subtasks for optimal AGI pipeline performance.

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          PHASE EXECUTION FLOW                                │
└─────────────────────────────────────────────────────────────────────────────┘

    ┌─────────────┐
    │   PHASE 1   │  ◄── BLOCKING: Must complete before any other phase
    │  Cell Type  │      
    │ Annotation  │      Output: celltype_annotated.h5ad
    └──────┬──────┘
           │
           ▼
    ┌──────┴──────┐
    │             │
    ▼             ▼
┌─────────┐   ┌─────────┐
│ PHASE 2 │   │ PHASE 3 │  ◄── Can run in PARALLEL after P1
│Spatial &│   │Viral &  │
│  TCR    │   │Mutation │
└────┬────┘   └────┬────┘
     │             │
     │             ▼
     │       ┌─────────┐
     │       │ PHASE 4 │  ◄── Requires P3 outputs
     │       │Neoantigen│
     │       │Prediction│
     │       └────┬────┘
     │             │
     └──────┬──────┘
            │
            ▼
     ┌─────────────┐
     │   PHASE 5   │  ◄── Requires P2 + P3 + P4 outputs
     │ Integration │
     │ & Analysis  │
     └──────┬──────┘
            │
            ▼
     ┌─────────────┐
     │   PHASE 6   │  ◄── Requires ALL previous phases
     │ Publication │
     │ & Synthesis │
     └─────────────┘
```

---

# Phase Definitions

## PHASE 1: Cell Type Annotation (BLOCKING)
**Priority:** CRITICAL - All downstream phases depend on this  
**Timeline:** March 18-20, 2026  
**Estimated Subtasks:** 5-8  

### Tasks Included
| Task ID | Name | Description |
|---------|------|-------------|
| 1.3 | Cell Type Annotation | popV with Tabula Sapiens + cluster-based assignment |

### Input Files
- `data/outputs/analysis/processed/all_pucks_processed.h5ad`

### Output Files
- `data/outputs/analysis/annotated/all_pucks_celltype_annotated.h5ad`
- `data/outputs/analysis/figures/celltype_umap.pdf`
- `data/outputs/analysis/figures/celltype_spatial_*.pdf`
- `reports/celltype_annotation_summary.md`

### Key Deliverables
- adata.obs columns: `popv_prediction`, `popv_prediction_score`, `clusters`, `final_annotation`
- Cell type composition tables per puck
- Marker gene validation plots

### Tools & Packages
- Python: popV, scanpy, squidpy, bbknn
- Model: `popV/tabula_sapiens_All_Cells` (HuggingFace)

### Success Criteria
- >80% cells with popv_prediction_score > 0.5
- All expected cell types detected (T cells, B cells, myeloid, epithelial, fibroblasts)
- Spatial plots showing biologically plausible distributions

---

## PHASE 2: Spatial Analysis & TCR Foundation
**Priority:** HIGH  
**Timeline:** March 25-31, 2026  
**Estimated Subtasks:** 10-14  
**Dependencies:** PHASE 1 complete  

### Tasks Included
| Task ID | Name | Description |
|---------|------|-------------|
| 1.4 | Spatial Neighborhood Analysis | Neighborhood enrichment, co-occurrence |
| 1.5 | TCR Repertoire Characterization | MiXCR clonotype calling, diversity metrics |
| 1.6 | Spatially Variable Genes & LR | Moran's I, ligand-receptor analysis |
| 4.1* | TLS Identification & Characterization | Full TLS analysis (merged from duplicate) |

*Note: Task 4.1 merged here to consolidate all spatial landmark identification

### Input Files
- `data/outputs/analysis/annotated/all_pucks_celltype_annotated.h5ad` (from P1)
- `data/inputs/fastq/2022-01-28_Puck_*/Puck_*.matched.bam`
- FASTQ files for TCR calling (Read 2 contains CDR3)

### Output Files
- `data/outputs/analysis/tcr/clonotype_catalog.tsv`
- `data/outputs/analysis/tcr/diversity_metrics.csv`
- `data/outputs/analysis/tls/tls_characterization.csv`
- `data/outputs/analysis/figures/nhood_enrichment.pdf`
- `data/outputs/analysis/figures/spatial_tcr_map.pdf`
- `data/outputs/analysis/figures/tls_spatial_map.pdf`
- `data/outputs/analysis/figures/ligrec_heatmap.pdf`

### Key Deliverables
- TCR clonotypes linked to spatial barcodes
- adata.obs columns: `has_tcr`, `clonotype_id`, `clone_size`, `is_expanded`
- adata.obs column: `tls_status`, `tls_region_id`
- Neighborhood enrichment matrices
- Spatially variable gene list
- Ligand-receptor interaction predictions

### Tools & Packages
- Python: scanpy, squidpy, scipy.stats
- External: MiXCR v4.x
- Reference: https://github.com/soph-liu/Slide-TCR-seq

### Success Criteria
- TCR detected in >5% of T cell spots
- TLS regions identified with B/T cell co-localization
- Significant L-R pairs identified (PD-L1/PD-1, CXCL13/CXCR5)

---

## PHASE 3: Viral & Mutation Detection
**Priority:** HIGH  
**Timeline:** April 7-24, 2026  
**Estimated Subtasks:** 10-12  
**Dependencies:** PHASE 1 complete (needs cell type labels for SComatic)  
**Note:** Can run in PARALLEL with PHASE 2  

### Tasks Included
| Task ID | Name | Description |
|---------|------|-------------|
| 2.1 | HPV Detection | Kraken2 viral detection, load quantification |
| 2.2 | Somatic Mutation Calling | SComatic + APOBEC signature analysis |
| 2.3 | Fusion Detection | Arriba + STAR-Fusion consensus calling |

### Input Files
- `data/outputs/analysis/annotated/all_pucks_celltype_annotated.h5ad` (from P1)
- `data/inputs/fastq/2022-01-28_Puck_*/Puck_*.matched.bam`
- `data/inputs/ref/GRCh38/` (reference genome)

### Output Files
- `data/outputs/analysis/viral/hpv_status_annotation.csv`
- `data/outputs/analysis/viral/hpv_spatial_heatmap.pdf`
- `data/outputs/analysis/mutations/scomatic_filtered.vcf`
- `data/outputs/analysis/mutations/apobec_scores.csv`
- `data/outputs/analysis/mutations/tmb_by_region.csv`
- `data/outputs/analysis/fusions/consensus_fusions.tsv`

### Key Deliverables
- adata.obs columns: `hpv_status`, `hpv_viral_load`
- Filtered somatic variant catalog (VCF)
- APOBEC signature enrichment scores
- Tumor mutational burden (TMB) per region
- High-confidence fusion catalog

### Tools & Packages
- Python: pysam, pandas, squidpy
- External: Kraken2 + Bracken, SComatic, VEP, STAR-Fusion, Arriba
- Databases: gnomAD, COSMIC signatures

### Success Criteria
- HPV+ regions clearly defined (E6/E7 expression)
- Somatic variants detected in tumor regions
- APOBEC signature enrichment in HPV+ regions
- Consensus fusions (detected by both tools)

---

## PHASE 4: Neoantigen Prediction
**Priority:** MEDIUM-HIGH  
**Timeline:** May 1-19, 2026  
**Estimated Subtasks:** 10-12  
**Dependencies:** PHASE 3 complete (needs VCF, fusions)  

### Tasks Included
| Task ID | Name | Description |
|---------|------|-------------|
| 3.1 | SNV Neoantigen Prediction | pVACseq MHC binding prediction |
| 3.2+3.3 | Comprehensive Database Cross-Reference | TSNAdb, NEPdb, VDJdb, ANNOVAR, COSMIC |
| 3.4 | Fusion Neoantigen Prediction | Junction peptide MHC binding |
| 3.5 | HLA Binding for Novel Variants | NetMHCpan de novo prediction |

*Note: Tasks 3.2 and 3.3 combined into comprehensive database annotation step

### Input Files
- `data/outputs/analysis/mutations/scomatic_filtered.vcf` (from P3)
- `data/outputs/analysis/fusions/consensus_fusions.tsv` (from P3)
- `data/outputs/analysis/tcr/clonotype_catalog.tsv` (from P2)
- HLA types (infer from RNA-seq or request from Sophia)

### Output Files
- `data/outputs/analysis/neoantigens/pvacseq_results.tsv`
- `data/outputs/analysis/neoantigens/database_matches.csv`
- `data/outputs/analysis/neoantigens/fusion_neoantigens.tsv`
- `data/outputs/analysis/neoantigens/hla_binding_predictions.tsv`
- `data/outputs/analysis/neoantigens/annotated_variants.tsv`

### Key Deliverables
- Predicted neoantigen catalog with binding affinities
- Validated neoantigens from TSNAdb/NEPdb
- TCR-epitope matches from VDJdb
- COSMIC/gnomAD annotated variants
- Neoantigen load per spot annotation

### Tools & Packages
- Python: pandas, requests (API queries)
- External: pVACseq, NetMHCpan 4.1, ANNOVAR
- Databases: TSNAdb v2.0, NEPdb, VDJdb, COSMIC, gnomAD, dbSNP

### Success Criteria
- Neoantigens predicted with IC50 < 500 nM
- Cross-reference matches found in public databases
- HLA typing successful (or common alleles used)

---

## PHASE 5: Integration & Synthesis
**Priority:** MEDIUM  
**Timeline:** April 28-30, 2026 (can start after P2+P3) + May continuation  
**Estimated Subtasks:** 8-10  
**Dependencies:** PHASE 2, 3, and 4 complete  

### Tasks Included
| Task ID | Name | Description |
|---------|------|-------------|
| 2.4 | Mutation-TCR Integration | Spatial co-localization analysis |
| 2.5 | Clonal Expansion Analysis | TCR diversity, expansion patterns |

### Input Files
- `data/outputs/analysis/annotated/all_pucks_celltype_annotated.h5ad`
- `data/outputs/analysis/tcr/clonotype_catalog.tsv`
- `data/outputs/analysis/mutations/scomatic_filtered.vcf`
- `data/outputs/analysis/neoantigens/pvacseq_results.tsv`
- `data/outputs/analysis/tls/tls_characterization.csv`

### Output Files
- `data/outputs/analysis/integrated/mutation_tcr_colocalization.csv`
- `data/outputs/analysis/integrated/tcr_neoantigen_pairs.csv`
- `data/outputs/analysis/integrated/clonal_expansion_patterns.csv`
- `data/outputs/analysis/figures/mutation_tcr_spatial.pdf`
- `data/outputs/analysis/figures/clone_distribution.pdf`

### Key Deliverables
- Distance metrics: mutation spots to TCR+ spots
- Candidate TCR-neoantigen pairs (prioritized)
- Clonal expansion in tumor vs stroma vs TLS
- Spatial clustering of expanded clones

### Tools & Packages
- Python: scanpy, squidpy, scipy.stats, pandas

### Success Criteria
- At least 3 candidate TCR-neoantigen pairs identified
- Clear spatial patterns of clonal expansion
- Statistical significance of co-localization

---

## PHASE 6: Publication & Final Synthesis
**Priority:** MEDIUM  
**Timeline:** May 29 - June 10, 2026  
**Estimated Subtasks:** 6-8  
**Dependencies:** ALL previous phases complete  

### Tasks Included
| Task ID | Name | Description |
|---------|------|-------------|
| 4.2 | Comprehensive Spatial Correlation | Statistical analysis of all spatial relationships |
| 4.3 | Publication Figures & Manuscript | Final figures, manuscript draft |

### Input Files
- ALL outputs from Phases 1-5

### Output Files
- `data/outputs/analysis/figures/main_figure_*.pdf` (8 main figures)
- `data/outputs/analysis/figures/supp_figure_*.pdf` (4-6 supplemental)
- `reports/spatial_correlation_report.md`
- `reports/manuscript_draft.docx`
- `reports/figure_legends.docx`
- `reports/executive_summary.md`

### Main Figures Planned
1. Study overview & QC
2. HPV detection & viral landscape
3. Mutation landscape & APOBEC signature
4. Neoantigen prediction & distribution
5. TCR repertoire & clonality
6. TLS characterization
7. Spatial correlations (key findings)
8. HPV+ vs HPV- comparison

### Key Deliverables
- Publication-ready figures (300 dpi PDF)
- Correlation matrices with p-values
- Statistical test results
- Manuscript draft (~4000 words)
- GitHub repo with reproducible code

### Tools & Packages
- Python: matplotlib, seaborn, scanpy, squidpy
- statsmodels for regression

### Success Criteria
- All main figures generated
- Statistical analyses complete
- Manuscript outline filled

---

# Task Dependency Matrix

| Task | Depends On | Blocks |
|------|------------|--------|
| 1.3 Cell Types | - | ALL other tasks |
| 1.4 Spatial Nhood | 1.3 | 4.2, 4.3 |
| 1.5 TCR | 1.3 (for cell type context) | 2.4, 2.5, 3.2 |
| 1.6 SVG/LR | 1.3 | 4.2, 4.3 |
| 4.1 TLS | 1.3, 1.4 | 2.5, 4.2, 4.3 |
| 2.1 HPV | 1.3 | 2.2, 4.2 |
| 2.2 SComatic | 1.3 | 2.3, 3.1, 3.2+3.3 |
| 2.3 Fusions | 2.2 | 3.4 |
| 3.1 pVACseq | 2.2 | 3.2+3.3, 3.5 |
| 3.2+3.3 DB CrossRef | 2.2, 3.1, 1.5 | 2.4, 4.2 |
| 3.4 Fusion Neo | 2.3 | 4.2 |
| 3.5 HLA Binding | 3.1 | 4.2 |
| 2.4 Mut-TCR | 1.5, 2.2, 3.1 | 4.2 |
| 2.5 Clonal Exp | 1.5, 4.1 | 4.2 |
| 4.2 Spatial Corr | ALL P1-P5 | 4.3 |
| 4.3 Publication | 4.2 | - |

---

# Timeline Summary

| Phase | Tasks | Start | End | Duration |
|-------|-------|-------|-----|----------|
| P1 | 1.3 | Mar 18 | Mar 20 | 3 days |
| P2 | 1.4, 1.5, 1.6, 4.1 | Mar 21 | Mar 31 | 11 days |
| P3 | 2.1, 2.2, 2.3 | Mar 21* | Apr 24 | ~5 weeks |
| P4 | 3.1, 3.2+3.3, 3.4, 3.5 | Apr 25 | May 19 | ~4 weeks |
| P5 | 2.4, 2.5 | May 20 | May 28 | 9 days |
| P6 | 4.2, 4.3 | May 29 | Jun 10 | 13 days |

*P2 and P3 can start simultaneously after P1 completes

---

# Critical Path

The minimum time to completion follows this critical path:

```
P1 (3 days) → P3 (5 weeks) → P4 (4 weeks) → P5 (9 days) → P6 (13 days)
```

Total critical path: ~14 weeks (March 18 - June 10)

P2 runs in parallel with P3, so it's not on the critical path unless it takes longer.

---

# Checkpoint Files

Each phase produces checkpoint files that enable restart and verify completion:

| Phase | Checkpoint File | Verification |
|-------|----------------|--------------|
| P1 | `data/outputs/analysis/annotated/all_pucks_celltype_annotated.h5ad` | Contains `final_annotation` column |
| P2 | `data/outputs/analysis/tcr/clonotype_catalog.tsv` | >100 clonotypes detected |
| P2 | `data/outputs/analysis/tls/tls_characterization.csv` | TLS regions identified |
| P3 | `data/outputs/analysis/mutations/scomatic_filtered.vcf` | >100 somatic variants |
| P3 | `data/outputs/analysis/viral/hpv_status_annotation.csv` | HPV status per spot |
| P4 | `data/outputs/analysis/neoantigens/pvacseq_results.tsv` | Predicted neoantigens |
| P5 | `data/outputs/analysis/integrated/tcr_neoantigen_pairs.csv` | Candidate pairs |
| P6 | `data/outputs/analysis/figures/main_figure_1.pdf` | Figures generated |

---

# Environment Specifications

## Core Python Environment
```yaml
name: slide_tcr_seq
channels:
  - conda-forge
  - bioconda
dependencies:
  - python>=3.10
  - scanpy>=1.9.0
  - squidpy>=1.2.0
  - anndata>=0.9.0
  - popv
  - bbknn
  - pandas
  - numpy
  - scipy
  - matplotlib
  - seaborn
  - pysam
  - statsmodels
```

## External Tools (Separate Installation)
- MiXCR v4.x
- Kraken2 + Bracken
- SComatic
- VEP (Ensembl Variant Effect Predictor)
- STAR-Fusion v1.13+
- Arriba v2.4+
- pVACseq
- NetMHCpan 4.1
- ANNOVAR

---

# Collaborators & Contact

| Name | Institution | Role |
|------|-------------|------|
| Jake Lehle | Texas Biomed | Lead analyst |
| Dr. Sophia Liu | Ragon Institute | Data provider |
| Dr. Diako Ebrahimi | Texas Biomed | Primary mentor |
| Dr. Reuben Harris | UT Health San Antonio | APOBEC expertise |

**GitHub Repository:** https://github.com/JakeLehle/HNSCC_slide-TCR-seq

---

# Notes for AGI Pipeline Execution

1. **CRITICAL**: Phase 1 MUST complete successfully before any other phase starts
2. Reference scripts in `scripts/example_reference_scripts/` contain working code - READ THEM
3. All tools are Python-based unless explicitly noted - DO NOT use R
4. Checkpoint intermediate results to enable restart on failure
5. Save all generated scripts to `scripts/` directory
6. Web search may be needed for tool documentation
7. Each phase prompt should be self-contained with explicit file paths

---

# Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-02-13 | Initial comprehensive master prompt |
| 1.1 | 2026-02-13 | Reorganized into 6 phases, merged duplicate TLS task, combined 3.2+3.3 |
