#!/bin/bash

#SBATCH -J REF_GRCh38_v49                                              # Job name
#SBATCH -o /work/sdz852/WORKING/LOGS/GRCh38v49.o.log              # Name of the stdout output file
#SBATCH -e /work/sdz852/WORKING/LOGS/GRCh38v49.e.log              # Name of the stderr error file
#SBATCH -t 10-00:00:00                                                   # Time of job
#SBATCH -p compute2                                              # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80                                                    # Total number of cores 80 max

# Startup scripts

module load anaconda3
conda activate slide-TCR-seq

#===============================================================================
# DOWNLOAD AND BUILD GRCh38 STAR REFERENCE GENOME
#
# This script:
#   1. Creates the reference directory
#   2. Downloads GRCh38 genome FASTA from GENCODE
#   3. Downloads GENCODE gene annotation (GTF)
#   4. Builds STAR genome index
#
# Author: Jake Lehle
# Project: HPV+ HNSCC Spatial Transcriptomics
# Date: December 2024
#
# Prerequisites:
#   - STAR installed (latest version, e.g., 2.7.11b)
#   - wget or curl
#   - ~40GB disk space for reference files
#   - ~32GB RAM for index generation
#===============================================================================

set -e
set -o pipefail

#===============================================================================
# CONFIGURATION
#===============================================================================

# Reference directory
REF_DIR="/work/sdz852/WORKING/slide-TCR-seq/ref/GRCh38"

# GENCODE release (v44 is recent and well-supported)
GENCODE_RELEASE="49"

# URLs for GENCODE GRCh38 (primary assembly - no patches/haplotypes)
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"

# STAR parameters
THREADS=$(nproc)

# sjdbOverhang should be ReadLength - 1. For 42bp reads, use 41
# For flexibility with different read lengths, 100 is a good default
SJDB_OVERHANG=100

#===============================================================================
# LOGGING
#===============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

#===============================================================================
# STEP 1: CREATE DIRECTORY STRUCTURE
#===============================================================================

log "=============================================="
log "STEP 1: Creating directory structure"
log "=============================================="

mkdir -p "${REF_DIR}/downloads"
mkdir -p "${REF_DIR}/star"

log "Reference directory: ${REF_DIR}"
log "Downloads directory: ${REF_DIR}/downloads"
log "STAR index directory: ${REF_DIR}/star"

#===============================================================================
# STEP 2: CHECK STAR VERSION
#===============================================================================

log "=============================================="
log "STEP 2: Checking STAR version"
log "=============================================="

if ! command -v STAR &> /dev/null; then
    log "ERROR: STAR not found in PATH"
    log "Please install STAR first: conda install -c bioconda star"
    exit 1
fi

STAR_VERSION=$(STAR --version)
log "STAR version: ${STAR_VERSION}"

# Save version for reference
echo "${STAR_VERSION}" > "${REF_DIR}/star_version.txt"

#===============================================================================
# STEP 3: DOWNLOAD GENOME FASTA
#===============================================================================

log "=============================================="
log "STEP 3: Downloading GRCh38 genome FASTA"
log "=============================================="

GENOME_GZ="${REF_DIR}/downloads/GRCh38.primary_assembly.genome.fa.gz"
GENOME_FA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"

if [ -f "${GENOME_FA}" ]; then
    log "Genome FASTA already exists, skipping download..."
else
    if [ -f "${GENOME_GZ}" ]; then
        log "Compressed genome already downloaded, skipping download..."
    else
        log "Downloading from GENCODE release ${GENCODE_RELEASE}..."
        log "URL: ${GENOME_URL}"
        
        wget -c -O "${GENOME_GZ}" "${GENOME_URL}"
        
        if [ $? -ne 0 ]; then
            log "ERROR: Failed to download genome FASTA"
            exit 1
        fi
    fi
    
    log "Decompressing genome FASTA..."
    gunzip -c "${GENOME_GZ}" > "${GENOME_FA}"
    
    log "Genome FASTA ready: ${GENOME_FA}"
fi

# Verify file
GENOME_SIZE=$(du -h "${GENOME_FA}" | cut -f1)
log "Genome file size: ${GENOME_SIZE}"

#===============================================================================
# STEP 4: DOWNLOAD GTF ANNOTATION
#===============================================================================

log "=============================================="
log "STEP 4: Downloading GENCODE GTF annotation"
log "=============================================="

GTF_GZ="${REF_DIR}/downloads/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz"
GTF_FILE="${REF_DIR}/gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf"

if [ -f "${GTF_FILE}" ]; then
    log "GTF annotation already exists, skipping download..."
else
    if [ -f "${GTF_GZ}" ]; then
        log "Compressed GTF already downloaded, skipping download..."
    else
        log "Downloading from GENCODE release ${GENCODE_RELEASE}..."
        log "URL: ${GTF_URL}"
        
        wget -c -O "${GTF_GZ}" "${GTF_URL}"
        
        if [ $? -ne 0 ]; then
            log "ERROR: Failed to download GTF annotation"
            exit 1
        fi
    fi
    
    log "Decompressing GTF annotation..."
    gunzip -c "${GTF_GZ}" > "${GTF_FILE}"
    
    log "GTF annotation ready: ${GTF_FILE}"
fi

# Verify file
GTF_SIZE=$(du -h "${GTF_FILE}" | cut -f1)
GTF_GENES=$(grep -c "gene_type" "${GTF_FILE}" || echo "unknown")
log "GTF file size: ${GTF_SIZE}"
log "GTF entries with gene_type: ${GTF_GENES}"

#===============================================================================
# STEP 5: BUILD STAR INDEX
#===============================================================================

log "=============================================="
log "STEP 5: Building STAR genome index"
log "=============================================="

STAR_INDEX_DIR="${REF_DIR}/star"

# Check if index already exists
if [ -f "${STAR_INDEX_DIR}/SA" ] && [ -f "${STAR_INDEX_DIR}/Genome" ]; then
    log "STAR index already exists!"
    log "To rebuild, delete ${STAR_INDEX_DIR} and re-run this script"
    
    # Verify version compatibility
    if [ -f "${STAR_INDEX_DIR}/genomeParameters.txt" ]; then
        INDEX_VERSION=$(grep "versionGenome" "${STAR_INDEX_DIR}/genomeParameters.txt" | cut -f2)
        log "Existing index version: ${INDEX_VERSION}"
    fi
else
    log "Building STAR index with ${THREADS} threads..."
    log "This will take 30-60 minutes and requires ~32GB RAM..."
    log ""
    log "Parameters:"
    log "  Genome FASTA: ${GENOME_FA}"
    log "  GTF annotation: ${GTF_FILE}"
    log "  sjdbOverhang: ${SJDB_OVERHANG}"
    log "  Output directory: ${STAR_INDEX_DIR}"
    log ""
    
    # Record start time
    START_TIME=$(date +%s)
    
    STAR \
        --runMode genomeGenerate \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_INDEX_DIR}" \
        --genomeFastaFiles "${GENOME_FA}" \
        --sjdbGTFfile "${GTF_FILE}" \
        --sjdbOverhang ${SJDB_OVERHANG}
    
    if [ $? -ne 0 ]; then
        log "ERROR: STAR genome generation failed"
        exit 1
    fi
    
    # Record end time
    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    ELAPSED_MIN=$((ELAPSED / 60))
    
    log "STAR index generation completed in ${ELAPSED_MIN} minutes"
fi

#===============================================================================
# STEP 6: VERIFY INDEX
#===============================================================================

log "=============================================="
log "STEP 6: Verifying STAR index"
log "=============================================="

REQUIRED_FILES=("Genome" "SA" "SAindex" "chrLength.txt" "chrName.txt" "chrNameLength.txt" "chrStart.txt" "genomeParameters.txt")

ALL_PRESENT=true
for file in "${REQUIRED_FILES[@]}"; do
    if [ -f "${STAR_INDEX_DIR}/${file}" ]; then
        SIZE=$(du -h "${STAR_INDEX_DIR}/${file}" | cut -f1)
        log "  ✓ ${file} (${SIZE})"
    else
        log "  ✗ ${file} MISSING"
        ALL_PRESENT=false
    fi
done

if [ "${ALL_PRESENT}" = true ]; then
    log ""
    log "All required index files present!"
else
    log ""
    log "ERROR: Some index files are missing"
    exit 1
fi

# Show total index size
INDEX_SIZE=$(du -sh "${STAR_INDEX_DIR}" | cut -f1)
log "Total index size: ${INDEX_SIZE}"

#===============================================================================
# STEP 7: CREATE SUMMARY FILE
#===============================================================================

log "=============================================="
log "STEP 7: Creating reference summary"
log "=============================================="

SUMMARY_FILE="${REF_DIR}/reference_info.txt"

cat > "${SUMMARY_FILE}" << EOF
GRCh38 Reference Genome Summary
===============================
Created: $(date)

STAR Version: ${STAR_VERSION}
GENCODE Release: ${GENCODE_RELEASE}

Files:
------
Genome FASTA: ${GENOME_FA}
GTF Annotation: ${GTF_FILE}
STAR Index: ${STAR_INDEX_DIR}

Source URLs:
------------
Genome: ${GENOME_URL}
GTF: ${GTF_URL}

STAR Index Parameters:
----------------------
sjdbOverhang: ${SJDB_OVERHANG}
Threads used: ${THREADS}

Directory Structure:
--------------------
${REF_DIR}/
├── GRCh38.primary_assembly.genome.fa
├── gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf
├── star/
│   ├── Genome
│   ├── SA
│   ├── SAindex
│   └── ...
├── downloads/
│   ├── GRCh38.primary_assembly.genome.fa.gz
│   └── gencode.v${GENCODE_RELEASE}.primary_assembly.annotation.gtf.gz
├── star_version.txt
└── reference_info.txt

Usage in alignment scripts:
---------------------------
STAR_REF="${STAR_INDEX_DIR}"
EOF

log "Summary saved to: ${SUMMARY_FILE}"

#===============================================================================
# COMPLETE
#===============================================================================

log "=============================================="
log "REFERENCE GENOME BUILD COMPLETE"
log "=============================================="
log ""
log "STAR index location: ${STAR_INDEX_DIR}"
log ""
log "Update your alignment script to use:"
log "  STAR_REF=\"${STAR_INDEX_DIR}\""
log ""
log "You can now run the Slide-seq alignment pipeline!"
