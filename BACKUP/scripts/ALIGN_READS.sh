#!/bin/bash

#SBATCH -J ALIGN_SLIDE-TCR-SEQ
#SBATCH -o /work/sdz852/WORKING/LOGS/ALIGN_SPATIAL.o.log
#SBATCH -e /work/sdz852/WORKING/LOGS/ALIGN_SPATIAL.e.log
#SBATCH --mail-user=jake.lehle@utsa.edu
#SBATCH --mail-type=ALL
#SBATCH -t 10-00:00:00
#SBATCH -p compute2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 80

#===============================================================================
# SLIDE-TCR-SEQ COMPLETE ALIGNMENT PIPELINE v5
#
# This pipeline processes Sophia Liu's Slide-TCR-seq data from raw BCL files
# through to AnnData-compatible outputs with spatial coordinates.
#
# Key features:
#   1. BCL to FASTQ conversion (bcl2fastq)
#   2. Uses CORRECTED barcodes from barcode_matching.txt.gz as whitelist
#   3. Uses 9bp UMI (confirmed from Sophia's BAM tags)
#   4. Expects ~50,000 cells per puck
#
# The barcode_matching.txt.gz file contains:
#   Col1: Observed barcode (from reads, may have errors)
#   Col2: Corrected barcode (matched to spatial position)
#   Col3: X coordinate
#   Col4: Y coordinate
#
# Author: Jake Lehle
# Project: HPV+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
#===============================================================================

set -o pipefail

# Startup
module load anaconda3
conda activate slide-TCR-seq

#===============================================================================
# SYSTEM RESOURCE SETTINGS
#===============================================================================

# Increase file handle limit for STAR
CURRENT_ULIMIT=$(ulimit -n)
TARGET_ULIMIT=65535

if [ "$CURRENT_ULIMIT" -lt "$TARGET_ULIMIT" ]; then
    ulimit -n $TARGET_ULIMIT 2>/dev/null
    NEW_ULIMIT=$(ulimit -n)
    if [ "$NEW_ULIMIT" -ge "$TARGET_ULIMIT" ]; then
        echo "Increased file handle limit: $CURRENT_ULIMIT -> $NEW_ULIMIT"
    else
        echo "Warning: Could not increase ulimit to $TARGET_ULIMIT (current: $NEW_ULIMIT)"
    fi
else
    echo "File handle limit sufficient: $CURRENT_ULIMIT"
fi

#===============================================================================
# CONFIGURATION
#===============================================================================

# Directories
WORKING_DIR="/work/sdz852/WORKING/slide-TCR-seq/working"
FASTQ_BASE="/work/sdz852/WORKING/slide-TCR-seq/fastq"
STAR_REF="/work/sdz852/WORKING/slide-TCR-seq/ref/GRCh38/star"
BCL_DIR="${FASTQ_BASE}/220116_NB501164_1345_AHLGH2BGXK"
WHITELIST_DIR="${WORKING_DIR}/whitelists"

# CRITICAL: Barcode parameters based on Sophia's processed data
CB_LEN=14       # Cell barcode length (confirmed)
UMI_START=15    # UMI starts at position 15
UMI_LEN=9       # UMI length is 9bp (confirmed from XM tag in Sophia's BAM)

# Processing parameters
CURRENT_ULIMIT=$(ulimit -n)
if [ "$CURRENT_ULIMIT" -ge 65535 ]; then
    THREADS=$(nproc)
    BAM_SORT_BINS=200
    echo "Using high-performance settings (THREADS=$THREADS, BAM_SORT_BINS=$BAM_SORT_BINS)"
else
    THREADS=16
    BAM_SORT_BINS=50
    echo "Using conservative settings (THREADS=$THREADS, BAM_SORT_BINS=$BAM_SORT_BINS)"
fi

# Sample info
SAMPLES=("Puck_211214_29" "Puck_211214_37" "Puck_211214_40")
SAMPLE_BARCODES=("AGATTTAA" "GGCGTCGA" "ATCACTCG")

#===============================================================================
# SETUP
#===============================================================================

mkdir -p "${WORKING_DIR}/logs"
mkdir -p "${WORKING_DIR}/fastq"
mkdir -p "${WORKING_DIR}/aligned"
mkdir -p "${WORKING_DIR}/dge"
mkdir -p "${WORKING_DIR}/checkpoints"
mkdir -p "${WHITELIST_DIR}"

for sample in "${SAMPLES[@]}"; do
    mkdir -p "${WORKING_DIR}/aligned/${sample}"
    mkdir -p "${WORKING_DIR}/dge/${sample}"
done

#===============================================================================
# LOGGING
#===============================================================================

LOG_FILE="${WORKING_DIR}/logs/pipeline_v5_$(date +%Y%m%d_%H%M%S).log"

log() {
    local level=$1
    shift
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level}] $@" | tee -a "${LOG_FILE}"
}

log_info() { log "INFO" "$@"; }
log_warn() { log "WARN" "$@"; }
log_error() { log "ERROR" "$@"; }
log_success() { log "SUCCESS" "$@"; }

section_header() {
    echo "" | tee -a "${LOG_FILE}"
    echo "========================================" | tee -a "${LOG_FILE}"
    echo "$1" | tee -a "${LOG_FILE}"
    echo "========================================" | tee -a "${LOG_FILE}"
}

#===============================================================================
# CHECKPOINT FUNCTIONS
#===============================================================================

CHECKPOINT_DIR="${WORKING_DIR}/checkpoints"

set_checkpoint() {
    touch "${CHECKPOINT_DIR}/${1}.done"
    log_info "Checkpoint set: ${1}"
}

check_checkpoint() {
    [ -f "${CHECKPOINT_DIR}/${1}.done" ]
}

clear_all_checkpoints() {
    rm -f "${CHECKPOINT_DIR}"/*.done
    log_info "All checkpoints cleared"
}

#===============================================================================
# STEP 0: BCL TO FASTQ CONVERSION
#===============================================================================

step0_bcl_to_fastq() {
    section_header "STEP 0: BCL TO FASTQ CONVERSION"
    
    if check_checkpoint "step0_bcl2fastq"; then
        log_info "Step 0 already completed, skipping..."
        return 0
    fi
    
    # Check if FASTQs already exist
    local fastq_count=$(ls ${WORKING_DIR}/fastq/Puck_*_R1_001.fastq.gz 2>/dev/null | wc -l)
    if [ ${fastq_count} -ge 3 ]; then
        log_info "FASTQ files already exist (${fastq_count} found), skipping bcl2fastq..."
        set_checkpoint "step0_bcl2fastq"
        return 0
    fi
    
    # Check BCL directory exists
    if [ ! -d "${BCL_DIR}" ]; then
        log_error "BCL directory not found: ${BCL_DIR}"
        return 1
    fi
    
    # Create sample sheet
    cat > "${WORKING_DIR}/SampleSheet.csv" << 'EOF'
[Header]
IEMFileVersion,4
Date,2022-01-16
Workflow,GenerateFASTQ
Application,NextSeq FASTQ Only

[Reads]
42
50

[Settings]
CreateFastqForIndexReads,1
MinimumTrimmedReadLength,0
MaskShortAdapterReads,0

[Data]
Sample_ID,Sample_Name,index
Puck_211214_29,Puck_211214_29,AGATTTAA
Puck_211214_37,Puck_211214_37,GGCGTCGA
Puck_211214_40,Puck_211214_40,ATCACTCG
EOF

    log_info "Created sample sheet"
    log_info "Running bcl2fastq..."
    log_info "This may take 30-60 minutes..."
    
    bcl2fastq \
        --runfolder-dir "${BCL_DIR}" \
        --output-dir "${WORKING_DIR}/fastq" \
        --sample-sheet "${WORKING_DIR}/SampleSheet.csv" \
        --no-lane-splitting \
        --processing-threads ${THREADS} \
        --barcode-mismatches 1 \
        2>&1 | tee -a "${LOG_FILE}"
    
    if [ $? -ne 0 ]; then
        log_error "bcl2fastq failed"
        return 1
    fi
    
    # Verify output
    local fastq_count=$(ls ${WORKING_DIR}/fastq/Puck_*_R1_001.fastq.gz 2>/dev/null | wc -l)
    if [ ${fastq_count} -lt 3 ]; then
        log_error "Expected 3 FASTQ sets, found ${fastq_count}"
        return 1
    fi
    
    log_success "BCL to FASTQ conversion complete"
    set_checkpoint "step0_bcl2fastq"
    return 0
}

#===============================================================================
# STEP 1: EXTRACT CORRECTED BARCODES FROM BARCODE_MATCHING FILES
#===============================================================================

step1_extract_corrected_barcodes() {
    section_header "STEP 1: EXTRACT CORRECTED BARCODES"
    
    if check_checkpoint "step1_corrected_whitelists"; then
        log_info "Corrected whitelists already extracted, skipping..."
        return 0
    fi
    
    for SAMPLE in "${SAMPLES[@]}"; do
        log_info "Processing ${SAMPLE}..."
        
        # Find the barcode_matching file
        local BC_MATCH_FILE="${FASTQ_BASE}/2022-01-28_${SAMPLE}/barcode_matching/${SAMPLE}_barcode_matching.txt.gz"
        local OUTPUT_FILE="${WHITELIST_DIR}/corrected_barcodes_${SAMPLE}.txt"
        
        if [ ! -f "${BC_MATCH_FILE}" ]; then
            log_error "Barcode matching file not found: ${BC_MATCH_FILE}"
            return 1
        fi
        
        # Extract unique corrected barcodes (column 2), remove the -1 suffix
        log_info "  Extracting corrected barcodes from ${BC_MATCH_FILE}..."
        
        zcat "${BC_MATCH_FILE}" | \
            cut -f2 | \
            sed 's/-1$//' | \
            sort -u > "${OUTPUT_FILE}"
        
        # Get statistics
        local TOTAL_OBSERVED=$(zcat "${BC_MATCH_FILE}" | wc -l)
        local UNIQUE_CORRECTED=$(wc -l < "${OUTPUT_FILE}")
        local BC_LENGTH=$(head -1 "${OUTPUT_FILE}" | tr -d '\n' | wc -c)
        
        log_info "  Total observed barcodes: ${TOTAL_OBSERVED}"
        log_info "  Unique corrected barcodes: ${UNIQUE_CORRECTED}"
        log_info "  Barcode length: ${BC_LENGTH}bp"
        
        # Verify barcode length
        if [ "${BC_LENGTH}" -ne "${CB_LEN}" ]; then
            log_warn "  Barcode length (${BC_LENGTH}) differs from expected (${CB_LEN})!"
        fi
        
        # Show first few barcodes
        log_info "  First 5 corrected barcodes:"
        head -5 "${OUTPUT_FILE}" | while read bc; do
            log_info "    ${bc}"
        done
    done
    
    log_success "Corrected barcode extraction complete"
    
    # Summary
    echo "" | tee -a "${LOG_FILE}"
    echo "Whitelist summary:" | tee -a "${LOG_FILE}"
    for f in "${WHITELIST_DIR}"/corrected_barcodes_*.txt; do
        if [ -f "$f" ]; then
            local count=$(wc -l < "$f")
            echo "  $(basename $f): ${count} barcodes" | tee -a "${LOG_FILE}"
        fi
    done
    
    set_checkpoint "step1_corrected_whitelists"
    return 0
}

#===============================================================================
# STEP 2: STAR ALIGNMENT WITH CORRECTED BARCODES
#===============================================================================

step2_star_alignment() {
    local SAMPLE=$1
    
    section_header "STEP 2: STAR ALIGNMENT - ${SAMPLE}"
    
    if check_checkpoint "step2_star_${SAMPLE}"; then
        log_info "STAR alignment already completed for ${SAMPLE}, skipping..."
        return 0
    fi
    
    local OUT_DIR="${WORKING_DIR}/aligned/${SAMPLE}"
    local WHITELIST="${WHITELIST_DIR}/corrected_barcodes_${SAMPLE}.txt"
    
    # Check if BAM already exists
    if [ -f "${OUT_DIR}/Aligned.sortedByCoord.out.bam" ] && [ -f "${OUT_DIR}/Aligned.sortedByCoord.out.bam.bai" ]; then
        log_info "BAM file already exists, skipping alignment..."
        set_checkpoint "step2_star_${SAMPLE}"
        return 0
    fi
    
    # Find FASTQ files
    local R1=$(ls ${WORKING_DIR}/fastq/${SAMPLE}_S*_R1_001.fastq.gz 2>/dev/null | head -1)
    local R2=$(ls ${WORKING_DIR}/fastq/${SAMPLE}_S*_R2_001.fastq.gz 2>/dev/null | head -1)
    
    if [ -z "$R1" ] || [ -z "$R2" ]; then
        log_error "FASTQ files not found for ${SAMPLE}"
        log_error "Expected in: ${WORKING_DIR}/fastq/"
        return 1
    fi
    
    if [ ! -f "${WHITELIST}" ]; then
        log_error "Whitelist not found: ${WHITELIST}"
        log_error "Run step 1 first"
        return 1
    fi
    
    local NUM_BC=$(wc -l < "${WHITELIST}")
    
    log_info "Input files:"
    log_info "  R1: ${R1}"
    log_info "  R2: ${R2}"
    log_info "  Whitelist: ${WHITELIST} (${NUM_BC} barcodes)"
    log_info ""
    log_info "Barcode structure:"
    log_info "  Cell Barcode: ${CB_LEN}bp (positions 1-${CB_LEN})"
    log_info "  UMI: ${UMI_LEN}bp (positions ${UMI_START}-$((UMI_START + UMI_LEN - 1)))"
    
    # Clean up any previous failed runs
    rm -rf "${OUT_DIR}"/_STARtmp 2>/dev/null
    rm -f "${OUT_DIR}"/_STAR* 2>/dev/null
    
    log_info ""
    log_info "Starting STAR alignment..."
    log_info "Expected runtime: 1-3 hours per sample"
    
    STAR \
        --runThreadN ${THREADS} \
        --genomeDir "${STAR_REF}" \
        --readFilesIn "${R2}" "${R1}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${OUT_DIR}/" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --outBAMsortingBinsN ${BAM_SORT_BINS} \
        --limitBAMsortRAM 60000000000 \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "${WHITELIST}" \
        --soloCBstart 1 \
        --soloCBlen ${CB_LEN} \
        --soloUMIstart ${UMI_START} \
        --soloUMIlen ${UMI_LEN} \
        --soloBarcodeReadLength 0 \
        --soloFeatures Gene GeneFull \
        --soloUMIdedup 1MM_All \
        --soloCBmatchWLtype 1MM_multi \
        --soloCellFilter EmptyDrops_CR \
        --soloOutFileNames Solo.out/ genes.tsv barcodes.tsv matrix.mtx \
        --outReadsUnmapped Fastx \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        2>&1 | tee -a "${LOG_FILE}"
    
    local exit_code=$?
    
    if [ ${exit_code} -ne 0 ]; then
        log_error "STAR failed with exit code ${exit_code}"
        return 1
    fi
    
    if [ ! -f "${OUT_DIR}/Aligned.sortedByCoord.out.bam" ]; then
        log_error "BAM file not created"
        return 1
    fi
    
    log_info "Indexing BAM file..."
    samtools index "${OUT_DIR}/Aligned.sortedByCoord.out.bam"
    
    log_info "Compressing unmapped reads..."
    [ -f "${OUT_DIR}/Unmapped.out.mate1" ] && gzip -f "${OUT_DIR}/Unmapped.out.mate1"
    [ -f "${OUT_DIR}/Unmapped.out.mate2" ] && gzip -f "${OUT_DIR}/Unmapped.out.mate2"
    
    log_success "STAR alignment complete for ${SAMPLE}"
    set_checkpoint "step2_star_${SAMPLE}"
    return 0
}

#===============================================================================
# STEP 3: PREPARE ANNDATA INPUTS
#===============================================================================

step3_prepare_anndata() {
    local SAMPLE=$1
    
    section_header "STEP 3: PREPARE ANNDATA INPUTS - ${SAMPLE}"
    
    if check_checkpoint "step3_anndata_${SAMPLE}"; then
        log_info "AnnData preparation already completed for ${SAMPLE}, skipping..."
        return 0
    fi
    
    local OUT_DIR="${WORKING_DIR}/aligned/${SAMPLE}"
    local DGE_DIR="${WORKING_DIR}/dge/${SAMPLE}"
    local SOLO_DIR="${OUT_DIR}/Solo.out/Gene"
    
    if [ ! -d "${SOLO_DIR}" ]; then
        log_error "STARsolo output not found: ${SOLO_DIR}"
        return 1
    fi
    
    # Determine which matrix to use
    local MATRIX_DIR=""
    if [ -d "${SOLO_DIR}/filtered" ] && [ -f "${SOLO_DIR}/filtered/matrix.mtx" ]; then
        MATRIX_DIR="${SOLO_DIR}/filtered"
        log_info "Using filtered matrix"
    elif [ -d "${SOLO_DIR}/raw" ] && [ -f "${SOLO_DIR}/raw/matrix.mtx" ]; then
        MATRIX_DIR="${SOLO_DIR}/raw"
        log_info "Using raw matrix (filtered not available)"
    else
        log_error "No matrix found in ${SOLO_DIR}"
        return 1
    fi
    
    log_info "Copying matrix files..."
    
    for file in matrix.mtx barcodes.tsv features.tsv genes.tsv; do
        if [ -f "${MATRIX_DIR}/${file}" ]; then
            cp "${MATRIX_DIR}/${file}" "${DGE_DIR}/"
            gzip -f "${DGE_DIR}/${file}"
            log_info "  ${file} -> ${file}.gz"
        fi
    done
    
    # Ensure genes.tsv exists (some STAR versions output features.tsv)
    if [ -f "${DGE_DIR}/features.tsv.gz" ] && [ ! -f "${DGE_DIR}/genes.tsv.gz" ]; then
        cp "${DGE_DIR}/features.tsv.gz" "${DGE_DIR}/genes.tsv.gz"
        log_info "  Copied features.tsv.gz -> genes.tsv.gz"
    fi
    
    # Copy spatial coordinates from barcode_matching file
    local SOPHIA_DIR="${FASTQ_BASE}/2022-01-28_${SAMPLE}/barcode_matching"
    local BC_MATCH="${SOPHIA_DIR}/${SAMPLE}_barcode_matching.txt.gz"
    local BC_XY="${SOPHIA_DIR}/${SAMPLE}_barcode_xy.txt.gz"
    
    if [ -f "${BC_MATCH}" ]; then
        cp "${BC_MATCH}" "${DGE_DIR}/"
        log_info "  Copied: $(basename ${BC_MATCH})"
    fi
    
    if [ -f "${BC_XY}" ]; then
        cp "${BC_XY}" "${DGE_DIR}/"
        log_info "  Copied: $(basename ${BC_XY})"
    fi
    
    # Create a spatial coordinates file in standard format
    # Format: barcode, x, y
    log_info "Creating spatial_coords.csv..."
    echo "barcode,x,y" > "${DGE_DIR}/spatial_coords.csv"
    zcat "${BC_MATCH}" | \
        cut -f2,3,4 | \
        sort -u | \
        sed 's/\t/,/g' >> "${DGE_DIR}/spatial_coords.csv"
    gzip -f "${DGE_DIR}/spatial_coords.csv"
    log_info "  Created spatial_coords.csv.gz"
    
    log_success "AnnData inputs prepared for ${SAMPLE}"
    set_checkpoint "step3_anndata_${SAMPLE}"
    return 0
}

#===============================================================================
# STEP 4: SUMMARY AND QC
#===============================================================================

step4_summary() {
    section_header "PIPELINE SUMMARY"
    
    echo "" | tee -a "${LOG_FILE}"
    echo "Pipeline Parameters:" | tee -a "${LOG_FILE}"
    echo "  Cell Barcode: ${CB_LEN}bp (positions 1-${CB_LEN})" | tee -a "${LOG_FILE}"
    echo "  UMI: ${UMI_LEN}bp (positions ${UMI_START}-$((UMI_START + UMI_LEN - 1)))" | tee -a "${LOG_FILE}"
    echo "  Whitelist source: barcode_matching.txt.gz (corrected barcodes)" | tee -a "${LOG_FILE}"
    echo "" | tee -a "${LOG_FILE}"
    
    for SAMPLE in "${SAMPLES[@]}"; do
        echo "=== ${SAMPLE} ===" | tee -a "${LOG_FILE}"
        
        local OUT_DIR="${WORKING_DIR}/aligned/${SAMPLE}"
        local DGE_DIR="${WORKING_DIR}/dge/${SAMPLE}"
        
        # Whitelist stats
        local WL="${WHITELIST_DIR}/corrected_barcodes_${SAMPLE}.txt"
        if [ -f "${WL}" ]; then
            echo "Whitelist barcodes: $(wc -l < ${WL})" | tee -a "${LOG_FILE}"
        fi
        
        # STAR alignment stats
        if [ -f "${OUT_DIR}/Log.final.out" ]; then
            echo "" | tee -a "${LOG_FILE}"
            echo "STAR Alignment:" | tee -a "${LOG_FILE}"
            grep -E "Number of input reads|Uniquely mapped reads|% of reads mapped" \
                "${OUT_DIR}/Log.final.out" | tee -a "${LOG_FILE}"
        fi
        
        # STARsolo stats
        if [ -f "${OUT_DIR}/Solo.out/Gene/Summary.csv" ]; then
            echo "" | tee -a "${LOG_FILE}"
            echo "STARsolo Summary:" | tee -a "${LOG_FILE}"
            cat "${OUT_DIR}/Solo.out/Gene/Summary.csv" | tee -a "${LOG_FILE}"
        fi
        
        # Cell count from filtered matrix
        if [ -f "${DGE_DIR}/barcodes.tsv.gz" ]; then
            local cell_count=$(zcat "${DGE_DIR}/barcodes.tsv.gz" | wc -l)
            echo "" | tee -a "${LOG_FILE}"
            echo "Cells in filtered matrix: ${cell_count}" | tee -a "${LOG_FILE}"
        fi
        
        echo "" | tee -a "${LOG_FILE}"
    done
    
    echo "========================================" | tee -a "${LOG_FILE}"
    echo "Output Locations:" | tee -a "${LOG_FILE}"
    echo "  FASTQ files: ${WORKING_DIR}/fastq/" | tee -a "${LOG_FILE}"
    echo "  Corrected whitelists: ${WHITELIST_DIR}/" | tee -a "${LOG_FILE}"
    echo "  Aligned BAMs: ${WORKING_DIR}/aligned/{sample}/" | tee -a "${LOG_FILE}"
    echo "  DGE matrices: ${WORKING_DIR}/dge/{sample}/" | tee -a "${LOG_FILE}"
    echo "========================================" | tee -a "${LOG_FILE}"
}

#===============================================================================
# VALIDATION
#===============================================================================

validate_inputs() {
    section_header "VALIDATING INPUTS"
    
    local errors=0
    
    # Check STAR
    if command -v STAR &> /dev/null; then
        log_info "STAR: $(STAR --version 2>&1 | head -1)"
    else
        log_error "STAR not found"
        ((errors++))
    fi
    
    # Check bcl2fastq
    if command -v bcl2fastq &> /dev/null; then
        log_info "bcl2fastq: OK"
    else
        log_error "bcl2fastq not found"
        ((errors++))
    fi
    
    # Check STAR reference
    if [ -f "${STAR_REF}/SA" ]; then
        log_info "STAR reference: OK"
    else
        log_error "STAR reference not found: ${STAR_REF}"
        ((errors++))
    fi
    
    # Check BCL directory exists
    if [ -d "${BCL_DIR}" ]; then
        log_info "BCL directory: OK"
    else
        log_error "BCL directory not found: ${BCL_DIR}"
        ((errors++))
    fi
    
    # Check barcode_matching files exist
    for SAMPLE in "${SAMPLES[@]}"; do
        local bc_file="${FASTQ_BASE}/2022-01-28_${SAMPLE}/barcode_matching/${SAMPLE}_barcode_matching.txt.gz"
        if [ -f "${bc_file}" ]; then
            log_info "Barcode matching ${SAMPLE}: OK"
        else
            log_error "Barcode matching file not found: ${bc_file}"
            ((errors++))
        fi
    done
    
    # Check samtools
    if command -v samtools &> /dev/null; then
        log_info "samtools: OK"
    else
        log_error "samtools not found"
        ((errors++))
    fi
    
    if [ ${errors} -gt 0 ]; then
        log_error "Validation failed with ${errors} error(s)"
        return 1
    fi
    
    log_success "All inputs validated"
    return 0
}

#===============================================================================
# MAIN
#===============================================================================

main() {
    section_header "SLIDE-TCR-SEQ COMPLETE PIPELINE v5"
    
    log_info "Started at $(date)"
    log_info "Working directory: ${WORKING_DIR}"
    log_info ""
    log_info "Key parameters:"
    log_info "  Cell Barcode: ${CB_LEN}bp"
    log_info "  UMI: ${UMI_LEN}bp"
    log_info "  Whitelist: Corrected barcodes from barcode_matching.txt.gz"
    
    # Handle --clean flag
    if [ "$1" == "--clean" ]; then
        log_warn "Cleaning all checkpoints..."
        clear_all_checkpoints
    fi
    
    # Validate
    if ! validate_inputs; then
        exit 1
    fi
    
    local success=true
    
    # Step 0: BCL to FASTQ
    if ! step0_bcl_to_fastq; then
        log_error "Step 0 (bcl2fastq) failed"
        success=false
    fi
    
    # Step 1: Extract corrected barcodes
    if [ "$success" = true ] && ! step1_extract_corrected_barcodes; then
        log_error "Step 1 (extract barcodes) failed"
        success=false
    fi
    
    # Steps 2-3: Process each sample
    if [ "$success" = true ]; then
        for i in "${!SAMPLES[@]}"; do
            SAMPLE="${SAMPLES[$i]}"
            log_info "Processing ${SAMPLE} ($((i+1))/${#SAMPLES[@]})"
            
            if ! step2_star_alignment "${SAMPLE}"; then
                log_error "Step 2 (STAR) failed for ${SAMPLE}"
                success=false
                continue
            fi
            
            if ! step3_prepare_anndata "${SAMPLE}"; then
                log_error "Step 3 (AnnData prep) failed for ${SAMPLE}"
                success=false
                continue
            fi
        done
    fi
    
    # Summary
    step4_summary
    
    section_header "PIPELINE COMPLETE"
    
    if [ "$success" = true ]; then
        log_success "Pipeline completed successfully at $(date)"
        exit 0
    else
        log_error "Pipeline completed with errors at $(date)"
        exit 1
    fi
}

# Help
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    cat << EOF
Slide-TCR-seq Complete Alignment Pipeline v5

This pipeline processes Sophia Liu's Slide-TCR-seq data from raw BCL files
through to AnnData-compatible outputs with spatial coordinates.

Key features:
  - BCL to FASTQ conversion (bcl2fastq)
  - Uses CORRECTED barcodes from barcode_matching.txt.gz
  - Uses 9bp UMI (confirmed from Sophia's processed data)
  - Expects ~50,000 cells per puck

Usage: $(basename $0) [OPTIONS]

Options:
    --clean     Clear all checkpoints and start fresh
    --help      Show this help message

Pipeline steps:
    0. BCL to FASTQ conversion (bcl2fastq)
    1. Extract corrected barcodes from barcode_matching.txt.gz
    2. STAR alignment with corrected whitelist (per sample)
    3. Prepare AnnData-compatible outputs (per sample)
    4. Summary and QC

Expected outputs:
    - ~50,000 cells per puck (vs ~150 with incorrect whitelist)
    - Reads with valid barcodes: ~50-70% (vs ~1% previously)

Output locations:
    - FASTQ files: ${WORKING_DIR}/fastq/
    - Whitelists: ${WORKING_DIR}/whitelists/
    - Aligned BAMs: ${WORKING_DIR}/aligned/{sample}/
    - DGE matrices: ${WORKING_DIR}/dge/{sample}/

EOF
    exit 0
fi

main "$@"
