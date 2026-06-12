#!/bin/bash

#SBATCH -J SPATIAL_01_ALIGN
#SBATCH -o /master/jlehle/WORKING/LOGS/Step01_Align.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step01_Align.e.%j.log
#SBATCH --mail-user=jlehle@txbiomed.org
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00
#SBATCH -p normal
#SBATCH --mem=900G
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 80

#===============================================================================
# STEP 01: SLIDE-TCR-SEQ ALIGNMENT PIPELINE
#
# Processes Sophia Liu's Slide-TCR-seq data from raw BCL files through to
# AnnData-compatible outputs with spatial coordinates.
#
# Pipeline:
#   0. BCL to FASTQ conversion (bcl2fastq)
#   1. Extract CORRECTED barcodes from barcode_matching.txt.gz as whitelist
#   2. STAR alignment with STARsolo (per puck)
#   3. Prepare DGE matrices + spatial coordinates (per puck)
#   4. Summary and QC
#
# Key parameters:
#   - Cell barcode: 14bp (positions 1-14)
#   - UMI: 9bp (positions 15-23, confirmed from Sophia's BAM XM tag)
#   - Expected: ~50,000 cells per puck
#
# The barcode_matching.txt.gz file contains:
#   Col1: Observed barcode (from reads, may have errors)
#   Col2: Corrected barcode (matched to spatial position)
#   Col3: X coordinate (microns)
#   Col4: Y coordinate (microns)
#
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# Author: Jake Lehle, Texas Biomedical Research Institute
# Server: Zeus / Titan (Texas Biomed HPC)
#===============================================================================

set -o pipefail

# Startup
source ~/anaconda3/bin/activate
conda activate slide-TCR-seq

#===============================================================================
# CONFIGURATION
# All paths relative to the unified project structure.
# Python scripts import equivalent paths from spatial_config.py.
#===============================================================================

# --- Project root ---
PROJECT_ROOT="/master/jlehle/WORKING/slide-TCR-seq-working"

# --- Input directories (read-only) ---
INPUT_DIR="${PROJECT_ROOT}/data/inputs"
FASTQ_BASE="${INPUT_DIR}/fastq"
STAR_REF="${INPUT_DIR}/ref/GRCh38/star"
GENOME_FA="${INPUT_DIR}/ref/GRCh38/GRCh38.primary_assembly.genome.fa"
GTF_FILE="${INPUT_DIR}/ref/GRCh38/gencode.v49.primary_assembly.annotation.gtf"
BCL_DIR="${FASTQ_BASE}/220116_NB501164_1345_AHLGH2BGXK"

# --- Output directories ---
OUTPUT_BASE="${PROJECT_ROOT}/data/outputs/01_alignment"
DEMUX_FASTQ_DIR="${OUTPUT_BASE}/demux_fastq"
WHITELIST_DIR="${OUTPUT_BASE}/whitelists"
ALIGNED_DIR="${OUTPUT_BASE}/aligned"
DGE_DIR="${OUTPUT_BASE}/dge"
LOG_DIR="${OUTPUT_BASE}/logs"
CHECKPOINT_DIR="${OUTPUT_BASE}/checkpoints"

# --- Barcode parameters (confirmed from Sophia's processed data) ---
CB_LEN=14       # Cell barcode length
UMI_START=15    # UMI starts at position 15
UMI_LEN=9       # UMI length (confirmed from XM tag in Sophia's BAM)

# --- Processing parameters ---
THREADS=$(nproc)
BAM_SORT_BINS=200
SORT_RAM=60000000000   # 60GB for BAM sorting

# --- Sample definitions ---
SAMPLES=("Puck_211214_29" "Puck_211214_37" "Puck_211214_40")
SAMPLE_BARCODES=("AGATTTAA" "GGCGTCGA" "ATCACTCG")
SAMPLE_INPUT_DIRS=("2022-01-28_Puck_211214_29" "2022-01-28_Puck_211214_37" "2022-01-28_Puck_211214_40")

# --- System tuning ---
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
# SETUP: Create output directory tree
#===============================================================================

mkdir -p "${LOG_DIR}"
mkdir -p "${DEMUX_FASTQ_DIR}"
mkdir -p "${WHITELIST_DIR}"
mkdir -p "${CHECKPOINT_DIR}"

for sample in "${SAMPLES[@]}"; do
    mkdir -p "${ALIGNED_DIR}/${sample}"
    mkdir -p "${DGE_DIR}/${sample}"
done

#===============================================================================
# LOGGING
#===============================================================================

LOG_FILE="${LOG_DIR}/Step01_pipeline_$(date +%Y%m%d_%H%M%S).log"

log() {
    local level=$1
    shift
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level}] $@" | tee -a "${LOG_FILE}"
}

log_info()    { log "INFO" "$@"; }
log_warn()    { log "WARN" "$@"; }
log_error()   { log "ERROR" "$@"; }
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
    local fastq_count=$(ls ${DEMUX_FASTQ_DIR}/Puck_*_R1_001.fastq.gz 2>/dev/null | wc -l)
    if [ ${fastq_count} -ge 3 ]; then
        log_info "FASTQ files already exist (${fastq_count} found), skipping bcl2fastq..."
        set_checkpoint "step0_bcl2fastq"
        return 0
    fi

    # Check BCL directory
    if [ ! -d "${BCL_DIR}" ]; then
        log_error "BCL directory not found: ${BCL_DIR}"
        return 1
    fi

    # Create sample sheet
    cat > "${OUTPUT_BASE}/SampleSheet.csv" << 'EOF'
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

    log_info "Created sample sheet: ${OUTPUT_BASE}/SampleSheet.csv"
    log_info "Running bcl2fastq..."
    log_info "Expected runtime: 30-60 minutes"

    bcl2fastq \
        --runfolder-dir "${BCL_DIR}" \
        --output-dir "${DEMUX_FASTQ_DIR}" \
        --sample-sheet "${OUTPUT_BASE}/SampleSheet.csv" \
        --no-lane-splitting \
        --processing-threads ${THREADS} \
        --barcode-mismatches 1 \
        2>&1 | tee -a "${LOG_FILE}"

    if [ $? -ne 0 ]; then
        log_error "bcl2fastq failed"
        return 1
    fi

    # Verify output
    local fastq_count=$(ls ${DEMUX_FASTQ_DIR}/Puck_*_R1_001.fastq.gz 2>/dev/null | wc -l)
    if [ ${fastq_count} -lt 3 ]; then
        log_error "Expected 3 FASTQ sets, found ${fastq_count}"
        return 1
    fi

    log_success "BCL to FASTQ conversion complete (${fastq_count} sample sets)"
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

    for i in "${!SAMPLES[@]}"; do
        local SAMPLE="${SAMPLES[$i]}"
        local INPUT_SUBDIR="${SAMPLE_INPUT_DIRS[$i]}"
        log_info "Processing ${SAMPLE}..."

        # Find the barcode_matching file
        local BC_MATCH_FILE="${FASTQ_BASE}/${INPUT_SUBDIR}/barcode_matching/${SAMPLE}_barcode_matching.txt.gz"
        local OUTPUT_FILE="${WHITELIST_DIR}/corrected_barcodes_${SAMPLE}.txt"

        if [ ! -f "${BC_MATCH_FILE}" ]; then
            log_error "Barcode matching file not found: ${BC_MATCH_FILE}"
            return 1
        fi

        # Extract unique corrected barcodes (column 2), remove the -1 suffix
        log_info "  Extracting corrected barcodes from: $(basename ${BC_MATCH_FILE})"

        zcat "${BC_MATCH_FILE}" | \
            cut -f2 | \
            sed 's/-1$//' | \
            sort -u > "${OUTPUT_FILE}"

        # Statistics
        local TOTAL_OBSERVED=$(zcat "${BC_MATCH_FILE}" | wc -l)
        local UNIQUE_CORRECTED=$(wc -l < "${OUTPUT_FILE}")
        local BC_LENGTH=$(head -1 "${OUTPUT_FILE}" | tr -d '\n' | wc -c)

        log_info "  Total observed barcodes: ${TOTAL_OBSERVED}"
        log_info "  Unique corrected barcodes: ${UNIQUE_CORRECTED}"
        log_info "  Barcode length: ${BC_LENGTH}bp"

        # Verify barcode length matches expected
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
# STEP 2: STAR ALIGNMENT WITH CORRECTED BARCODES (per puck)
#===============================================================================

step2_star_alignment() {
    local SAMPLE=$1

    section_header "STEP 2: STAR ALIGNMENT - ${SAMPLE}"

    if check_checkpoint "step2_star_${SAMPLE}"; then
        log_info "STAR alignment already completed for ${SAMPLE}, skipping..."
        return 0
    fi

    local OUT_DIR="${ALIGNED_DIR}/${SAMPLE}"
    local WHITELIST="${WHITELIST_DIR}/corrected_barcodes_${SAMPLE}.txt"

    # Check if BAM already exists
    if [ -f "${OUT_DIR}/Aligned.sortedByCoord.out.bam" ] && \
       [ -f "${OUT_DIR}/Aligned.sortedByCoord.out.bam.bai" ]; then
        log_info "BAM file already exists, skipping alignment..."
        set_checkpoint "step2_star_${SAMPLE}"
        return 0
    fi

    # Find demultiplexed FASTQ files
    local R1=$(ls ${DEMUX_FASTQ_DIR}/${SAMPLE}_S*_R1_001.fastq.gz 2>/dev/null | head -1)
    local R2=$(ls ${DEMUX_FASTQ_DIR}/${SAMPLE}_S*_R2_001.fastq.gz 2>/dev/null | head -1)

    if [ -z "$R1" ] || [ -z "$R2" ]; then
        log_error "FASTQ files not found for ${SAMPLE}"
        log_error "Expected in: ${DEMUX_FASTQ_DIR}/"
        return 1
    fi

    if [ ! -f "${WHITELIST}" ]; then
        log_error "Whitelist not found: ${WHITELIST}"
        log_error "Run step 1 first"
        return 1
    fi

    local NUM_BC=$(wc -l < "${WHITELIST}")

    log_info "Input files:"
    log_info "  R1 (barcodes): ${R1}"
    log_info "  R2 (cDNA):     ${R2}"
    log_info "  Whitelist:     ${WHITELIST} (${NUM_BC} barcodes)"
    log_info ""
    log_info "Barcode structure:"
    log_info "  Cell Barcode: ${CB_LEN}bp (positions 1-${CB_LEN})"
    log_info "  UMI: ${UMI_LEN}bp (positions ${UMI_START}-$((UMI_START + UMI_LEN - 1)))"

    # Clean up any previous failed runs
    rm -rf "${OUT_DIR}"/_STARtmp 2>/dev/null
    rm -f "${OUT_DIR}"/_STAR* 2>/dev/null

    log_info ""
    log_info "Starting STAR alignment..."
    log_info "Threads: ${THREADS}"
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
        --limitBAMsortRAM ${SORT_RAM} \
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
# STEP 3: PREPARE DGE MATRICES + SPATIAL COORDINATES (per puck)
#===============================================================================

step3_prepare_dge() {
    local SAMPLE=$1
    local IDX=$2

    section_header "STEP 3: PREPARE DGE + SPATIAL - ${SAMPLE}"

    if check_checkpoint "step3_dge_${SAMPLE}"; then
        log_info "DGE preparation already completed for ${SAMPLE}, skipping..."
        return 0
    fi

    local STAR_OUT="${ALIGNED_DIR}/${SAMPLE}"
    local OUT="${DGE_DIR}/${SAMPLE}"
    local SOLO_DIR="${STAR_OUT}/Solo.out/Gene"
    local INPUT_SUBDIR="${SAMPLE_INPUT_DIRS[$IDX]}"

    if [ ! -d "${SOLO_DIR}" ]; then
        log_error "STARsolo output not found: ${SOLO_DIR}"
        return 1
    fi

    # Determine which matrix to use (prefer filtered)
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

    log_info "Copying matrix files to: ${OUT}/"

    for file in matrix.mtx barcodes.tsv features.tsv genes.tsv; do
        if [ -f "${MATRIX_DIR}/${file}" ]; then
            cp "${MATRIX_DIR}/${file}" "${OUT}/"
            gzip -f "${OUT}/${file}"
            log_info "  ${file} -> ${file}.gz"
        fi
    done

    # Ensure genes.tsv exists (some STAR versions output features.tsv)
    if [ -f "${OUT}/features.tsv.gz" ] && [ ! -f "${OUT}/genes.tsv.gz" ]; then
        cp "${OUT}/features.tsv.gz" "${OUT}/genes.tsv.gz"
        log_info "  Copied features.tsv.gz -> genes.tsv.gz"
    fi

    # Copy spatial coordinate files from Sophia's input data
    local SOPHIA_BC_DIR="${FASTQ_BASE}/${INPUT_SUBDIR}/barcode_matching"
    local BC_MATCH="${SOPHIA_BC_DIR}/${SAMPLE}_barcode_matching.txt.gz"
    local BC_XY="${SOPHIA_BC_DIR}/${SAMPLE}_barcode_xy.txt.gz"

    if [ -f "${BC_MATCH}" ]; then
        cp "${BC_MATCH}" "${OUT}/"
        log_info "  Copied: $(basename ${BC_MATCH})"
    fi

    if [ -f "${BC_XY}" ]; then
        cp "${BC_XY}" "${OUT}/"
        log_info "  Copied: $(basename ${BC_XY})"
    fi

    # Create standardized spatial coordinates CSV
    # Format: barcode,x,y (corrected barcodes with -1 suffix stripped)
    log_info "Creating spatial_coords.csv..."
    echo "barcode,x,y" > "${OUT}/spatial_coords.csv"
    zcat "${BC_MATCH}" | \
        awk -F'\t' '{gsub(/-1$/, "", $2); print $2","$3","$4}' | \
        sort -u >> "${OUT}/spatial_coords.csv"
    gzip -f "${OUT}/spatial_coords.csv"
    log_info "  Created spatial_coords.csv.gz"

    # Report cell count
    if [ -f "${OUT}/barcodes.tsv.gz" ]; then
        local cell_count=$(zcat "${OUT}/barcodes.tsv.gz" | wc -l)
        log_info "  Cells in DGE matrix: ${cell_count}"
    fi

    log_success "DGE + spatial preparation complete for ${SAMPLE}"
    set_checkpoint "step3_dge_${SAMPLE}"
    return 0
}

#===============================================================================
# STEP 4: SUMMARY AND QC
#===============================================================================

step4_summary() {
    section_header "STEP 01 PIPELINE SUMMARY"

    echo "" | tee -a "${LOG_FILE}"
    echo "Pipeline Parameters:" | tee -a "${LOG_FILE}"
    echo "  Cell Barcode: ${CB_LEN}bp (positions 1-${CB_LEN})" | tee -a "${LOG_FILE}"
    echo "  UMI: ${UMI_LEN}bp (positions ${UMI_START}-$((UMI_START + UMI_LEN - 1)))" | tee -a "${LOG_FILE}"
    echo "  Whitelist source: barcode_matching.txt.gz (corrected barcodes)" | tee -a "${LOG_FILE}"
    echo "  STAR reference: ${STAR_REF}" | tee -a "${LOG_FILE}"
    echo "" | tee -a "${LOG_FILE}"

    for SAMPLE in "${SAMPLES[@]}"; do
        echo "=== ${SAMPLE} ===" | tee -a "${LOG_FILE}"

        local OUT="${ALIGNED_DIR}/${SAMPLE}"
        local DGE="${DGE_DIR}/${SAMPLE}"

        # Whitelist stats
        local WL="${WHITELIST_DIR}/corrected_barcodes_${SAMPLE}.txt"
        if [ -f "${WL}" ]; then
            echo "  Whitelist barcodes: $(wc -l < ${WL})" | tee -a "${LOG_FILE}"
        fi

        # STAR alignment stats
        if [ -f "${OUT}/Log.final.out" ]; then
            echo "" | tee -a "${LOG_FILE}"
            echo "  STAR Alignment:" | tee -a "${LOG_FILE}"
            grep -E "Number of input reads|Uniquely mapped reads|% of reads mapped" \
                "${OUT}/Log.final.out" | sed 's/^/    /' | tee -a "${LOG_FILE}"
        fi

        # STARsolo stats
        if [ -f "${OUT}/Solo.out/Gene/Summary.csv" ]; then
            echo "" | tee -a "${LOG_FILE}"
            echo "  STARsolo Summary:" | tee -a "${LOG_FILE}"
            cat "${OUT}/Solo.out/Gene/Summary.csv" | sed 's/^/    /' | tee -a "${LOG_FILE}"
        fi

        # Cell count from filtered matrix
        if [ -f "${DGE}/barcodes.tsv.gz" ]; then
            local cell_count=$(zcat "${DGE}/barcodes.tsv.gz" | wc -l)
            echo "" | tee -a "${LOG_FILE}"
            echo "  Cells in DGE matrix: ${cell_count}" | tee -a "${LOG_FILE}"
        fi

        echo "" | tee -a "${LOG_FILE}"
    done

    echo "========================================" | tee -a "${LOG_FILE}"
    echo "Output Locations:" | tee -a "${LOG_FILE}"
    echo "  Demux FASTQs:    ${DEMUX_FASTQ_DIR}/" | tee -a "${LOG_FILE}"
    echo "  Whitelists:      ${WHITELIST_DIR}/" | tee -a "${LOG_FILE}"
    echo "  Aligned BAMs:    ${ALIGNED_DIR}/{puck}/" | tee -a "${LOG_FILE}"
    echo "  DGE matrices:    ${DGE_DIR}/{puck}/" | tee -a "${LOG_FILE}"
    echo "  Logs:            ${LOG_DIR}/" | tee -a "${LOG_FILE}"
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
        log_error "STAR not found in PATH"
        ((errors++))
    fi

    # Check bcl2fastq
    if command -v bcl2fastq &> /dev/null; then
        log_info "bcl2fastq: OK"
    else
        log_warn "bcl2fastq not found (only needed if demux FASTQs missing)"
    fi

    # Check samtools
    if command -v samtools &> /dev/null; then
        log_info "samtools: $(samtools --version | head -1)"
    else
        log_error "samtools not found in PATH"
        ((errors++))
    fi

    # Check STAR reference
    if [ -f "${STAR_REF}/SA" ]; then
        log_info "STAR reference: OK (${STAR_REF})"
    else
        log_error "STAR reference not found: ${STAR_REF}/SA"
        ((errors++))
    fi

    # Check BCL directory
    if [ -d "${BCL_DIR}" ]; then
        log_info "BCL directory: OK"
    else
        log_warn "BCL directory not found: ${BCL_DIR}"
        log_warn "  (only needed if demux FASTQs missing)"
    fi

    # Check barcode_matching files
    for i in "${!SAMPLES[@]}"; do
        local SAMPLE="${SAMPLES[$i]}"
        local INPUT_SUBDIR="${SAMPLE_INPUT_DIRS[$i]}"
        local bc_file="${FASTQ_BASE}/${INPUT_SUBDIR}/barcode_matching/${SAMPLE}_barcode_matching.txt.gz"
        if [ -f "${bc_file}" ]; then
            log_info "Barcode matching ${SAMPLE}: OK"
        else
            log_error "Barcode matching not found: ${bc_file}"
            ((errors++))
        fi
    done

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
    section_header "STEP 01: SLIDE-TCR-SEQ ALIGNMENT PIPELINE"

    log_info "Started at $(date)"
    log_info "Server: $(hostname)"
    log_info "Project root: ${PROJECT_ROOT}"
    log_info ""
    log_info "Key parameters:"
    log_info "  Cell Barcode: ${CB_LEN}bp"
    log_info "  UMI: ${UMI_LEN}bp"
    log_info "  Threads: ${THREADS}"
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

    # Steps 2-3: Process each puck
    if [ "$success" = true ]; then
        for i in "${!SAMPLES[@]}"; do
            SAMPLE="${SAMPLES[$i]}"
            log_info "Processing ${SAMPLE} ($((i+1))/${#SAMPLES[@]})"

            if ! step2_star_alignment "${SAMPLE}"; then
                log_error "Step 2 (STAR) failed for ${SAMPLE}"
                success=false
                continue
            fi

            if ! step3_prepare_dge "${SAMPLE}" "${i}"; then
                log_error "Step 3 (DGE prep) failed for ${SAMPLE}"
                success=false
                continue
            fi
        done
    fi

    # Summary
    step4_summary

    section_header "STEP 01 COMPLETE"

    if [ "$success" = true ]; then
        log_success "Pipeline completed successfully at $(date)"
        log_info ""
        log_info "Next: Run Step02_AnnData_Prep.py to load DGE matrices into AnnData"
        exit 0
    else
        log_error "Pipeline completed with errors at $(date)"
        exit 1
    fi
}

# Help
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    cat << EOF
Step 01: Slide-TCR-seq Alignment Pipeline

Processes Sophia Liu's Slide-TCR-seq data from raw BCL files through to
AnnData-compatible outputs with spatial coordinates.

Usage: sbatch $(basename $0)
       $(basename $0) --clean    # Clear checkpoints, start fresh
       $(basename $0) --help     # Show this message

Pipeline:
    0. BCL to FASTQ (bcl2fastq)
    1. Extract corrected barcodes (whitelist)
    2. STAR alignment with STARsolo (per puck)
    3. DGE matrix + spatial coordinates (per puck)
    4. Summary

Outputs (all under data/outputs/01_alignment/):
    demux_fastq/                Demultiplexed FASTQ files
    whitelists/                 Corrected barcode lists
    aligned/{puck}/             STAR BAM + Solo.out
    dge/{puck}/                 DGE matrix + spatial coords
    logs/                       Pipeline logs

Expected:
    ~50,000 cells per puck
    ~50-70% reads with valid barcodes

EOF
    exit 0
fi

main "$@"
