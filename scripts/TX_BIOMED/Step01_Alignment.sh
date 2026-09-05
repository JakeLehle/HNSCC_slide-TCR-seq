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
# STEP 01: SLIDE-TCR-SEQ ALIGNMENT PIPELINE  (v2)
#
# Processes Sophia Liu's Slide-TCR-seq data from raw BCL files through to
# AnnData-compatible outputs with spatial coordinates.
#
# Pipeline:
#   0. BCL to FASTQ conversion (bcl2fastq)
#   1. Rewrite R1 to correct the split-barcode geometry        [NEW in v2]
#   2. Extract OBSERVED barcodes as whitelist                  [CHANGED in v2]
#   3. STAR alignment with STARsolo (per puck)
#   4. Validate barcode recovery, hard gate                    [NEW in v2]
#   5. Prepare DGE matrices + spatial coordinates (per puck)
#   6. Summary and QC
#
#-------------------------------------------------------------------------------
# CHANGELOG (v2) - two independent bugs, both confirmed empirically
#-------------------------------------------------------------------------------
# (a) R1 GEOMETRY. v1 ran --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15,
#     which assumes a contiguous 14bp barcode followed by the UMI. Slide-seq V2
#     R1 is 42bp with a SPLIT barcode:
#         1-8    bead barcode part 1
#         9-26   UP linker  TCTTCAGCGTTCCCGAGA
#         27-32  bead barcode part 2
#         33-41  UMI
#         42     spare cycle
#     v1 therefore built each CB as 8 real bases + TCTTCA (linker), and read the
#     UMI out of more linker. Confirmed: CR[9:14]=="TCTTCA" in 74.6% of reads,
#     UR=="GCGTTCCCG" in 66%, R1[9:26]==linker in 55.5%.
#     Result: 1.21% valid barcodes, 153 estimated cells (Puck_211214_29).
#     FIX: step1 rewrites R1 into a synthetic 23bp read (1-8 + 27-32 + 33-41),
#     which makes the solo parameters below correct AS WRITTEN. Do not "simplify"
#     by pointing STAR at the raw FASTQs again.
#
# (b) WHITELIST COLUMN. v1 used barcode_matching column 2 (CORRECTED) as the
#     STARsolo whitelist. A whitelist must contain barcodes AS THEY APPEAR IN
#     READS, which is column 1 (OBSERVED). Measured on 500k reads of puck 29:
#         corrected column (v1):  31.82% exact
#         observed  column (v2):  48.94% exact
#     Column 2 is also where all the N characters live (8,359 / 55,936 entries);
#     column 1 has zero. Switching columns removes the N problem entirely.
#     FIX: step2 extracts column 1, unique, no suffix strip (col 1 has none).
#
#-------------------------------------------------------------------------------
# CONSEQUENCE FOR DOWNSTREAM (read this before using Solo.out)
#-------------------------------------------------------------------------------
# Because the whitelist is OBSERVED barcodes, the CB tag written to the BAM is
# an observed VARIANT, not a bead. There are ~2.8 variants per bead (155,263
# observed -> 55,936 corrected for puck 29). Therefore:
#   - Solo.out will report ~155k "cells". This is expected, not a failure.
#   - Solo.out UMI dedup runs per variant, not per bead.
#   - The Solo.out matrix is NOT a drop-in replacement for Sophia's expression
#     matrix. Collapsing it by summing columns would double-count UMIs seen
#     under two variants of the same bead. A correct from-scratch DGE needs
#     UMI-aware collapse from the BAM (CB + UB). Separate script, not here.
#   - Step05a translates CB observed -> corrected -> +puck_id suffix before
#     SComatic. SComatic counts reads, not UMIs, so the variant-level CB is
#     harmless there.
# step5 writes obs2corr_{puck}.tsv, the map every downstream step should use.
#
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# Author: Jake Lehle, Texas Biomedical Research Institute
# Server: Zeus / Titan (Texas Biomed HPC)
#===============================================================================

set -o pipefail

source ~/anaconda3/bin/activate
conda activate slide-TCR-seq

#===============================================================================
# CONFIGURATION
#===============================================================================

PROJECT_ROOT="/master/jlehle/WORKING/slide-TCR-seq-working"

INPUT_DIR="${PROJECT_ROOT}/data/inputs"
FASTQ_BASE="${INPUT_DIR}/fastq"
STAR_REF="${INPUT_DIR}/ref/GRCh38/star"
GENOME_FA="${INPUT_DIR}/ref/GRCh38/GRCh38.primary_assembly.genome.fa"
GTF_FILE="${INPUT_DIR}/ref/GRCh38/gencode.v49.primary_assembly.annotation.gtf"
BCL_DIR="${FASTQ_BASE}/220116_NB501164_1345_AHLGH2BGXK"

OUTPUT_BASE="${PROJECT_ROOT}/data/outputs/01_alignment"
DEMUX_FASTQ_DIR="${OUTPUT_BASE}/demux_fastq"
FIXED_FASTQ_DIR="${OUTPUT_BASE}/demux_fastq_fixed"
WHITELIST_DIR="${OUTPUT_BASE}/whitelists"
ALIGNED_DIR="${OUTPUT_BASE}/aligned"
DGE_DIR="${OUTPUT_BASE}/dge"
LOG_DIR="${OUTPUT_BASE}/logs"
CHECKPOINT_DIR="${OUTPUT_BASE}/checkpoints"

# --- R1 geometry (Slide-seq V2 split barcode) ---
R1_RAW_LEN=42
BC1_START=1;  BC1_LEN=8       # barcode part 1
UP_START=9;   UP_LEN=18       # UP linker
BC2_START=27; BC2_LEN=6       # barcode part 2
RAW_UMI_START=33              # UMI in the RAW read
UP_LINKER="TCTTCAGCGTTCCCGAGA"
UP_MIN_FRAC=0.30              # abort if exact linker below this

# --- Synthetic R1 geometry (what STARsolo sees after step1) ---
CB_LEN=14                     # 8 + 6
UMI_START=15
UMI_LEN=9

# --- Gate: minimum acceptable valid-barcode rate from Solo.out ---
MIN_VALID_BC_RATE=0.30        # measured floor was 0.4894 exact, pre-1MM

# --- Barcode match type. Whitelist is dense (many entries mutually Hamming-1),
#     so multi-matches are expected. Most are variants of the SAME bead and are
#     harmless. If step4 reports high noTooManyWLmatches, set this to Exact. ---
CB_MATCH_TYPE="1MM_multi"

THREADS=$(nproc)
BAM_SORT_BINS=200
SORT_RAM=60000000000

SAMPLES=("Puck_211214_29" "Puck_211214_37" "Puck_211214_40")
SAMPLE_BARCODES=("AGATTTAA" "GGCGTCGA" "ATCACTCG")
SAMPLE_INPUT_DIRS=("2022-01-28_Puck_211214_29" "2022-01-28_Puck_211214_37" "2022-01-28_Puck_211214_40")

# --- Optional single-sample mode: Step01.sh --only Puck_211214_29 ---
if [ "${1:-}" == "--only" ] && [ -n "${2:-}" ]; then
    for i in "${!SAMPLES[@]}"; do
        if [ "${SAMPLES[$i]}" == "$2" ]; then
            SAMPLES=("${SAMPLES[$i]}")
            SAMPLE_BARCODES=("${SAMPLE_BARCODES[$i]}")
            SAMPLE_INPUT_DIRS=("${SAMPLE_INPUT_DIRS[$i]}")
            FOUND=1; break
        fi
    done
    [ -z "${FOUND:-}" ] && { echo "FATAL: unknown sample $2"; exit 1; }
fi

CURRENT_ULIMIT=$(ulimit -n)
TARGET_ULIMIT=65535
if [ "$CURRENT_ULIMIT" -lt "$TARGET_ULIMIT" ]; then
    ulimit -n $TARGET_ULIMIT 2>/dev/null
fi

command -v pigz >/dev/null && ZIP="pigz -p 8" || ZIP="gzip"

#===============================================================================
# SETUP
#===============================================================================

mkdir -p "${LOG_DIR}" "${DEMUX_FASTQ_DIR}" "${FIXED_FASTQ_DIR}" \
         "${WHITELIST_DIR}" "${CHECKPOINT_DIR}"
for sample in "${SAMPLES[@]}"; do
    mkdir -p "${ALIGNED_DIR}/${sample}" "${DGE_DIR}/${sample}"
done

LOG_FILE="${LOG_DIR}/Step01_pipeline_$(date +%Y%m%d_%H%M%S).log"

log() { local level=$1; shift; echo "[$(date '+%Y-%m-%d %H:%M:%S')] [${level}] $@" | tee -a "${LOG_FILE}"; }
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

set_checkpoint()   { touch "${CHECKPOINT_DIR}/${1}.done"; log_info "Checkpoint set: ${1}"; }
check_checkpoint() { [ -f "${CHECKPOINT_DIR}/${1}.done" ]; }
clear_all_checkpoints() { rm -f "${CHECKPOINT_DIR}"/*.done; log_info "All checkpoints cleared"; }

bc_match_file() {
    echo "${FASTQ_BASE}/${SAMPLE_INPUT_DIRS[$1]}/barcode_matching/${SAMPLES[$1]}_barcode_matching.txt.gz"
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

    local fastq_count=$(ls ${DEMUX_FASTQ_DIR}/Puck_*_R1_001.fastq.gz 2>/dev/null | wc -l)
    if [ ${fastq_count} -ge 3 ]; then
        log_info "FASTQ files already exist (${fastq_count} found), skipping bcl2fastq..."
        set_checkpoint "step0_bcl2fastq"
        return 0
    fi

    if [ ! -d "${BCL_DIR}" ]; then
        log_error "BCL directory not found: ${BCL_DIR}"
        return 1
    fi

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

    log_info "Running bcl2fastq (expected 30-60 min)..."
    bcl2fastq \
        --runfolder-dir "${BCL_DIR}" \
        --output-dir "${DEMUX_FASTQ_DIR}" \
        --sample-sheet "${OUTPUT_BASE}/SampleSheet.csv" \
        --no-lane-splitting \
        --processing-threads ${THREADS} \
        --barcode-mismatches 1 \
        2>&1 | tee -a "${LOG_FILE}"
    [ $? -ne 0 ] && { log_error "bcl2fastq failed"; return 1; }

    fastq_count=$(ls ${DEMUX_FASTQ_DIR}/Puck_*_R1_001.fastq.gz 2>/dev/null | wc -l)
    [ ${fastq_count} -lt 3 ] && { log_error "Expected 3 FASTQ sets, found ${fastq_count}"; return 1; }

    log_success "BCL to FASTQ complete (${fastq_count} sample sets)"
    set_checkpoint "step0_bcl2fastq"
    return 0
}

#===============================================================================
# STEP 1: FIX R1 SPLIT-BARCODE GEOMETRY  [NEW in v2]
#
# Rewrites 42bp R1 -> 23bp synthetic R1:  [1-8] + [27-32] + [33-41]
# Record count is preserved exactly so R1/R2 stay in lockstep (STAR pairs by
# record order). R2 and I1 are symlinked through so STAR sees one directory.
#===============================================================================

step1_fix_r1_geometry() {
    local SAMPLE=$1
    section_header "STEP 1: FIX R1 GEOMETRY - ${SAMPLE}"

    if check_checkpoint "step1_r1fix_${SAMPLE}"; then
        log_info "R1 geometry fix already completed for ${SAMPLE}, skipping..."
        return 0
    fi

    local IN_R1=$(ls ${DEMUX_FASTQ_DIR}/${SAMPLE}_S*_R1_001.fastq.gz 2>/dev/null | head -1)
    [ -z "${IN_R1}" ] && { log_error "Raw R1 not found for ${SAMPLE}"; return 1; }
    local OUT_R1="${FIXED_FASTQ_DIR}/$(basename ${IN_R1})"

    # --- Guard A: read length ---
    local LENS=$(zcat "${IN_R1}" | awk 'NR%4==2{print length($0)} NR>=400000{exit}' | sort -u | tr '\n' ',')
    log_info "  R1 lengths (first 100k reads): ${LENS}"
    local NSHORT=$(zcat "${IN_R1}" | awk -v m=$((RAW_UMI_START+UMI_LEN-1)) \
        'NR%4==2 && length($0)<m {c++} NR>=4000000{exit} END{print c+0}')
    log_info "  reads shorter than $((RAW_UMI_START+UMI_LEN-1))bp in first 1M: ${NSHORT}"

    # --- Guard B: UP linker at the expected position ---
    local NLINK=$(zcat "${IN_R1}" | awk -v L="${UP_LINKER}" -v s=${UP_START} -v n=${UP_LEN} \
        'NR%4==2{t++; if(substr($0,s,n)==L) c++} NR>=800000{exit} END{printf "%.4f", c/t}')
    log_info "  exact UP linker at ${UP_START}-$((UP_START+UP_LEN-1)): ${NLINK}"
    awk -v v="${NLINK}" -v m="${UP_MIN_FRAC}" \
        'BEGIN{ if (v+0 < m+0) exit 1 }' || {
        log_error "  FATAL: UP linker not found at expected position (${NLINK} < ${UP_MIN_FRAC})."
        log_error "  R1 geometry differs from Slide-seq V2. Do NOT proceed; re-derive offsets."
        return 1
    }

    # --- Rewrite (sequence and quality sliced identically) ---
    log_info "  Rewriting R1 -> ${OUT_R1}"
    zcat "${IN_R1}" | awk \
        -v b1s=${BC1_START} -v b1l=${BC1_LEN} \
        -v b2s=${BC2_START} -v b2l=${BC2_LEN} \
        -v us=${RAW_UMI_START} -v ul=${UMI_LEN} '
        NR%4==1 || NR%4==3 { print; next }
        {
            if (length($0) >= us+ul-1)
                print substr($0,b1s,b1l) substr($0,b2s,b2l) substr($0,us,ul)
            else { print ""; short++ }
        }
        END { if (short) print "WARN: " short " short reads emitted blank" > "/dev/stderr" }
    ' | ${ZIP} > "${OUT_R1}"

    # --- Verify: record count preserved, uniform 23bp ---
    local NIN=$(zcat "${IN_R1}"  | awk 'END{print NR/4}')
    local NOU=$(zcat "${OUT_R1}" | awk 'END{print NR/4}')
    local OLEN=$(zcat "${OUT_R1}" | awk 'NR%4==2{print length($0)}' | head -100000 | sort -u | tr '\n' ',')
    log_info "  reads in/out: ${NIN} / ${NOU}   synthetic lengths: ${OLEN}"
    [ "${NIN}" != "${NOU}" ] && { log_error "  FATAL: record count changed"; return 1; }

    # --- Symlink R2 / I1 through unchanged ---
    for M in R2 I1; do
        local SRC=$(ls ${DEMUX_FASTQ_DIR}/${SAMPLE}_S*_${M}_001.fastq.gz 2>/dev/null | head -1)
        [ -n "${SRC}" ] && ln -sf "${SRC}" "${FIXED_FASTQ_DIR}/$(basename ${SRC})"
    done

    log_success "R1 geometry fix complete for ${SAMPLE}"
    set_checkpoint "step1_r1fix_${SAMPLE}"
    return 0
}

#===============================================================================
# STEP 2: EXTRACT OBSERVED BARCODES AS WHITELIST  [CHANGED in v2]
#
# Column 1 (observed), not column 2 (corrected). See CHANGELOG (b).
# Also writes obs2corr_{puck}.tsv, the map for the Step05a retag.
#===============================================================================

step2_extract_observed_whitelist() {
    section_header "STEP 2: EXTRACT OBSERVED BARCODE WHITELIST"

    if check_checkpoint "step2_observed_whitelists"; then
        log_info "Observed whitelists already extracted, skipping..."
        return 0
    fi

    for i in "${!SAMPLES[@]}"; do
        local SAMPLE="${SAMPLES[$i]}"
        local BC_MATCH_FILE=$(bc_match_file $i)
        local WL_FILE="${WHITELIST_DIR}/observed_barcodes_${SAMPLE}.txt"
        local MAP_FILE="${WHITELIST_DIR}/obs2corr_${SAMPLE}.tsv"

        log_info "Processing ${SAMPLE}..."
        [ ! -f "${BC_MATCH_FILE}" ] && { log_error "Not found: ${BC_MATCH_FILE}"; return 1; }

        # Whitelist: column 1, unique. No suffix strip (col 1 carries none).
        zcat "${BC_MATCH_FILE}" | cut -f1 | sort -u > "${WL_FILE}"

        # Map: observed -> corrected (-1 stripped, matching h5ad obs_names).
        { printf "observed\tcorrected\n"
          zcat "${BC_MATCH_FILE}" | awk -F'\t' '{sub(/-1$/,"",$2); print $1"\t"$2}' | sort -u
        } > "${MAP_FILE}"

        local N_ROWS=$(zcat "${BC_MATCH_FILE}" | wc -l)
        local N_OBS=$(wc -l < "${WL_FILE}")
        local N_CORR=$(tail -n +2 "${MAP_FILE}" | cut -f2 | sort -u | wc -l)
        local BC_LENS=$(awk '{print length($0)}' "${WL_FILE}" | sort -u | tr '\n' ',')
        local N_WITH_N=$(grep -c N "${WL_FILE}" || true)
        local N_AMBIG=$(tail -n +2 "${MAP_FILE}" | cut -f1 | sort | uniq -d | wc -l)

        log_info "  rows in barcode_matching: ${N_ROWS}"
        log_info "  unique observed (whitelist): ${N_OBS}"
        log_info "  unique corrected (beads):    ${N_CORR}"
        log_info "  variants per bead:           $(awk -v a=${N_OBS} -v b=${N_CORR} 'BEGIN{printf "%.2f", a/b}')"
        log_info "  barcode lengths: ${BC_LENS}   containing N: ${N_WITH_N}"

        # --- Assertions ---
        [ "${BC_LENS}" != "${CB_LEN}," ] && {
            log_error "  FATAL: whitelist not uniformly ${CB_LEN}bp (got ${BC_LENS})"; return 1; }
        [ "${N_WITH_N}" -ne 0 ] && {
            log_error "  FATAL: ${N_WITH_N} observed barcodes contain N. Expected 0."
            log_error "  STARsolo rejects N in CB. Did the columns get swapped?"; return 1; }
        [ "${N_AMBIG}" -ne 0 ] && {
            log_error "  FATAL: ${N_AMBIG} observed barcodes map to >1 bead."
            log_error "  The retag translation is not a function. Needs a tie-break rule."; return 1; }

        log_info "  map written: ${MAP_FILE}"
    done

    log_success "Observed whitelist extraction complete"
    set_checkpoint "step2_observed_whitelists"
    return 0
}

#===============================================================================
# STEP 3: STAR ALIGNMENT (per puck)
#
# Solo parameters below are correct ONLY against the synthetic R1 from step1.
#===============================================================================

step3_star_alignment() {
    local SAMPLE=$1
    section_header "STEP 3: STAR ALIGNMENT - ${SAMPLE}"

    if check_checkpoint "step3_star_${SAMPLE}"; then
        log_info "STAR alignment already completed for ${SAMPLE}, skipping..."
        return 0
    fi

    local OUT_DIR="${ALIGNED_DIR}/${SAMPLE}"
    local WHITELIST="${WHITELIST_DIR}/observed_barcodes_${SAMPLE}.txt"
    local R1=$(ls ${FIXED_FASTQ_DIR}/${SAMPLE}_S*_R1_001.fastq.gz 2>/dev/null | head -1)
    local R2=$(ls ${FIXED_FASTQ_DIR}/${SAMPLE}_S*_R2_001.fastq.gz 2>/dev/null | head -1)

    [ -z "$R1" ] || [ -z "$R2" ] && { log_error "FASTQs not found in ${FIXED_FASTQ_DIR}"; return 1; }
    [ ! -f "${WHITELIST}" ] && { log_error "Whitelist not found: ${WHITELIST}"; return 1; }

    # Guard: R1 must be the SYNTHETIC read, not the raw 42bp one.
    local R1LEN=$(zcat "${R1}" | awk 'NR%4==2{print length($0)} NR>=40000{exit}' | sort -u | tr '\n' ',')
    [ "${R1LEN}" != "$((CB_LEN+UMI_LEN))," ] && {
        log_error "FATAL: R1 is ${R1LEN}, expected $((CB_LEN+UMI_LEN)),"
        log_error "Pointing STAR at raw R1 reproduces the v1 bug. Run step1 first."; return 1; }

    local NUM_BC=$(wc -l < "${WHITELIST}")
    log_info "  R1 (synthetic): ${R1}  [${R1LEN%,}bp]"
    log_info "  R2 (cDNA):      ${R2}"
    log_info "  Whitelist:      ${WHITELIST} (${NUM_BC} observed barcodes)"
    log_info "  CB ${CB_LEN}bp at 1-${CB_LEN}, UMI ${UMI_LEN}bp at ${UMI_START}-$((UMI_START+UMI_LEN-1))"
    log_info "  CB match type:  ${CB_MATCH_TYPE}"

    rm -rf "${OUT_DIR}"/_STARtmp "${OUT_DIR}"/Solo.out 2>/dev/null
    rm -f  "${OUT_DIR}"/_STAR* 2>/dev/null

    log_info "Starting STAR (threads ${THREADS}, expect 1-3h)..."
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
        --soloCBmatchWLtype ${CB_MATCH_TYPE} \
        --soloCellFilter None \
        --soloOutFileNames Solo.out/ genes.tsv barcodes.tsv matrix.mtx \
        --outReadsUnmapped Fastx \
        --outFilterScoreMinOverLread 0.3 \
        --outFilterMatchNminOverLread 0.3 \
        2>&1 | tee -a "${LOG_FILE}"

    local exit_code=$?
    [ ${exit_code} -ne 0 ] && { log_error "STAR failed (exit ${exit_code})"; return 1; }
    [ ! -f "${OUT_DIR}/Aligned.sortedByCoord.out.bam" ] && { log_error "BAM not created"; return 1; }

    log_info "Indexing BAM..."
    samtools index -@ ${THREADS} "${OUT_DIR}/Aligned.sortedByCoord.out.bam"

    log_info "Compressing unmapped reads..."
    [ -f "${OUT_DIR}/Unmapped.out.mate1" ] && ${ZIP} -f "${OUT_DIR}/Unmapped.out.mate1"
    [ -f "${OUT_DIR}/Unmapped.out.mate2" ] && ${ZIP} -f "${OUT_DIR}/Unmapped.out.mate2"

    rm -rf "${OUT_DIR}"/_STARtmp 2>/dev/null

    log_success "STAR alignment complete for ${SAMPLE}"
    set_checkpoint "step3_star_${SAMPLE}"
    return 0
}

#===============================================================================
# STEP 4: VALIDATE BARCODE RECOVERY  [NEW in v2]
#
# Hard gate. The v1 run produced a complete, well-formed, useless BAM: 50.9M
# mapped reads and 1.21% valid barcodes. Nothing downstream noticed. This stops
# that from happening silently again.
#===============================================================================

step4_validate_barcodes() {
    local SAMPLE=$1
    section_header "STEP 4: VALIDATE BARCODE RECOVERY - ${SAMPLE}"

    local SUMMARY="${ALIGNED_DIR}/${SAMPLE}/Solo.out/Gene/Summary.csv"
    local BCSTATS="${ALIGNED_DIR}/${SAMPLE}/Solo.out/Barcodes.stats"
    [ ! -f "${SUMMARY}" ] && { log_error "Summary.csv not found: ${SUMMARY}"; return 1; }

    local RATE=$(awk -F',' '/Reads With Valid Barcodes/{print $2}' "${SUMMARY}")
    local NREADS=$(awk -F',' '/^Number of Reads/{print $2}' "${SUMMARY}")
    local UNIQ=$(awk -F',' '/Reads Mapped to Genome: Unique,/{print $2}' "${SUMMARY}")

    log_info "  Number of Reads:           ${NREADS}"
    log_info "  Reads With Valid Barcodes: ${RATE}"
    log_info "  Reads Mapped Unique:       ${UNIQ}"
    log_info "  (v1 broken run for reference: 0.0120532)"

    if [ -f "${BCSTATS}" ]; then
        log_info "  Barcodes.stats:"
        sed 's/^/    /' "${BCSTATS}" | tee -a "${LOG_FILE}"
        local MULTWL=$(awk '/yesMultWLmatchWithMM/{print $2}' "${BCSTATS}")
        local TOOMANY=$(awk '/noTooManyWLmatches/{print $2}' "${BCSTATS}")
        log_info "  multi-WL-match with MM: ${MULTWL}   rejected too-many-matches: ${TOOMANY}"
        log_info "  (high values here => dense whitelist ambiguity; consider CB_MATCH_TYPE=Exact)"
    fi

    awk -v r="${RATE}" -v m="${MIN_VALID_BC_RATE}" 'BEGIN{ if (r+0 < m+0) exit 1 }' || {
        log_error "  FATAL: valid barcode rate ${RATE} below floor ${MIN_VALID_BC_RATE}"
        log_error "  Exact-match floor measured offline was 0.4894. Something regressed."
        log_error "  Check: is R1 the synthetic read? is the whitelist column 1?"
        return 1; }

    log_success "Barcode recovery acceptable for ${SAMPLE} (${RATE})"
    return 0
}

#===============================================================================
# STEP 5: PREPARE DGE + SPATIAL COORDINATES (per puck)
#
# NOTE: the Solo.out matrix is VARIANT-level, not bead-level. See the header
# CONSEQUENCE block. Raw matrix is used (soloCellFilter None); cell calling
# happens downstream against the annotated bead set, not here.
#===============================================================================

step5_prepare_dge() {
    local SAMPLE=$1
    local IDX=$2
    section_header "STEP 5: PREPARE DGE + SPATIAL - ${SAMPLE}"

    if check_checkpoint "step5_dge_${SAMPLE}"; then
        log_info "DGE preparation already completed for ${SAMPLE}, skipping..."
        return 0
    fi

    local SOLO_DIR="${ALIGNED_DIR}/${SAMPLE}/Solo.out/Gene"
    local OUT="${DGE_DIR}/${SAMPLE}"
    local BC_MATCH=$(bc_match_file ${IDX})
    local BC_XY="${FASTQ_BASE}/${SAMPLE_INPUT_DIRS[$IDX]}/barcode_matching/${SAMPLE}_barcode_xy.txt.gz"

    [ ! -d "${SOLO_DIR}/raw" ] && { log_error "Raw matrix not found: ${SOLO_DIR}/raw"; return 1; }

    log_info "Copying raw (variant-level) matrix to ${OUT}/"
    for file in matrix.mtx barcodes.tsv features.tsv genes.tsv; do
        [ -f "${SOLO_DIR}/raw/${file}" ] && { cp "${SOLO_DIR}/raw/${file}" "${OUT}/"; ${ZIP} -f "${OUT}/${file}"; }
    done
    [ -f "${OUT}/features.tsv.gz" ] && [ ! -f "${OUT}/genes.tsv.gz" ] && \
        cp "${OUT}/features.tsv.gz" "${OUT}/genes.tsv.gz"

    # Carry the map and the raw matching file alongside the matrix
    cp "${WHITELIST_DIR}/obs2corr_${SAMPLE}.tsv" "${OUT}/"
    [ -f "${BC_MATCH}" ] && cp "${BC_MATCH}" "${OUT}/"
    [ -f "${BC_XY}" ]    && cp "${BC_XY}" "${OUT}/"

    # Spatial coords keyed on CORRECTED barcode (matches h5ad obs_names)
    log_info "Creating spatial_coords.csv..."
    { echo "barcode,x,y"
      zcat "${BC_MATCH}" | awk -F'\t' '{sub(/-1$/,"",$2); print $2","$3","$4}' | sort -u
    } > "${OUT}/spatial_coords.csv"
    ${ZIP} -f "${OUT}/spatial_coords.csv"

    local n_var=$(zcat "${OUT}/barcodes.tsv.gz" | wc -l)
    log_info "  variant-level barcodes in matrix: ${n_var}"
    log_info "  NOTE: this is variants, not beads. Do not treat as a cell count."

    log_success "DGE + spatial preparation complete for ${SAMPLE}"
    set_checkpoint "step5_dge_${SAMPLE}"
    return 0
}

#===============================================================================
# STEP 6: SUMMARY
#===============================================================================

step6_summary() {
    section_header "STEP 01 PIPELINE SUMMARY"
    {
    echo ""
    echo "Geometry (Slide-seq V2 split barcode):"
    echo "  raw R1:       ${R1_RAW_LEN}bp"
    echo "  barcode:      ${BC1_START}-$((BC1_START+BC1_LEN-1)) + ${BC2_START}-$((BC2_START+BC2_LEN-1))  (${CB_LEN}bp)"
    echo "  UP linker:    ${UP_START}-$((UP_START+UP_LEN-1))  ${UP_LINKER}"
    echo "  UMI (raw):    ${RAW_UMI_START}-$((RAW_UMI_START+UMI_LEN-1))  (${UMI_LEN}bp)"
    echo "  synthetic R1: $((CB_LEN+UMI_LEN))bp, CB 1-${CB_LEN}, UMI ${UMI_START}-$((UMI_START+UMI_LEN-1))"
    echo "  whitelist:    barcode_matching column 1 (OBSERVED)"
    echo "  CB match:     ${CB_MATCH_TYPE}"
    echo "  STAR ref:     ${STAR_REF}"
    echo ""
    } | tee -a "${LOG_FILE}"

    for i in "${!SAMPLES[@]}"; do
        local SAMPLE="${SAMPLES[$i]}"
        local OUT="${ALIGNED_DIR}/${SAMPLE}"
        echo "=== ${SAMPLE} ===" | tee -a "${LOG_FILE}"
        local WL="${WHITELIST_DIR}/observed_barcodes_${SAMPLE}.txt"
        [ -f "${WL}" ] && echo "  Whitelist (observed): $(wc -l < ${WL})" | tee -a "${LOG_FILE}"
        [ -f "${WHITELIST_DIR}/obs2corr_${SAMPLE}.tsv" ] && \
            echo "  Beads (corrected):    $(tail -n +2 ${WHITELIST_DIR}/obs2corr_${SAMPLE}.tsv | cut -f2 | sort -u | wc -l)" | tee -a "${LOG_FILE}"
        [ -f "${OUT}/Log.final.out" ] && {
            echo "  STAR:" | tee -a "${LOG_FILE}"
            grep -E "Number of input reads|Uniquely mapped reads" "${OUT}/Log.final.out" \
                | sed 's/^/    /' | tee -a "${LOG_FILE}"; }
        [ -f "${OUT}/Solo.out/Gene/Summary.csv" ] && {
            echo "  STARsolo:" | tee -a "${LOG_FILE}"
            grep -E "Number of Reads|Valid Barcodes|Mapped to Genome: Unique,|Sequencing Saturation" \
                "${OUT}/Solo.out/Gene/Summary.csv" | sed 's/^/    /' | tee -a "${LOG_FILE}"; }
        echo "" | tee -a "${LOG_FILE}"
    done

    {
    echo "========================================"
    echo "Output Locations:"
    echo "  Raw FASTQs:      ${DEMUX_FASTQ_DIR}/"
    echo "  Fixed FASTQs:    ${FIXED_FASTQ_DIR}/"
    echo "  Whitelists+map:  ${WHITELIST_DIR}/"
    echo "  Aligned BAMs:    ${ALIGNED_DIR}/{puck}/"
    echo "  DGE matrices:    ${DGE_DIR}/{puck}/   (VARIANT-level)"
    echo "  Logs:            ${LOG_DIR}/"
    echo "========================================"
    } | tee -a "${LOG_FILE}"
}

#===============================================================================
# VALIDATION
#===============================================================================

validate_inputs() {
    section_header "VALIDATING INPUTS"
    local errors=0

    command -v STAR &> /dev/null && log_info "STAR: $(STAR --version 2>&1 | head -1)" \
        || { log_error "STAR not found"; ((errors++)); }
    command -v samtools &> /dev/null && log_info "samtools: $(samtools --version | head -1)" \
        || { log_error "samtools not found"; ((errors++)); }
    command -v bcl2fastq &> /dev/null && log_info "bcl2fastq: OK" \
        || log_warn "bcl2fastq not found (only needed if demux FASTQs missing)"
    command -v pigz &> /dev/null && log_info "pigz: OK (parallel gzip)" \
        || log_warn "pigz not found, falling back to gzip (slower)"

    [ -f "${STAR_REF}/SA" ] && log_info "STAR reference: OK" \
        || { log_error "STAR reference not found: ${STAR_REF}/SA"; ((errors++)); }
    [ -f "${GENOME_FA}" ] && log_info "Genome FASTA: OK" \
        || log_warn "Genome FASTA not found: ${GENOME_FA} (needed later by SComatic)"
    [ -d "${BCL_DIR}" ] && log_info "BCL directory: OK" \
        || log_warn "BCL directory not found (only needed if demux FASTQs missing)"

    for i in "${!SAMPLES[@]}"; do
        local f=$(bc_match_file $i)
        [ -f "${f}" ] && log_info "Barcode matching ${SAMPLES[$i]}: OK" \
            || { log_error "Not found: ${f}"; ((errors++)); }
    done

    [ ${errors} -gt 0 ] && { log_error "Validation failed with ${errors} error(s)"; return 1; }
    log_success "All inputs validated"
    return 0
}

#===============================================================================
# MAIN
#===============================================================================

main() {
    section_header "STEP 01: SLIDE-TCR-SEQ ALIGNMENT PIPELINE (v2)"
    log_info "Started at $(date)  on $(hostname)"
    log_info "Samples: ${SAMPLES[*]}"
    log_info "Threads: ${THREADS}"

    [ "${1:-}" == "--clean" ] && { log_warn "Cleaning all checkpoints..."; clear_all_checkpoints; }

    validate_inputs || exit 1

    local success=true

    step0_bcl_to_fastq || { log_error "Step 0 failed"; success=false; }

    if [ "$success" = true ]; then
        for SAMPLE in "${SAMPLES[@]}"; do
            step1_fix_r1_geometry "${SAMPLE}" || { log_error "Step 1 failed for ${SAMPLE}"; success=false; break; }
        done
    fi

    [ "$success" = true ] && { step2_extract_observed_whitelist || { log_error "Step 2 failed"; success=false; }; }

    if [ "$success" = true ]; then
        for i in "${!SAMPLES[@]}"; do
            SAMPLE="${SAMPLES[$i]}"
            log_info "Processing ${SAMPLE} ($((i+1))/${#SAMPLES[@]})"
            step3_star_alignment "${SAMPLE}"    || { log_error "Step 3 failed for ${SAMPLE}"; success=false; continue; }
            step4_validate_barcodes "${SAMPLE}" || { log_error "Step 4 GATE failed for ${SAMPLE}"; success=false; continue; }
            step5_prepare_dge "${SAMPLE}" "${i}" || { log_error "Step 5 failed for ${SAMPLE}"; success=false; continue; }
        done
    fi

    step6_summary
    section_header "STEP 01 COMPLETE"

    if [ "$success" = true ]; then
        log_success "Pipeline completed successfully at $(date)"
        log_info "Next: Step05a retag (CB observed -> corrected -> +puck_id) for SComatic"
        exit 0
    else
        log_error "Pipeline completed with errors at $(date)"
        exit 1
    fi
}

if [ "${1:-}" == "--help" ] || [ "${1:-}" == "-h" ]; then
    cat << EOF
Step 01 v2: Slide-TCR-seq Alignment Pipeline

Usage: sbatch $(basename $0)
       $(basename $0) --only Puck_211214_29   # single puck
       $(basename $0) --clean                 # clear checkpoints
       $(basename $0) --help

Pipeline:
    0. BCL to FASTQ (bcl2fastq)
    1. Fix R1 split-barcode geometry  -> demux_fastq_fixed/
    2. Extract OBSERVED barcode whitelist + obs2corr map
    3. STAR alignment with STARsolo (per puck)
    4. Validate barcode recovery (HARD GATE, floor ${MIN_VALID_BC_RATE})
    5. DGE matrix + spatial coordinates (per puck)
    6. Summary

v2 fixes two bugs in v1: R1 was parsed with contiguous-barcode geometry
(Slide-seq V2 is split around an 18bp UP linker), and the whitelist used
the corrected barcode column instead of the observed one. v1 yielded
1.21% valid barcodes; v2 should yield 49% or better.

The Solo.out matrix is VARIANT-level (~2.8 variants per bead), not
bead-level. Use obs2corr_{puck}.tsv to translate.
EOF
    exit 0
fi

main "$@"
