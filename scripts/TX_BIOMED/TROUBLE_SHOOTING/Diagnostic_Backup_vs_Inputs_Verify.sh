#!/usr/bin/env bash
# =============================================================================
# Diagnostic_Backup_vs_Inputs_Verify.sh   v2
# =============================================================================
# READ-ONLY. Creates no links, copies nothing, deletes nothing.
#
# v2 change: the lane-BAM integrity sweep is promoted to BEAT 1, because it is
# now the question that decides the shape of the project. All 18 unmapped.bam
# files failed `samtools quickcheck` in the original delivery, which killed the
# HPV16 arm. An intact BGZF EOF marker was observed on one of them, so the
# failure was never truncation. This script settles whether they are corrupt,
# whether the re-delivery fixes them, or whether they were always fine and
# quickcheck was rejecting a headerless pre-alignment BAM.
#
#   BEAT 0  disk space and tree footprints
#   BEAT 1  lane BAM integrity, all 36 files, BACKUP vs inputs
#   BEAT 2  what unmapped.bam actually is  (interpretation A vs B)
#   BEAT 3  BACKUP vs inputs, path + size, per puck
#   BEAT 4  are our alignment/ subdirectories populated
#   BEAT 5  what are the six TCR files
#   BEAT 6  head-hash spot check on large equal-size files
#
#   bash Diagnostic_Backup_vs_Inputs_Verify.sh              # fast, no read counts
#   bash Diagnostic_Backup_vs_Inputs_Verify.sh --counts     # add samtools view -c
#                                                           # (slow, GB-scale BAMs)
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
BK="$PROOT/BACKUP"
IN="$PROOT/data/inputs/fastq"
OUT="$PROOT/data/outputs/00_inventory"
REPORT="$OUT/backup_vs_inputs_report.txt"
BAMTSV="$OUT/lane_bam_integrity.tsv"

PUCKS=(29 37 40)
DO_COUNTS=0
for a in "$@"; do [[ "$a" == "--counts" ]] && DO_COUNTS=1; done

mkdir -p "$OUT"
: > "$REPORT"

log() { echo "$@" | tee -a "$REPORT"; }
banner() {
    log ""
    log "========================================================================"
    log "$1"
    log "========================================================================"
}

# system samtools is broken on this cluster, use the conda build
if ! command -v samtools >/dev/null 2>&1 || ! samtools --version >/dev/null 2>&1; then
    source ~/anaconda3/bin/activate slide-TCR-seq 2>/dev/null || true
fi
SAMTOOLS=$(command -v samtools || echo "")

log "Diagnostic_Backup_vs_Inputs_Verify.sh  v2"
log "run:      $(date)"
log "host:     $(hostname)"
log "samtools: ${SAMTOOLS:-NOT FOUND}"
[[ -n "$SAMTOOLS" ]] && log "version:  $("$SAMTOOLS" --version 2>/dev/null | head -1)"
log "counts:   $DO_COUNTS"

# -----------------------------------------------------------------------------
banner "BEAT 0  disk space and tree footprints"
# -----------------------------------------------------------------------------
df -h "$PROOT" | tee -a "$REPORT"
log ""
log "BACKUP total:  $(du -sh "$BK" 2>/dev/null | cut -f1)"
log "inputs total:  $(du -sh "$PROOT/data/inputs" 2>/dev/null | cut -f1)"
for P in "${PUCKS[@]}"; do
    D="2022-01-28_Puck_211214_${P}"
    log "  puck $P  BACKUP $(du -sh "$BK/$D" 2>/dev/null | cut -f1)   inputs $(du -sh "$IN/$D" 2>/dev/null | cut -f1)"
done

# -----------------------------------------------------------------------------
banner "BEAT 1  lane BAM integrity, all 36 files, BACKUP vs inputs"
# -----------------------------------------------------------------------------
# The canonical BGZF end-of-file block. Present = the writer closed the file
# cleanly, so a quickcheck failure is a header problem, not truncation.
BGZF_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"

bam_probe() {
    local B="$1" TREE="$2"
    local rel="${B#$BK/}"; rel="${rel#$IN/}"

    if [[ ! -f "$B" ]]; then
        printf "%s\t%s\tABSENT\t-\t-\t-\t-\t-\t-\n" "$TREE" "$rel" >> "$BAMTSV"
        return
    fi

    local size magic eof qc nsq nrec firstlen
    size=$(stat -c%s "$B")
    magic=$(head -c 4 "$B" | xxd -p)
    if [[ "$(tail -c 28 "$B" | xxd -p | tr -d '\n')" == "$BGZF_EOF" ]]; then eof=EOF_OK; else eof=EOF_MISSING; fi

    if "$SAMTOOLS" quickcheck "$B" 2>/dev/null; then qc=PASS; else qc=FAIL; fi

    # the discriminating test: can the header and a record be read out of a file
    # that quickcheck rejects?
    nsq=$("$SAMTOOLS" view -H "$B" 2>/dev/null | grep -c '^@SQ')
    [[ -z "$nsq" ]] && nsq="ERR"
    firstlen=$("$SAMTOOLS" view "$B" 2>/dev/null | head -1 | awk '{print length($10)}')
    [[ -z "$firstlen" ]] && firstlen="NO_RECORD"

    if (( DO_COUNTS )); then
        nrec=$("$SAMTOOLS" view -c "$B" 2>/dev/null || echo ERR)
        [[ -z "$nrec" ]] && nrec=ERR
    else
        nrec="skipped"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$TREE" "$rel" "$size" "$magic" "$eof" "$qc" "$nsq" "$firstlen" "$nrec" >> "$BAMTSV"
}

if [[ -z "$SAMTOOLS" ]]; then
    log "  samtools not available, skipping BEAT 1 and BEAT 2"
else
    printf "tree\tpath\tsize\tmagic\tbgzf_eof\tquickcheck\tn_SQ\tread1_len\tn_records\n" > "$BAMTSV"

    for P in "${PUCKS[@]}"; do
        D="2022-01-28_Puck_211214_${P}"
        for FC in H52J2DMXY HLGH2BGXK; do
            LANES=(L001 L002); [[ "$FC" == HLGH2BGXK ]] && LANES=(L001 L002 L003 L004)
            for L in "${LANES[@]}"; do
                for KIND in final unmapped; do
                    bam_probe "$BK/$D/$FC/$L/Puck_211214_${P}.${KIND}.bam" BACKUP
                    bam_probe "$IN/$D/$FC/$L/Puck_211214_${P}.${KIND}.bam" inputs
                done
            done
        done
    done

    log "  full table: $BAMTSV"
    log ""
    log "  --- quickcheck result, unmapped.bam, by tree ---"
    awk -F'\t' 'NR>1 && $2 ~ /unmapped/ {c[$1" "$6]++} END{for(x in c) print "    "x"  n="c[x]}' "$BAMTSV" | sort | tee -a "$REPORT"
    log "  --- quickcheck result, final.bam, by tree ---"
    awk -F'\t' 'NR>1 && $2 ~ /final/ {c[$1" "$6]++} END{for(x in c) print "    "x"  n="c[x]}' "$BAMTSV" | sort | tee -a "$REPORT"
    log "  --- BGZF EOF marker ---"
    awk -F'\t' 'NR>1 {c[$1" "$5]++} END{for(x in c) print "    "x"  n="c[x]}' "$BAMTSV" | sort | tee -a "$REPORT"
    log "  --- @SQ line counts (0 means a pre-alignment unaligned BAM) ---"
    awk -F'\t' 'NR>1 {k=($2 ~ /unmapped/ ? "unmapped":"final"); c[$1" "k" nSQ="$7]++} END{for(x in c) print "    "x"  n="c[x]}' "$BAMTSV" | sort | tee -a "$REPORT"
    log ""
    log "  --- files where quickcheck FAILED but a record still read out ---"
    awk -F'\t' 'NR>1 && $6=="FAIL" && $8!="NO_RECORD" {print "    "$1"  "$2"  read1_len="$8}' "$BAMTSV" | tee -a "$REPORT"
    log "      (any row above means the file is readable and quickcheck is the wrong test)"
    log ""
    log "  --- size agreement, BACKUP vs inputs ---"
    awk -F'\t' 'NR>1 {s[$2"|"$1]=$3; seen[$2]=1}
        END{ok=0;bad=0;new=0
            for (p in seen){b=s[p"|BACKUP"];i=s[p"|inputs"]
              if(i=="ABSENT"||i==""){new++; print "    NEW_IN_BACKUP  "p}
              else if(b==i){ok++} else {bad++; print "    SIZE_DIFFER    "p"  backup="b"  inputs="i}}
            print "    equal="ok"  differ="bad"  new="new}' "$BAMTSV" | tee -a "$REPORT"

    log ""
    log "  DECISION TABLE"
    log "    inputs FAIL, BACKUP PASS                  -> transfer damage, re-delivery fixes it"
    log "    both FAIL, same size, records readable    -> never corrupt, quickcheck semantics"
    log "    both FAIL, different size                 -> bad at source, no delivery fixes it"
    log "    both PASS                                 -> original finding was a tooling artifact"
fi

# -----------------------------------------------------------------------------
banner "BEAT 2  what unmapped.bam actually is  (interpretation A vs B)"
# -----------------------------------------------------------------------------
# A = pre-alignment tagged unaligned BAM. Every read in the library, XC and XM
#     already set, zero @SQ lines. FASTQ recoverable, HPV16 reads are in here.
# B = post-alignment unaligned-only. Small fraction of final.bam, @SQ present.

if [[ -n "$SAMTOOLS" ]]; then
    D="2022-01-28_Puck_211214_29"
    for B in "$BK/$D/H52J2DMXY/L001/Puck_211214_29.final.bam" \
             "$BK/$D/H52J2DMXY/L001/Puck_211214_29.unmapped.bam"; do
        [[ -f "$B" ]] || { log "  missing: $B"; continue; }
        log ""
        log "--- $(basename "$B")   $(du -h "$B" | cut -f1)"
        log "  header, non-@SQ lines:"
        "$SAMTOOLS" view -H "$B" 2>/dev/null | grep -v '^@SQ' | cut -c1-220 | sed 's/^/      /' | tee -a "$REPORT"
        log "  @SQ count: $("$SAMTOOLS" view -H "$B" 2>/dev/null | grep -c '^@SQ')"
        log "  first record (300 char):"
        "$SAMTOOLS" view "$B" 2>/dev/null | head -1 | cut -c1-300 | sed 's/^/      /' | tee -a "$REPORT"
        log "  tags on record 1: $("$SAMTOOLS" view "$B" 2>/dev/null | head -1 | tr '\t' '\n' | grep -E '^[A-Za-z][A-Za-z0-9]:[AifZHB]:' | cut -d: -f1 | tr '\n' ' ')"
        log "  FLAG, first 5:       $("$SAMTOOLS" view "$B" 2>/dev/null | head -5 | awk '{printf "%s ", $2}')"
        log "  SEQ length, first 5: $("$SAMTOOLS" view "$B" 2>/dev/null | head -5 | awk '{printf "%s ", length($10)}')"
        (( DO_COUNTS )) && log "  read count: $("$SAMTOOLS" view -c "$B" 2>/dev/null)"
    done
    log ""
    log "  INTERPRETATION GUIDE"
    log "    @SQ == 0, FLAG == 4 throughout, XC and XM present, count >= final"
    log "       -> A. The whole library. samtools fastq recovers R2 with quals."
    log "             R1 must be synthesized from XC + XM, so it is derived, not raw."
    log "             HPV16 reads are in here and were never given a chance to align."
    log "    @SQ > 0 and count is a small fraction of final"
    log "       -> B. Post-alignment unaligned only. HPV16 reservoir, not a FASTQ source."
fi

# -----------------------------------------------------------------------------
banner "BEAT 3  path + size diff, per puck"
# -----------------------------------------------------------------------------
TOTAL_ONLY_BK=0; TOTAL_ONLY_IN=0; TOTAL_DIFF=0; TOTAL_SAME=0

for P in "${PUCKS[@]}"; do
    D="2022-01-28_Puck_211214_${P}"
    B_LIST="$OUT/backup_puck${P}.tsv"; I_LIST="$OUT/inputs_puck${P}.tsv"

    find "$BK/$D" -type f -printf "%P\t%s\n" 2>/dev/null | sort > "$B_LIST"
    find "$IN/$D" -type f -printf "%P\t%s\n" 2>/dev/null | sort > "$I_LIST"
    cut -f1 "$B_LIST" | sort > "$OUT/.bp_$P"; cut -f1 "$I_LIST" | sort > "$OUT/.ip_$P"

    comm -23 "$OUT/.bp_$P" "$OUT/.ip_$P" > "$OUT/only_in_backup_puck${P}.txt"
    comm -13 "$OUT/.bp_$P" "$OUT/.ip_$P" > "$OUT/only_in_inputs_puck${P}.txt"
    comm -12 "$OUT/.bp_$P" "$OUT/.ip_$P" > "$OUT/.shared_$P"

    : > "$OUT/size_mismatch_puck${P}.txt"
    SAME=0; DIFFN=0
    while IFS= read -r f; do
        bs=$(awk -F'\t' -v k="$f" '$1==k{print $2; exit}' "$B_LIST")
        is=$(awk -F'\t' -v k="$f" '$1==k{print $2; exit}' "$I_LIST")
        if [[ "$bs" == "$is" ]]; then SAME=$((SAME+1))
        else DIFFN=$((DIFFN+1))
             printf "%s\tbackup=%s\tinputs=%s\tdelta=%s\n" "$f" "$bs" "$is" "$((bs - is))" \
                >> "$OUT/size_mismatch_puck${P}.txt"; fi
    done < "$OUT/.shared_$P"

    OB=$(wc -l < "$OUT/only_in_backup_puck${P}.txt"); OI=$(wc -l < "$OUT/only_in_inputs_puck${P}.txt")
    log ""
    log "puck $P   BACKUP $(wc -l < "$B_LIST") files, inputs $(wc -l < "$I_LIST") files"
    log "  only in BACKUP (new):  $OB"
    log "  only in inputs:        $OI"
    log "  shared, equal size:    $SAME"
    log "  shared, SIZE MISMATCH: $DIFFN"
    (( DIFFN > 0 )) && { log "  *** first 10 mismatches ***"; head -10 "$OUT/size_mismatch_puck${P}.txt" | sed 's/^/      /' | tee -a "$REPORT"; }
    (( OB > 0 ))    && { log "  *** first 15 new ***";        head -15 "$OUT/only_in_backup_puck${P}.txt" | sed 's/^/      /' | tee -a "$REPORT"; }

    TOTAL_ONLY_BK=$((TOTAL_ONLY_BK+OB)); TOTAL_ONLY_IN=$((TOTAL_ONLY_IN+OI))
    TOTAL_DIFF=$((TOTAL_DIFF+DIFFN)); TOTAL_SAME=$((TOTAL_SAME+SAME))
    rm -f "$OUT/.bp_$P" "$OUT/.ip_$P" "$OUT/.shared_$P"
done

log ""
log "ACROSS ALL THREE PUCKS   new=$TOTAL_ONLY_BK  inputs_only=$TOTAL_ONLY_IN  equal=$TOTAL_SAME  mismatch=$TOTAL_DIFF"

# -----------------------------------------------------------------------------
banner "BEAT 4  are our alignment/ subdirectories populated"
# -----------------------------------------------------------------------------
for P in "${PUCKS[@]}"; do
    D="2022-01-28_Puck_211214_${P}"
    for FC in H52J2DMXY HLGH2BGXK; do
        for LANE in "$IN/$D/$FC"/L00*; do
            [[ -d "$LANE" ]] || continue
            N_IN=$(find "$LANE/alignment" -type f 2>/dev/null | wc -l)
            N_BK=$(find "$BK/$D/$FC/$(basename "$LANE")/alignment" -type f 2>/dev/null | wc -l)
            log "  puck$P/$FC/$(basename "$LANE")/alignment   inputs=$N_IN  backup=$N_BK"
        done
    done
done

log ""
log "  STAR summary, puck 29 H52J2DMXY L001 (confirms 60 nt and per-lane depth):"
SLOG="$BK/2022-01-28_Puck_211214_29/H52J2DMXY/L001/alignment/Puck_211214_29.star.Log.final.out"
if [[ -f "$SLOG" ]]; then
    grep -E "input reads|input read length|Uniquely mapped|too short|multiple loci" "$SLOG" | sed 's/^/      /' | tee -a "$REPORT"
else
    log "      not found"
fi
log ""
log "  cellular tagging summary (vendor statement of barcode geometry):"
CTAG="$BK/2022-01-28_Puck_211214_29/H52J2DMXY/L001/alignment/Puck_211214_29.cellular_tagging.summary.txt"
if [[ -f "$CTAG" ]]; then head -20 "$CTAG" | sed 's/^/      /' | tee -a "$REPORT"; else log "      not found"; fi
log ""
log "  adapter and polyA loss (reads in neither BAM, needed for the archive caveat):"
for S in adapter_trimming polyA_filtering; do
    F="$BK/2022-01-28_Puck_211214_29/H52J2DMXY/L001/alignment/Puck_211214_29.${S}.summary.txt"
    [[ -f "$F" ]] && { log "    $S:"; head -12 "$F" | sed 's/^/        /' | tee -a "$REPORT"; }
done

# -----------------------------------------------------------------------------
banner "BEAT 5  what are the six TCR files"
# -----------------------------------------------------------------------------
UP_LINKER="TCTTCAGCGTTCCCGAGA"

for F in "$BK"/B59_*_hTCR_tcr.csv "$BK"/TCR_*.gz; do
    [[ -f "$F" ]] || continue
    log ""
    log "--- $(basename "$F")   $(du -h "$F" | cut -f1)  ($(stat -c%s "$F") bytes)"
    log "  file: $(file -b "$F")"
    case "$F" in
        *.csv)
            log "  lines: $(wc -l < "$F")"
            log "  header: $(head -1 "$F" | cut -c1-400)"
            log "  rows 2-4:"
            sed -n '2,4p' "$F" | cut -c1-400 | sed 's/^/      /' | tee -a "$REPORT"
            log "  commas in header: $(head -1 "$F" | tr -cd ',' | wc -c)   tabs: $(head -1 "$F" | tr -cd '\t' | wc -c)"
            log "  rows with '-1':   $(grep -c -- '-1' "$F" 2>/dev/null || echo 0)"
            log "  rows with 'Puck': $(grep -c 'Puck' "$F" 2>/dev/null || echo 0)"
            log "  rows with 'TRA':  $(grep -c 'TRA' "$F" 2>/dev/null || echo 0)"
            log "  rows with 'TRB':  $(grep -c 'TRB' "$F" 2>/dev/null || echo 0)"
            log "  rows with 'CAS' (common CDR3 motif): $(grep -c 'CAS' "$F" 2>/dev/null || echo 0)"
            ;;
        *.gz)
            log "  first 400 bytes decompressed:"
            zcat "$F" 2>/dev/null | head -c 400 | cat -v | sed 's/^/      /' | tee -a "$REPORT"
            log ""
            log "  first 8 decompressed lines (160 char):"
            zcat "$F" 2>/dev/null | head -8 | cut -c1-160 | sed 's/^/      /' | tee -a "$REPORT"
            log "  line 1 first char (@ implies FASTQ): $(zcat "$F" 2>/dev/null | head -1 | cut -c1)"
            log "  tar signature at byte 257:           $(zcat "$F" 2>/dev/null | head -c 262 | tail -c 6 | tr -d '\0')"
            log ""
            log "  GEOMETRY PROBE, first 200k lines:"
            log "    UP linker found anywhere: $(zcat "$F" 2>/dev/null | head -200000 | grep -c "$UP_LINKER")"
            log "    read length distribution (top 5), assuming FASTQ:"
            zcat "$F" 2>/dev/null | head -200000 | awk 'NR%4==2{print length($0)}' \
                | sort -n | uniq -c | sort -rn | head -5 | sed 's/^/        /' | tee -a "$REPORT"
            log "    UP linker at bases 9-26 (Slide-seq R1 architecture):"
            zcat "$F" 2>/dev/null | head -200000 | awk -v L="$UP_LINKER" \
                'NR%4==2{n++; if(substr($0,9,18)==L) h++} END{if(n) printf "        %d of %d = %.1f%%\n", h, n, 100*h/n}' \
                | tee -a "$REPORT"
            ;;
    esac
done

log ""
log "  TCR GEOMETRY NOTE (we derive this ourselves, no sample sheet is coming)"
log "    rhTCRseq Fraction 2 is amplified from the same bead-barcoded cDNA pool"
log "    as Fraction 1, so the barcode-bearing read should carry the identical"
log "    split architecture: 8 nt + 18 nt UP linker + 6 nt + 9 nt UMI."
log "    A high UP-linker hit rate at bases 9-26 confirms it without metadata."
log "    Definitive check: exact-match rate of the reconstructed 14 nt barcode"
log "    against barcode_matching column 1 (observed), which is already on disk."

# -----------------------------------------------------------------------------
banner "BEAT 6  head-hash spot check on large equal-size files"
# -----------------------------------------------------------------------------
HASH_TARGETS=(
  "2022-01-28_Puck_211214_29/Puck_211214_29.matched.bam"
  "2022-01-28_Puck_211214_37/Puck_211214_37.matched.bam"
  "2022-01-28_Puck_211214_40/Puck_211214_40.matched.bam"
  "2022-01-28_Puck_211214_29/H52J2DMXY/L001/Puck_211214_29.final.bam"
  "2022-01-28_Puck_211214_29/H52J2DMXY/L001/Puck_211214_29.unmapped.bam"
)
for REL in "${HASH_TARGETS[@]}"; do
    B="$BK/$REL"; I="$IN/$REL"
    log ""
    log "--- $REL"
    [[ -f "$B" ]] || { log "    absent from BACKUP"; continue; }
    [[ -f "$I" ]] || { log "    absent from inputs (NEW)"; continue; }
    BS=$(stat -c%s "$B"); IS=$(stat -c%s "$I")
    log "    size backup=$BS  inputs=$IS  $( [[ $BS == $IS ]] && echo EQUAL || echo DIFFER )"
    BH=$(head -c 67108864 "$B" | md5sum | cut -d' ' -f1)
    IH=$(head -c 67108864 "$I" | md5sum | cut -d' ' -f1)
    log "    head64M backup=$BH"
    log "    head64M inputs=$IH  $( [[ $BH == $IH ]] && echo MATCH || echo MISMATCH )"
done

banner "DONE"
log "report:    $REPORT"
log "bam table: $BAMTSV"
log "lists:     $OUT/{only_in_backup,only_in_inputs,size_mismatch}_puck{29,37,40}.txt"
