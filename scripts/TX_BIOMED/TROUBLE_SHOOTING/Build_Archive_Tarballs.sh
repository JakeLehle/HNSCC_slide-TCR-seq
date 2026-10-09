#!/usr/bin/env bash
#SBATCH --job-name=archive_tar
#SBATCH --partition=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/master/jlehle/WORKING/LOGS/archive_tar_%j.out
#SBATCH --error=/master/jlehle/WORKING/LOGS/archive_tar_%j.err
# =============================================================================
# Build_Archive_Tarballs.sh
# =============================================================================
# Packages everything in data/inputs into per-component tar.gz files, with
# checksums and an integrity test on each, ready for rclone to Dropbox.
#
# COVERAGE GUARANTEE
#   The script enumerates every top-level entry under data/inputs, subtracts
#   what the named components cover, and sweeps the remainder into `misc`.
#   Nothing can be dropped without appearing in the coverage table in the log.
#   This exists because an earlier version silently omitted panglaodb, the
#   marker reference behind the unified_annotation calls.
#
# WHY PER-COMPONENT AND NOT ONE TARBALL
#   The full set is ~240 GB. A single stream that size means an interrupted
#   upload restarts from zero, one corrupted byte can render everything past it
#   unreadable, and restoring one puck means unpacking all of it.
#
# WHY pigz -1
#   BAMs are BGZF, the ONT files are gzip, the BCL cbcl files are compressed.
#   That is ~95% of the bytes and gzip gains a couple of percent on it. Expect
#   tarballs at 96-99% of raw. Packaging is the point, not compression.
#
# USAGE
#   bash   Build_Archive_Tarballs.sh                      # dry run, prints plan
#   sbatch Build_Archive_Tarballs.sh --go
#   sbatch Build_Archive_Tarballs.sh --go --include-ref
#   sbatch Build_Archive_Tarballs.sh --go --rclone dropbox:Backups/SlideTCRseq
#   bash   Build_Archive_Tarballs.sh --verify             # re-test tarballs
#
# Author: Jake Lehle, Texas Biomedical Research Institute
# Project: HPV16+ HNSCC Spatial Transcriptomics (Sophia Liu collaboration)
# =============================================================================

set -uo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
INPUTS="$PROOT/data/inputs"
INV="$PROOT/data/outputs/00_inventory"
STAMP=$(date +%Y-%m)
DEST="$PROOT/data/ARCHIVE_${STAMP}"
BCL_DIR=220116_NB501164_1345_AHLGH2BGXK

GO=0; VERIFY=0; INCLUDE_REF=0; RCLONE_TARGET=""
CPUS="${SLURM_CPUS_PER_TASK:-8}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --go)          GO=1; shift ;;
        --verify)      VERIFY=1; shift ;;
        --include-ref) INCLUDE_REF=1; shift ;;
        --dest)        DEST="$2"; shift 2 ;;
        --rclone)      RCLONE_TARGET="$2"; shift 2 ;;
        -h|--help)     sed -n '10,40p' "$0"; exit 0 ;;
        *) echo "unknown arg: $1"; exit 1 ;;
    esac
done

banner() { echo; echo "===================================================================="; echo "$1"; echo "===================================================================="; }

echo "=============================================================="
echo "Build_Archive_Tarballs"
echo "  started: $(date)"
echo "  host:    $(hostname)"
echo "  job:     ${SLURM_JOB_ID:-interactive}"
echo "  cpus:    $CPUS"
echo "  source:  $INPUTS"
echo "  dest:    $DEST"
echo "=============================================================="
(( GO || VERIFY )) || echo "*** DRY RUN. Add --go to act. ***"

if command -v pigz >/dev/null 2>&1; then
    ZIP="pigz -1 -p $CPUS"; UNZIP="pigz -dc -p $CPUS"
    echo "  compressor: pigz -1 on $CPUS threads"
else
    echo "  compressor: gzip -1 (pigz not found, this will be much slower)"
    echo "              consider: conda install -c conda-forge pigz"
    ZIP="gzip -1"; UNZIP="gzip -dc"
fi

# -----------------------------------------------------------------------------
# VERIFY MODE
# -----------------------------------------------------------------------------
if (( VERIFY )); then
banner "VERIFY"
    FAIL=0
    if [[ -f "$DEST/checksums.md5" ]]; then
        echo "  checksums:"
        ( cd "$DEST" && md5sum -c checksums.md5 ) || FAIL=1
    else
        echo "  no checksums.md5 at $DEST"; FAIL=1
    fi
    echo
    echo "  tarball integrity:"
    for T in "$DEST"/*.tar.gz; do
        [[ -f "$T" ]] || continue
        printf "    %-55s " "$(basename "$T")"
        if $UNZIP "$T" | tar -tf - > /dev/null 2>&1; then echo "OK"; else echo "CORRUPT"; FAIL=1; fi
    done
    echo
    (( FAIL )) && { echo "  *** VERIFICATION FAILED ***"; exit 1; }
    echo "  all components verified"
    exit 0
fi

# -----------------------------------------------------------------------------
# COMPONENTS
#   name | parent dir | space-separated subpaths relative to parent
# -----------------------------------------------------------------------------
COMPONENTS=()
for P in 29 37 40; do
    COMPONENTS+=("slideseq_Puck_211214_${P}|$INPUTS/fastq|2022-01-28_Puck_211214_${P}")
done
COMPONENTS+=("tcr|$INPUTS|tcr")
COMPONENTS+=("bcl_${BCL_DIR}|$INPUTS/fastq|$BCL_DIR")
(( INCLUDE_REF )) && COMPONENTS+=("ref|$INPUTS|ref")

# -----------------------------------------------------------------------------
banner "COVERAGE  every top-level entry under data/inputs"
# -----------------------------------------------------------------------------
# fastq/ is covered piecewise, so account for its children rather than itself.
declare -A COVERED=()
COVERED["tcr"]=1
(( INCLUDE_REF )) && COVERED["ref"]=1

declare -A FASTQ_COVERED=()
for P in 29 37 40; do FASTQ_COVERED["2022-01-28_Puck_211214_${P}"]=1; done
FASTQ_COVERED["$BCL_DIR"]=1

MISC_PATHS=()
MISC_FASTQ=()

printf "  %-46s %-10s %s\n" "ENTRY" "SIZE" "DISPOSITION"
for entry in "$INPUTS"/*; do
    [[ -e "$entry" ]] || continue
    name=$(basename "$entry")
    sz=$(du -sh "$entry" 2>/dev/null | cut -f1)
    if [[ "$name" == "fastq" ]]; then
        printf "  %-46s %-10s %s\n" "$name/" "$sz" "covered piecewise, see below"
        for sub in "$entry"/*; do
            [[ -e "$sub" ]] || continue
            sname=$(basename "$sub")
            ssz=$(du -sh "$sub" 2>/dev/null | cut -f1)
            if [[ -n "${FASTQ_COVERED[$sname]:-}" ]]; then
                printf "    %-44s %-10s %s\n" "$sname" "$ssz" "named component"
            else
                printf "    %-44s %-10s %s\n" "$sname" "$ssz" "-> misc"
                MISC_FASTQ+=("$sname")
            fi
        done
    elif [[ -n "${COVERED[$name]:-}" ]]; then
        printf "  %-46s %-10s %s\n" "$name" "$sz" "named component"
    elif [[ "$name" == "ref" ]]; then
        printf "  %-46s %-10s %s\n" "$name" "$sz" "EXCLUDED (use --include-ref)"
    else
        printf "  %-46s %-10s %s\n" "$name" "$sz" "-> misc"
        MISC_PATHS+=("$name")
    fi
done

if (( ${#MISC_PATHS[@]} )); then
    COMPONENTS+=("misc|$INPUTS|${MISC_PATHS[*]}")
    echo
    echo "  misc component (from data/inputs):  ${MISC_PATHS[*]}"
fi
if (( ${#MISC_FASTQ[@]} )); then
    COMPONENTS+=("misc_fastq|$INPUTS/fastq|${MISC_FASTQ[*]}")
    echo "  misc_fastq component:               ${MISC_FASTQ[*]}"
fi
if (( ! ${#MISC_PATHS[@]} && ! ${#MISC_FASTQ[@]} )); then
    echo
    echo "  nothing left over, named components cover everything"
fi

# -----------------------------------------------------------------------------
banner "PRECHECK"
# -----------------------------------------------------------------------------
FAIL=0

# The README must come from the inventory diagnostic so its numbers are measured
# rather than asserted. No README, no archive.
if [[ -f "$INV/ARCHIVE_README.md" ]]; then
    echo "  OK      $INV/ARCHIVE_README.md"
else
    echo "  MISSING $INV/ARCHIVE_README.md"
    echo "          run: sbatch Diagnostic_Full_BAM_Inventory.sh --counts"
    FAIL=1
fi

for d in processed ont; do
    n=$(find "$INPUTS/tcr/$d" -type f 2>/dev/null | wc -l)
    if (( n >= 3 )); then echo "  OK      data/inputs/tcr/$d ($n files)"
    else echo "  MISSING data/inputs/tcr/$d has $n files, expected 3"
         echo "          run: sbatch Stage_New_Delivery.sh --go"; FAIL=1; fi
done

RAW_TOTAL=0
echo
for C in "${COMPONENTS[@]}"; do
    IFS='|' read -r NAME PARENT SUBS <<< "$C"
    ok=1; sz=0
    for s in $SUBS; do
        if [[ -e "$PARENT/$s" ]]; then
            sz=$(( sz + $(du -sb "$PARENT/$s" 2>/dev/null | cut -f1) ))
        else
            echo "  MISSING $PARENT/$s"; ok=0; FAIL=1
        fi
    done
    (( ok )) && printf "  OK      %-40s %s\n" "$NAME" "$(numfmt --to=iec "$sz")"
    RAW_TOTAL=$((RAW_TOTAL + sz))
done

echo
echo "  raw total:  $(numfmt --to=iec "$RAW_TOTAL")"
AVAIL=$(df -B1 --output=avail "$PROOT" 2>/dev/null | tail -1)
echo "  free space: $(numfmt --to=iec "${AVAIL:-0}")"
# tarballs land beside the source, so require room for a near-full second copy
if [[ -n "${AVAIL:-}" ]] && (( AVAIL < RAW_TOTAL )); then
    echo "  INSUFFICIENT SPACE (tarballs are written alongside the source)"; FAIL=1
fi

(( FAIL )) && { echo; echo "STOP. Fix the above first."; exit 1; }
(( GO )) || { echo; echo "Dry run complete. Nothing written."; exit 0; }

# -----------------------------------------------------------------------------
banner "BUILD"
# -----------------------------------------------------------------------------
mkdir -p "$DEST"
MANIFEST="$DEST/MANIFEST.tsv"
[[ -f "$MANIFEST" ]] || printf "component\tsource\tn_files\traw_bytes\tpacked_bytes\tratio\tmd5\n" > "$MANIFEST"
touch "$DEST/checksums.md5"

for C in "${COMPONENTS[@]}"; do
    IFS='|' read -r NAME PARENT SUBS <<< "$C"
    OUTF="$DEST/${NAME}.tar.gz"

    if [[ -f "$OUTF" ]]; then
        echo "  SKIP (exists): $(basename "$OUTF")   delete it to rebuild"
        continue
    fi

    NF=0; RAW=0
    for s in $SUBS; do
        NF=$((  NF  + $(find "$PARENT/$s" -type f 2>/dev/null | wc -l) ))
        RAW=$(( RAW + $(du -sb "$PARENT/$s" | cut -f1) ))
    done

    echo
    echo "  building $NAME"
    echo "    source: $PARENT  [$SUBS]"
    echo "    files:  $NF   raw: $(numfmt --to=iec "$RAW")"
    echo "    start:  $(date +%H:%M:%S)"

    # -h dereferences symlinks so anything staged as a link from BACKUP is
    # archived as real content. Without it the archive would hold dangling
    # links into a path that will not exist on restore.
    # SUBS is intentionally unquoted: it may hold several paths.
    if ! tar -cf - -h -C "$PARENT" $SUBS | $ZIP > "$OUTF.partial"; then
        echo "    FAILED during tar or compress"; rm -f "$OUTF.partial"; exit 1
    fi
    mv "$OUTF.partial" "$OUTF"

    PACKED=$(stat -c%s "$OUTF")
    echo "    packed: $(numfmt --to=iec "$PACKED")  ($(awk -v a="$PACKED" -v b="$RAW" 'BEGIN{printf "%.1f%%", 100*a/b}') of raw)"

    printf "    testing integrity ... "
    if $UNZIP "$OUTF" | tar -tf - > /dev/null 2>&1; then echo "OK"; else echo "CORRUPT"; exit 1; fi

    printf "    md5 ... "
    MD5=$(md5sum "$OUTF" | cut -d' ' -f1)
    echo "$MD5"
    echo "$MD5  $(basename "$OUTF")" >> "$DEST/checksums.md5"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$NAME" "$PARENT [$SUBS]" "$NF" "$RAW" "$PACKED" \
        "$(awk -v a="$PACKED" -v b="$RAW" 'BEGIN{printf "%.3f", a/b}')" "$MD5" >> "$MANIFEST"
    echo "    done:   $(date +%H:%M:%S)"
done

# -----------------------------------------------------------------------------
banner "DOCUMENTATION"
# -----------------------------------------------------------------------------
{
    echo "# Slide-TCR-seq Archive"
    echo
    echo "Packaged $(date +%Y-%m-%d) from \`$INPUTS\`."
    echo
    echo "## Restoring"
    echo
    echo '```bash'
    echo "md5sum -c checksums.md5                     # verify first, always"
    echo "pigz -dc slideseq_Puck_211214_29.tar.gz | tar -xf -"
    echo '```'
    echo
    echo "Each component extracts to its original directory name. \`tcr.tar.gz\`,"
    echo "\`misc.tar.gz\` and \`ref.tar.gz\` belong under \`data/inputs/\`; the"
    echo "slideseq, bcl and misc_fastq components belong under \`data/inputs/fastq/\`."
    echo
    echo "Symlinks were dereferenced at packing time, so every file in these"
    echo "tarballs is real content, not a link into a path that no longer exists."
    echo
    echo "## Components"
    echo
    echo "| Component | Files | Raw bytes | Packed bytes |"
    echo "|---|---|---|---|"
    awk -F'\t' 'NR>1 {printf "| `%s.tar.gz` | %s | %s | %s |\n", $1, $3, $4, $5}' "$MANIFEST"
    echo
    if (( ! INCLUDE_REF )); then
        echo "The reference genome and STAR index are **not** included; both are"
        echo "regenerable from GENCODE GRCh38 primary assembly plus gencode.v49."
        echo "Note that the Broad aligned \`matched.bam\` against GRCh38.102 with"
        echo "Ensembl contig naming, which is a different file: 193 of 194 contigs"
        echo "match by M5, and chrY differs, almost certainly PAR masking."
        echo
    fi
    echo "---"
    echo
    cat "$INV/ARCHIVE_README.md"
} > "$DEST/README.md"

for f in bam_inventory.tsv read_accounting.tsv geometry_validation.tsv tcr_inventory.tsv; do
    [[ -f "$INV/$f" ]] && cp -p "$INV/$f" "$DEST/"
done

echo "  README.md      $(wc -l < "$DEST/README.md") lines"
echo "  MANIFEST.tsv   $(( $(wc -l < "$MANIFEST") - 1 )) components"
echo "  checksums.md5  $(wc -l < "$DEST/checksums.md5") entries"

# -----------------------------------------------------------------------------
banner "SUMMARY"
# -----------------------------------------------------------------------------
ls -la "$DEST" | sed 's/^/  /'
echo
echo "  total: $(du -sh "$DEST" | cut -f1)"

# archived file count against the source, the last chance to catch an omission
SRC_N=$(find "$INPUTS" -type f -o -type l | wc -l)
ARC_N=$(awk -F'\t' 'NR>1 {s+=$3} END{print s+0}' "$MANIFEST")
echo "  files in data/inputs: $SRC_N"
echo "  files archived:       $ARC_N"
if (( ARC_N < SRC_N )); then
    echo "  NOTE: $(( SRC_N - ARC_N )) fewer. Expected only if --include-ref was omitted."
fi

# -----------------------------------------------------------------------------
banner "RCLONE"
# -----------------------------------------------------------------------------
if [[ -n "$RCLONE_TARGET" ]]; then
    if ! command -v rclone >/dev/null 2>&1; then
        echo "  rclone not found on PATH, skipping upload"
    else
        echo "  pushing to $RCLONE_TARGET"
        rclone copy "$DEST" "$RCLONE_TARGET" \
            --progress --transfers 4 --checkers 8 \
            --multi-thread-streams 4 --retries 5 --low-level-retries 10 \
            --log-file "/master/jlehle/WORKING/LOGS/rclone_archive_${SLURM_JOB_ID:-manual}.log" \
            --log-level INFO \
            || { echo "  rclone reported an error, see the log"; exit 1; }
        echo "  verifying remote against local"
        rclone check "$DEST" "$RCLONE_TARGET" --one-way \
            || { echo "  *** rclone check found differences ***"; exit 1; }
        echo "  upload verified"
    fi
else
    echo "  no --rclone target given. To push manually:"
    echo
    echo "    rclone copy $DEST dropbox:Backups/SlideTCRseq \\"
    echo "      --progress --transfers 4 --checkers 8 --retries 5"
    echo "    rclone check $DEST dropbox:Backups/SlideTCRseq --one-way"
fi

banner "DONE"
echo "  finished: $(date)"
echo
echo "  NEXT"
echo "    bash Build_Archive_Tarballs.sh --verify"
echo "    After the upload verifies, BACKUP's per-puck trees are confirmed"
echo "    byte-identical duplicates and can be removed to reclaim ~140 GB."
echo "    Note BACKUP is read-only; chmod -R u+w before removing."
