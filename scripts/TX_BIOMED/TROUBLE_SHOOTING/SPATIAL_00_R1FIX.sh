#!/bin/bash
#SBATCH -J SPATIAL_00_R1FIX
#SBATCH -o /master/jlehle/WORKING/LOGS/Step00_R1Fix.o.%j.log
#SBATCH -e /master/jlehle/WORKING/LOGS/Step00_R1Fix.e.%j.log
#SBATCH -t 1-00:00:00
#SBATCH -p normal
#SBATCH --mem=32G
#SBATCH -c 8
#===============================================================================
# STEP 00: FIX READ 1 GEOMETRY (Slide-seq V2 split barcode)
#
# Raw R1 is 42 bp with a split bead barcode:
#     1-8    bead barcode part 1
#     9-26   UP linker  TCTTCAGCGTTCCCGAGA
#     27-32  bead barcode part 2
#     33-41  UMI
#     42     spare cycle
#
# Step01 ran STARsolo with --soloCBstart 1 --soloCBlen 14, which read the
# linker as barcode bases 9-14 and as the UMI. Result: 1.21% valid barcodes,
# 153 estimated cells. Confirmed by CR[9:14] == TCTTCA in 74.6% of reads.
#
# This step rewrites R1 into a synthetic 23 bp read (14 bp barcode + 9 bp UMI)
# so Step01's existing solo parameters become correct as written. R2 untouched.
#===============================================================================
set -euo pipefail

PROOT=/master/jlehle/WORKING/slide-TCR-seq-working
FQ=$PROOT/data/outputs/01_alignment/demux_fastq
OUT=$PROOT/data/outputs/01_alignment/demux_fastq_fixed
LINKER="TCTTCAGCGTTCCCGAGA"
SAMPLES=(Puck_211214_29 Puck_211214_37 Puck_211214_40)

mkdir -p "$OUT"
command -v pigz >/dev/null && ZIP="pigz -p 8" || ZIP="gzip"

for S in "${SAMPLES[@]}"; do
    IN_R1="$FQ/${S}_S*_R1_001.fastq.gz"
    IN_R1=$(ls $IN_R1)
    BASE=$(basename "$IN_R1")
    OUT_R1="$OUT/$BASE"

    echo "=== $S ==="

    # --- Guard: assert the geometry before touching anything ---
    NBAD=$(zcat "$IN_R1" | awk 'NR%4==2 && length($0)<41 {c++} NR>=4000000{exit} END{print c+0}')
    NLINK=$(zcat "$IN_R1" | awk -v L="$LINKER" 'NR%4==2{n++; if(substr($0,9,18)==L) c++}
                                                NR>=800000{exit} END{printf "%.4f", c/n}')
    echo "  reads <41bp in first 1M: $NBAD"
    echo "  exact linker at 9-26:    $NLINK"
    awk -v v="$NLINK" 'BEGIN{ if (v+0 < 0.30){ print "FATAL: linker not at 9-26"; exit 1 } }' || exit 1

    # --- Rewrite: 1-8 + 27-32 + 33-41, sequence and quality identically ---
    zcat "$IN_R1" | awk '
        NR%4==1 || NR%4==3 { print; next }
        { if (length($0) >= 41)
              print substr($0,1,8) substr($0,27,6) substr($0,33,9)
          else { print ""; short++ }
        }
        END { if (short) print "WARN: " short " short reads emitted blank" > "/dev/stderr" }
    ' | $ZIP > "$OUT_R1"

    # --- Verify: record count preserved, length now 23 ---
    NIN=$(zcat "$IN_R1"  | awk 'END{print NR/4}')
    NOU=$(zcat "$OUT_R1" | awk 'END{print NR/4}')
    LEN=$(zcat "$OUT_R1" | awk 'NR%4==2{print length($0)}' | head -100000 | sort -u | tr '\n' ',')
    echo "  reads in/out: $NIN / $NOU   lengths: $LEN"
    [ "$NIN" = "$NOU" ] || { echo "FATAL: read count changed for $S"; exit 1; }

    # --- Symlink R2 and I1 through unchanged so Step01 sees one directory ---
    for M in R2 I1; do
        SRC=$(ls $FQ/${S}_S*_${M}_001.fastq.gz 2>/dev/null) || continue
        ln -sf "$SRC" "$OUT/$(basename "$SRC")"
    done
done
echo "DONE. Point Step01 FASTQ_DIR at $OUT"
