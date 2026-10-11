#!/bin/bash
#SBATCH --job-name=hfd_collide
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/collide_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/collide_%j.err

# Step 2.6: does any cassette sequence collide with the mouse genome?
#
# Every 90 bp window of the four transgene contigs is aligned against the
# MOUSE-ONLY index (the stock 10x refdata-gex-GRCm39-2024-A, i.e. exactly the
# mouse half of our reference with no transgene contigs to self-match). Filters
# are deliberately permissive - a superset of Cell Ranger's internal STAR
# settings, which are compiled into the binary and not recoverable from the
# pipestance - so a clean result closes the question whatever Cell Ranger uses.
#
# See the header of 2_6_cassette_mouse_collision.py for why one alignment run
# answers both directions of read diversion.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
STAR=$PROJ/Soft/cellranger-10.1.0/lib/bin/STAR
MOUSE_IDX=$PROJ/Data/Reference/refdata-gex-GRCm39-2024-A/star
PY=$PROJ/Code/2_6_cassette_mouse_collision.py
OUT=$PROJ/Results/2_6_cassette_mouse_collision

[[ -x $STAR ]]           || { echo "FATAL: no STAR at $STAR"; exit 1; }
[[ -s $MOUSE_IDX/SA ]]   || { echo "FATAL: no mouse STAR index at $MOUSE_IDX"; exit 1; }

mkdir -p "$OUT"
cd "$OUT"

echo "[$(date)] tiling cassettes into 90 bp windows"
python3 "$PY" tile --outdir "$OUT"

echo "[$(date)] aligning tiles against mouse-only index"
# Permissive on purpose: report a hit even when it is poor, so the mismatch
# distribution can be inspected rather than silently filtered away.
"$STAR" \
  --runThreadN "${SLURM_CPUS_PER_TASK:-8}" \
  --genomeDir "$MOUSE_IDX" \
  --readFilesIn "$OUT/tiles.fq" \
  --outFileNamePrefix "$OUT/tiles_" \
  --outSAMtype SAM \
  --outSAMunmapped None \
  --outFilterMultimapNmax 50 \
  --outFilterScoreMinOverLread 0.30 \
  --outFilterMatchNminOverLread 0.30 \
  --outFilterMismatchNmax 20 \
  --outFilterMismatchNoverLmax 0.3 \
  --alignIntronMax 1000000

echo "[$(date)] parsing alignments"
python3 "$PY" report --outdir "$OUT" --sam "$OUT/tiles_Aligned.out.sam"

# The SAM is regenerable and only clutters Results/; the parsed table stays.
rm -f "$OUT/tiles_Aligned.out.sam" "$OUT/tiles.fq"
rm -rf "$OUT/tiles__STARtmp"

echo "[$(date)] done"
