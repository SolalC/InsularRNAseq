#!/bin/bash
#SBATCH --job-name=hfd_unmapped
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/unmapped_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/unmapped_%A_%a.err

# Step 2.2 (diagnostic): do the UNMAPPED reads contain transgene sequence?
#
# Tests whether the near-total absence of EBFP2 (and of tdTomato in 5 mice) is a
# mapping failure - reads present but not aligned - or a genuine absence of the
# sequence from the libraries.
#
# Every unplaced-unmapped read is fetched via the BAM index (region '*', so the
# 29 GB BAM is not decompressed end to end) and scanned for 31-mers that are
# unique to one of the four transgene contigs. K-mers shared between EGFP and
# EBFP2 are excluded, so a hit is diagnostic for a single contig.

set -euo pipefail

module load samtools/1.18-gcc-12.3.0

PROJ=/scratch/user/uqschauq/Insular
MANIFEST=$PROJ/Data/sample_metadata.csv
COUNTS=$PROJ/Data/Counts
OUTDIR=$PROJ/Results/2_2_unmapped_check
MATCHER=$PROJ/Code/2_2_unmapped_kmer_match.py

mkdir -p "$OUTDIR"

# Sample list: fastq_prefix column of rows with fastq_ok==True. The file has
# CRLF line endings, so \r is stripped and columns are found by header name.
SAMPLE=$(tr -d '\r' < "$MANIFEST" | awk -F, -v n="$SLURM_ARRAY_TASK_ID" '
  NR==1 { for (i=1; i<=NF; i++) { if ($i=="fastq_prefix") p=i; if ($i=="fastq_ok") k=i } ; next }
  $k=="True" { c++; if (c==n) { print $p; exit } }')

if [[ -z "$SAMPLE" ]]; then
  echo "No sample for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

BAM=$COUNTS/$SAMPLE/outs/possorted_genome_bam.bam
OUT=$OUTDIR/${SAMPLE}_unmapped_transgene.tsv

echo "[$(date)] $SAMPLE: scanning unmapped reads in $BAM"

samtools view -@ 4 "$BAM" '*' \
  | python3 "$MATCHER" \
  > "$OUT"

echo "[$(date)] $SAMPLE done"
cat "$OUT"
