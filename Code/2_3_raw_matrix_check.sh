#!/bin/bash
#SBATCH --job-name=hfd_rawmat
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/rawmat_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/rawmat_%A_%a.err

# Step 2.3 (diagnostic): are transgene-positive nuclei being dropped by cell calling?
#
# Everything so far used filtered_feature_bc_matrix (the ~24-30k called cells).
# This walks the RAW matrix (~1.8M barcodes) and splits transgene-positive
# barcodes by whether cellranger called them cells, so a population sitting below
# the cell threshold becomes visible. Also measures ambient transgene UMIs in
# empty droplets, which is what a 1-UMI call must be specific against.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
MANIFEST=$PROJ/Data/sample_metadata.csv
COUNTS=$PROJ/Data/Counts
OUTDIR=$PROJ/Results/2_3_raw_matrix_check
SCRIPT=$PROJ/Code/2_3_raw_matrix_check.py

mkdir -p "$OUTDIR"

SAMPLE=$(tr -d '\r' < "$MANIFEST" | awk -F, -v n="$SLURM_ARRAY_TASK_ID" '
  NR==1 { for (i=1; i<=NF; i++) { if ($i=="fastq_prefix") p=i; if ($i=="fastq_ok") k=i } ; next }
  $k=="True" { c++; if (c==n) { print $p; exit } }')
[[ -n "$SAMPLE" ]] || { echo "No sample for array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }

RAW=$COUNTS/$SAMPLE/outs/raw_feature_bc_matrix
FILT=$COUNTS/$SAMPLE/outs/filtered_feature_bc_matrix

echo "[$(date)] $SAMPLE: walking raw matrix"
zcat "$RAW/matrix.mtx.gz" | tail -n +4 | python3 "$SCRIPT" \
    "$SAMPLE" "$RAW/features.tsv.gz" "$RAW/barcodes.tsv.gz" "$FILT/barcodes.tsv.gz" \
    > "$OUTDIR/${SAMPLE}_rawmat.tsv"
echo "[$(date)] $SAMPLE done"
cat "$OUTDIR/${SAMPLE}_rawmat.tsv"
