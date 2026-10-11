#!/bin/bash
#SBATCH --job-name=hfd_rescue
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/rescue_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/rescue_%A_%a.err

# Step 2.8: join the 2.7 BLAST hits back to cell barcodes and UMIs.
#
# One pass over each BAM's unplaced-unmapped reads (region '*', so the ~29 GB
# BAM is never decompressed end to end). Reads that 2.7 flagged are relocated by
# their position in the stream and their read name is verified, so a
# differently-ordered stream aborts instead of mis-joining.
#
# The same pass counts, over EVERY unmapped read, what fraction of barcodes are
# called cells. That is the ambient control: if rescued reads land in cells at
# the background rate they are ambient and worth nothing.
#
# See 2_8_rescue_assign.py for the filter set.

set -euo pipefail

module purge
module load samtools/1.18-gcc-12.3.0

PROJ=/scratch/user/uqschauq/Insular
MANIFEST=$PROJ/Data/sample_metadata.csv

# sample_metadata.csv is CRLF; strip \r and look the column up by header name.
SAMPLE=$(awk -F',' -v n="$SLURM_ARRAY_TASK_ID" '
  NR==1 { for (i=1;i<=NF;i++) { h=$i; gsub(/\r/,"",h); if (h=="fastq_prefix") c=i }; next }
  NR==n+1 { v=$c; gsub(/\r/,"",v); print v }' "$MANIFEST")

[[ -n $SAMPLE ]] || { echo "FATAL: no sample for array index $SLURM_ARRAY_TASK_ID"; exit 1; }
[[ -s $PROJ/Results/2_7_unmapped_blast/${SAMPLE}_hits.tsv ]] \
  || { echo "FATAL: no 2.7 hit table for $SAMPLE"; exit 1; }

echo "[$(date)] $SAMPLE: joining rescued reads to CB/UB"
python3 "$PROJ/Code/2_8_rescue_assign.py" --sample "$SAMPLE"
echo "[$(date)] done"
