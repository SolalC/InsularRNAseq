#!/bin/bash
#SBATCH --job-name=hfd_anti
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/anti_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/anti_%A_%a.err

# Step 2.10: harvest antisense reads already aligned to the transgene contigs.
# Indexed region query on the four small contigs - a few thousand records, not a
# full BAM stream. See 2_10_antisense_harvest.py for the controls.

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

echo "[$(date)] $SAMPLE: harvesting antisense transgene reads"
python3 "$PROJ/Code/2_10_antisense_harvest.py" --sample "$SAMPLE"
echo "[$(date)] done"
