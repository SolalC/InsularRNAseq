#!/bin/bash
#SBATCH --job-name=hfd_multiqc
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/multiqc_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/multiqc_%j.err

# Step 0 (continued): aggregate the 78 FastQC reports into one MultiQC report.

set -euo pipefail

module purge
module load multiqc/1.14-foss-2022a

PROJ=/scratch/user/uqschauq/Insular
IN=$PROJ/Results/0_fastqc
OUT=$PROJ/Results/0_multiqc
mkdir -p "$OUT"

n=$(find "$IN" -name '*_fastqc.zip' -not -name 'Undetermined_*' | wc -l)
[[ $n -eq 72 ]] || echo "WARNING: $n sample FastQC reports, expected 72"

# Undetermined_* is index-assignment failure, not sample data - keep it out of
# the aggregate so the per-sample statistics are not diluted.
multiqc --force --ignore '*Undetermined*' \
  --outdir "$OUT" --filename multiqc_raw_reads.html "$IN"
echo "[$(date)] done"
