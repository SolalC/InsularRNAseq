#!/bin/bash
#SBATCH --job-name=hfd_fastqc
#SBATCH --partition=general
#SBATCH --array=1-72
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/fastqc_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/fastqc_%A_%a.err

# Step 0: FastQC on the raw reads. Numbered 0 because it belongs BEFORE the
# reference build and count in pipeline order, even though it is being run after
# them here.
#
# R1 (28 bp: 16 bp barcode + 12 bp UMI) and R2 (90 bp cDNA) only. I1/I2 are
# sample indices and say nothing about read quality that matters for mapping.
#
# The instrument is a NextSeq 2000 (read names start VH00807), which uses
# two-colour chemistry: G is encoded as "no signal", so dark or low-signal cycles
# produce spurious poly-G tails. Those reads fail to align and would show up here
# as an overrepresented sequence and a rising %G towards the read end. That is
# the specific artefact worth looking for, given the unmapped records carry
# `uT:A:1` (alignment too short) with otherwise perfect partial alignments.

set -euo pipefail

module purge
module load fastqc/0.11.9-java-11

PROJ=/scratch/user/uqschauq/Insular
OUT=$PROJ/Results/0_fastqc
mkdir -p "$OUT"

# Deterministic file list: sorted, so the array index always names the same file.
# Undetermined_S0_* is excluded: those are reads that failed index assignment,
# not sample data, AND all three pools name them identically, so FastQC would
# write all three onto the same output basename and silently keep only the last.
# (Index hopping is assessed separately from the per-pool Index_Hopping_Counts.csv.)
mapfile -t FILES < <(find "$PROJ/Data" -maxdepth 2 \
  \( -name '*_R1_*.fastq.gz' -o -name '*_R2_*.fastq.gz' \) \
  -not -name 'Undetermined_*' | sort)

N=${#FILES[@]}
[[ $N -eq 72 ]] || { echo "FATAL: expected 72 sample R1/R2 files, found $N - fix --array"; exit 1; }

FQ=${FILES[$((SLURM_ARRAY_TASK_ID - 1))]}
[[ -s $FQ ]] || { echo "FATAL: missing $FQ"; exit 1; }

echo "[$(date)] $(basename "$FQ") ($(du -h "$FQ" | cut -f1))"
fastqc --threads "${SLURM_CPUS_PER_TASK:-2}" --outdir "$OUT" --quiet "$FQ"
echo "[$(date)] done"
