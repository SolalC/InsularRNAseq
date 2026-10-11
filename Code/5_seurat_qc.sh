#!/bin/bash
#SBATCH --job-name=hfd_qc
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/qc_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/qc_%j.err

# Peak RSS was 34.7 GB on the pre-reorder run; 120G is generous headroom for
# the extra doublet-bias pass and the subset copy.

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/5_seurat_qc.R"
echo "[$(date)] done"
