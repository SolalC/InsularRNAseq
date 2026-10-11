#!/bin/bash
#SBATCH --job-name=hfd_anno
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/anno_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/anno_%j.err

set -euo pipefail

module purge
module load r-bundle-cran/2024.11-foss-2024a

PROJ=/scratch/user/uqschauq/Insular

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/7_seurat_annotation.R"
echo "[$(date)] done"
