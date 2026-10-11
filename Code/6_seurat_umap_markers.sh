#!/bin/bash
#SBATCH --job-name=hfd_umap_markers
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/umap_markers_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/umap_markers_%j.err

# Both halves in one job. Pre-reorder these ran separately at 180G (umap) and
# 135 GB measured (markers, 51 min), so 250G covers the merged peak.

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/6_seurat_umap_markers.R"
echo "[$(date)] done"
