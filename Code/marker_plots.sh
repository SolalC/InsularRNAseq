#!/bin/bash
#SBATCH --job-name=hfd_marker_plots
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/marker_plots_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/marker_plots_%j.err

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/marker_plots.R" "${1:-10}" "${2:-200}"
echo "[$(date)] done"
