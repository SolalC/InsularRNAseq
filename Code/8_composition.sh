#!/bin/bash
#SBATCH --job-name=hfd_composition
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/composition_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/composition_%j.err

# Needs speckle in the user library: sbatch install_speckle.sh
# once, first.

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/8_composition.R"
echo "[$(date)] done"
