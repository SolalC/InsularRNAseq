#!/bin/bash
#SBATCH --job-name=hfd_subcluster
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/subcluster_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/subcluster_%j.err

# Subcluster the neuronal lineages independently.
#   sbatch 10_neuronal_subclustering.sh                 # both, primary res 0.6
#   sbatch 10_neuronal_subclustering.sh excitatory      # one lineage
#   sbatch 10_neuronal_subclustering.sh inhibitory 0.8  # one lineage, other res

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} args=${*:-none}"
Rscript --vanilla "$PROJ/Code/10_neuronal_subclustering.R" "${1:-both}" "${2:-}"
echo "[$(date)] done"
