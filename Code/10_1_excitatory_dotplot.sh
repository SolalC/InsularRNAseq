#!/bin/bash
#SBATCH --job-name=hfd_exc_dotplot
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/exc_dotplot_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/exc_dotplot_%j.err

# Marker dotplot for the excitatory subclustering at the 16-cluster solution.
#   sbatch 10_1_excitatory_dotplot.sh          # res 0.2, 5 genes per cluster
#   sbatch 10_1_excitatory_dotplot.sh 0.2 5    # explicit

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} args=${*:-none}"
Rscript --vanilla "$PROJ/Code/10_1_excitatory_dotplot.R" "${1:-0.2}" "${2:-5}"
echo "[$(date)] done"
