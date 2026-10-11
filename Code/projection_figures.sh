#!/bin/bash
#SBATCH --job-name=hfd_projfig
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/projfig_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/projfig_%j.err

# Projection UMAP + per-animal per-cell-type projection table, into
# Results/reporting. Housekeeping, not a pipeline step -- see the header of
# projection_figures.R. Reads the 17 GB annotated object once.
#
#   sbatch projection_figures.sh                      # >= 1 UMI
#   sbatch --export=ALL,MIN_UMI=10 projection_figures.sh   # >= 10 UMI, outputs tagged _ge10

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export MIN_UMI=${MIN_UMI:-1}
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} MIN_UMI=$MIN_UMI"
Rscript --vanilla "$PROJ/Code/projection_figures.R"
echo "[$(date)] done"
