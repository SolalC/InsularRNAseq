#!/bin/bash
#SBATCH --job-name=hfd_deg
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/deg_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/deg_%j.err

# Pseudobulk DE inside a neuronal subset, grouped by the Azimuth
# predicted.subclass written by 10_2_azimuth.sh.
#   sbatch 10_3_deg.sh                   # excitatory, res 0.2 predictions
#   sbatch 10_3_deg.sh inhibitory 0.2

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} args=${*:-none}"
Rscript --vanilla "$PROJ/Code/10_3_deg.R" "${1:-excitatory}" "${2:-0.2}"
echo "[$(date)] done"
