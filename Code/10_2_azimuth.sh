#!/bin/bash
#SBATCH --job-name=hfd_azimuth
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/azimuth_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/azimuth_%j.err

# Azimuth mapping of a neuronal subset onto the Allen mouse MOp reference.
# Needs Code/install_azimuth.sh to have run first.
#   sbatch 10_2_azimuth.sh                   # excitatory, res 0.2
#   sbatch 10_2_azimuth.sh inhibitory 0.2

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} args=${*:-none}"
Rscript --vanilla "$PROJ/Code/10_2_azimuth.R" "${1:-excitatory}" "${2:-0.2}" "${3:-}"
echo "[$(date)] done"
