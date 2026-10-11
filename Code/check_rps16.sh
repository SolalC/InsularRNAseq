#!/bin/bash
#SBATCH --job-name=hfd_chk_rps16
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/chk_rps16_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/chk_rps16_%j.err

# Diagnostic for the single DE hit from 10_3_excitatory_deg.R: Rps16 in L5 ET.
set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/check_rps16.R"
echo "[$(date)] done"
