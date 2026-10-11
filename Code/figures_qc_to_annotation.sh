#!/bin/bash
#SBATCH --job-name=hfd_figures
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/figures_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/figures_%j.err

# Every figure from QC through annotation, into Results/reporting.
#   sbatch figures_qc_to_annotation.sh          # everything
#   sbatch figures_qc_to_annotation.sh tables   # sections A+B only, no object

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} args=${*:-none}"
Rscript --vanilla "$PROJ/Code/figures_qc_to_annotation.R" "${1:-all}"
echo "[$(date)] done"
