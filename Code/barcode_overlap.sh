#!/bin/bash
#SBATCH --job-name=hfd_bcover
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --constraint=epyc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/bcover_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/bcover_%j.err

# Excitatory-neuron barcode overlap between animals -> a 12 x 12 csv in
# Results/reporting. See the header of barcode_overlap.R.
#
#   sbatch barcode_overlap.sh

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
Rscript --vanilla "$PROJ/Code/barcode_overlap.R"
echo "[$(date)] done"
