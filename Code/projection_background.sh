#!/bin/bash
#SBATCH --job-name=hfd_projbg
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/projbg_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/projbg_%j.err

# Ambient-background estimate for the projection calls across thresholds 1-10,
# using non-excitatory nuclei as the negative control. See the header of
# projection_background.R. Reads the 17 GB annotated object once.
#
#   sbatch projection_background.sh

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
Rscript --vanilla "$PROJ/Code/projection_background.R"
echo "[$(date)] done"
