#!/bin/bash
#SBATCH --job-name=hfd_proj_v2
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/proj_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/proj_%j.err

# Step 2.1: projection calls from the cellranger matrices.
#
# Runs 2_1_projection_validation.R, which streams only the four transgene rows
# out of each filtered_feature_bc_matrix and counts cells per projection region
# at the threshold set inside the R script (RUN <- "v2", MIN_UMI <- 1).
#
# r/4.5.2 is used because it carries data.table (1.18.0) in a non-interactive
# batch shell; the r/4.4.x modules were checked and 4.4.2 also works, 4.5.2 is
# the newest that does.

set -euo pipefail

module purge
module load r/4.5.2

PROJ=/scratch/user/uqschauq/Insular

Rscript --vanilla "$PROJ/Code/2_1_projection_validation.R"
