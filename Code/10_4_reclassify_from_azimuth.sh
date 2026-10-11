#!/bin/bash
#SBATCH --job-name=hfd_reclass
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/reclass_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/reclass_%j.err

# Revise `annot` where the Azimuth mapping disagrees with the manual call.
# Reads only the 10_2 prediction tables; the annotated object is not touched.
set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/10_4_reclassify_from_azimuth.R"
echo "[$(date)] done"
