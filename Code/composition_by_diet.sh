#!/bin/bash
#SBATCH --job-name=hfd_compdiet
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --constraint=epyc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/compdiet_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/compdiet_%j.err

# fig_23 summarised by diet and by diet x sex (animals averaged), into Results/reporting. See the
# header of composition_by_diet.R. Reads a small csv only, not the object.
#
#   sbatch composition_by_diet.sh

set -euo pipefail
module purge
module load r-bundle-cran/2024.11-foss-2024a
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
Rscript --vanilla "$PROJ/Code/composition_by_diet.R"
echo "[$(date)] done"
