#!/bin/bash
#SBATCH --job-name=hfd_percell_deg
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/percell_deg_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/percell_deg_%j.err

set -euo pipefail
module purge
module load r/4.4.2-heavy
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
# Runs: arguments if given (e.g. `sbatch THIS.sh excitatory`), otherwise both neuronal subsets in turn (pass `all` for the full object).
[ $# -eq 0 ] && set -- excitatory inhibitory
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} runs=$*"
for run in "$@"; do
  echo "[$(date)] ==== $run ===="
  Rscript --vanilla "$PROJ/Code/percell_deg_diet_sexfollowup.R" "$run"
done
echo "[$(date)] done"
