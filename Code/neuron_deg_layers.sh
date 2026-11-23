#!/bin/bash
#SBATCH --job-name=hfd_neuron_deg
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/neuron_deg_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/neuron_deg_%j.err

# Per-cell DEG (percell_deg_diet_sexfollowup.R) then the outputs (neuron_deg_outputs.R:
# volcano, all-DEG heatmap, sex Venns, KEGG/WikiPathways bar plots) for each run given.
# Default runs: excitatory_layer (excitatory nuclei pooled by layer) and inhibitory.
# `sbatch neuron_deg_layers.sh excitatory_layer` runs one. Set SKIP_DEG=1 to redo only the
# outputs from existing DEG tables (`sbatch --export=ALL,SKIP_DEG=1 neuron_deg_layers.sh`).
# Resources from job 33304077/33305013: peak RSS 25 GB (excitatory object, the larger one);
# excitatory DEG + outputs took 3 min, inhibitory under 20 min once writing to /scratch.
# R and presto's Wilcoxon are single-threaded, so 2 cores is enough. Figures only:
#   sbatch --mem=8G --time=00:15:00 --cpus-per-task=1 --export=ALL,SKIP_DEG=1 neuron_deg_layers.sh
set -euo pipefail
module purge
module load r/4.4.2-heavy
PROJ=/scratch/user/uqschauq/Insular
export R_LIBS_USER=$HOME/R/x86_64-pc-linux-gnu-library/4.4
# Write results to /scratch, not /QRISdata: jobs stalled (no CPU) while overwriting existing
# files there. Copy the folder to QRISdata afterwards if you want it there.
export PERCELL_OUT=${PERCELL_OUT:-$PROJ/Results/percell_wilcox_diet}
mkdir -p "$PERCELL_OUT"
[ $# -eq 0 ] && set -- excitatory_layer inhibitory
echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none} runs=$* skip_deg=${SKIP_DEG:-0} out=$PERCELL_OUT"
for run in "$@"; do
  if [ "${SKIP_DEG:-0}" != "1" ]; then
    echo "[$(date)] ==== $run: DEG ===="
    Rscript --vanilla "$PROJ/Code/percell_deg_diet_sexfollowup.R" "$run"
  fi
  echo "[$(date)] ==== $run: outputs ===="
  Rscript --vanilla "$PROJ/Code/neuron_deg_outputs.R" "$run"
done
echo "[$(date)] done"
