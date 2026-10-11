#!/bin/bash
#SBATCH --job-name=hfd_ambient_dbl
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/ambient_dbl_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/ambient_dbl_%A_%a.err

# ---------------------------------------------------------------------------
# Step 4: ambient RNA removal (CellBender) + doublet calling (scDblFinder).
#
# One file, two phases, selected by $PHASE. Resources differ so much between
# them that they are passed on the sbatch command line rather than pinned here:
#
#   PHASE=cellbender   (default)  array 1-12, one GPU and one sample per task
#   PHASE=gather                  single CPU task: combine + scDblFinder
#
#   j=$(sbatch --parsable --partition=gpu_cuda --qos=gpu --gres=gpu:1 \
#              --cpus-per-task=8 --mem=100G --time=04:00:00 --array=1-12 \
#              Code/4_ambient_doublets.sh)
#   sbatch --dependency=afterok:$j --partition=general --cpus-per-task=12 \
#          --mem=250G --time=04:00:00 --export=ALL,PHASE=gather \
#          Code/4_ambient_doublets.sh
#
# WHY AN ARRAY, AND WHY IT IS STILL 4-WIDE: the `gpu` QOS caps this account at
# MaxTRESPerUser gres/gpu=4 with MaxJobsPU=4, so only 4 tasks ever run at once
# -- an array does not lift that. What it buys is that tasks start as single
# GPUs free instead of waiting for a 4-GPU block, and the gather phase runs on
# a CPU node so no GPU is held during scDblFinder.
#
# WHY CELLBENDER IS UPSTREAM OF QC: it learns the ambient profile from the
# EMPTY droplets, so it reads raw_feature_bc_matrix.h5, not the filtered
# matrix. This is also why it cannot run on the step 3 aggregate -- cellranger
# aggr emits no raw matrix at all.
#
# Each array task derives its own sample from aggregation.csv and its own
# expected-cell count from that sample's metrics_summary.csv. There is no
# shared scratch file, so 12 tasks starting at once cannot race.
# ---------------------------------------------------------------------------

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
OUTDIR=$PROJ/Results/4_ambient_doublets
CBDIR=$OUTDIR/cellbender
AGGCSV=$PROJ/Data/Aggr_v2/HFDseq_aggr_v2/outs/aggregation.csv
VENV=$PROJ/Soft/cellbender-0.4.0/venv

# FPR 0, on the upstream instruction for a cohort jointly analysed for
# differential expression: at FPR > 0 the amount of true signal removed tracks
# each library's own ambient load, which differs measurably between these 12
# animals, and that bias lands straight in the diet/sex contrasts.
FPR=0
EPOCHS=150
DROPLET_HEADROOM=20000     # headroom above the called cell count, so the
                           # barcode-rank window reaches well into empty droplets

PHASE=${PHASE:-cellbender}
mkdir -p "$CBDIR"
echo "[$(date)] phase=$PHASE host=$(hostname) job=${SLURM_JOB_ID:-none} task=${SLURM_ARRAY_TASK_ID:-na}"

# ---- gather: combine the 12 corrected matrices + scDblFinder ---------------
if [[ "$PHASE" == "gather" ]]; then
  missing=0
  while IFS=, read -r S _; do
    [[ -s "$CBDIR/${S}_filtered.h5" ]] || { echo "MISSING: ${S}_filtered.h5" >&2; missing=1; }
  done < <(tail -n +2 "$AGGCSV")
  [[ "$missing" -eq 0 ]] || { echo "FATAL: some samples produced no CellBender output" >&2; exit 1; }

  module purge
  module load r-bundle-cran/2024.11-foss-2024a
  Rscript --vanilla "$PROJ/Code/4_ambient_doublets.R"
  echo "[$(date)] gather done"
  exit 0
fi

# ---- cellbender: one sample, one GPU --------------------------------------
idx=${SLURM_ARRAY_TASK_ID:?PHASE=cellbender must run as an array}
SAMPLE=$(tail -n +2 "$AGGCSV" | sed -n "${idx}p" | cut -d, -f1 | tr -d '\r')
[[ -n "$SAMPLE" ]] || { echo "FATAL: no sample at aggregation.csv row $idx" >&2; exit 1; }

RAWH5=$PROJ/Data/Counts_v2/$SAMPLE/outs/raw_feature_bc_matrix.h5
METRICS=$PROJ/Data/Counts_v2/$SAMPLE/outs/metrics_summary.csv
[[ -f "$RAWH5" ]]   || { echo "FATAL: missing $RAWH5" >&2; exit 1; }
[[ -f "$METRICS" ]] || { echo "FATAL: missing $METRICS" >&2; exit 1; }
[[ -x "$VENV/bin/cellbender" ]] || { echo "FATAL: no cellbender at $VENV" >&2; exit 1; }

module purge
module load python/3.10.4-gcccore-11.3.0
export PYTHONNOUSERSITE=1
# shellcheck disable=SC1091
source "$VENV/bin/activate"

EXPECTED=$(python - "$METRICS" <<'PY'
import csv, sys
with open(sys.argv[1], newline='') as f:
    row = next(csv.DictReader(f))
key = next(k for k in row if 'Estimated Number of Cells' in k)
print(int(row[key].replace(',', '').replace('"', '')))
PY
)
TOTALDROPS=$(( EXPECTED + DROPLET_HEADROOM ))

echo "[$(date +%H:%M:%S)] sample=$SAMPLE gpu=${CUDA_VISIBLE_DEVICES:-none} expected=$EXPECTED total_droplets=$TOTALDROPS fpr=$FPR"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || true

cellbender remove-background \
  --input  "$RAWH5" \
  --output "$CBDIR/${SAMPLE}.h5" \
  --cuda \
  --expected-cells "$EXPECTED" \
  --total-droplets-included "$TOTALDROPS" \
  --fpr "$FPR" \
  --epochs "$EPOCHS" \
  --checkpoint "$CBDIR/${SAMPLE}_ckpt.tar.gz"

[[ -s "$CBDIR/${SAMPLE}_filtered.h5" ]] \
  || { echo "FATAL: $SAMPLE produced no _filtered.h5" >&2; exit 1; }
echo "[$(date)] done: $SAMPLE"
