#!/bin/bash
#SBATCH --job-name=cb_smoke
#SBATCH --account=a_imb_cpdg
#SBATCH --partition=gpu_cuda
#SBATCH --qos=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=100G --time=01:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/cb_smoke_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/cb_smoke_%j.err

# Housekeeping: end-to-end validation of the CellBender venv on a GPU node
# before committing 12 x ~35 min of array time.
#
# It exists because two separate failures were each only visible at the END of
# a full run: an epyc4-built venv died with SIGILL on the epyc3 GPU nodes, and
# torch 2.13 broke checkpoint saving so the run failed 32 minutes in at
# load_or_compute_posterior_and_save(). Both would have been caught in 3
# minutes by this. Five epochs is enough to exercise train -> checkpoint ->
# posterior -> h5 output; the numbers it produces are meaningless.

set -euo pipefail
PROJ=/scratch/user/uqschauq/Insular
VENV=$PROJ/Soft/cellbender-0.4.0/venv
SAMPLE=Ech2FHFD
# Output is PRESERVED on purpose: a CellBender run is ~15 min even at 5 epochs,
# and keeping the h5 lets the R reader be iterated against a real file without
# paying for another one. Delete Results/4_ambient_doublets/.smoke by hand.
WORK=$PROJ/Results/4_ambient_doublets/.smoke
rm -rf "$WORK"; mkdir -p "$WORK"

echo "node=$(hostname -s) arch=$(scontrol show node "$(hostname -s)" | grep -oP 'ActiveFeatures=\K\S+')"
module purge; module load python/3.10.4-gcccore-11.3.0
export PYTHONNOUSERSITE=1
source "$VENV/bin/activate"
python -c "import torch; print('torch', torch.__version__, '| cuda', torch.cuda.is_available(), '|', torch.cuda.get_device_name(0))"

echo "--- 5-epoch run on $SAMPLE ---"
cellbender remove-background \
  --input  "$PROJ/Data/Counts_v2/$SAMPLE/outs/raw_feature_bc_matrix.h5" \
  --output "$WORK/${SAMPLE}.h5" \
  --cuda --expected-cells 25362 --total-droplets-included 45362 \
  --fpr 0 --epochs 5 --checkpoint "$WORK/${SAMPLE}_ckpt.tar.gz"

echo "--- CHECK 1: checkpoint tarball written? ---"
[[ -s "$WORK/${SAMPLE}_ckpt.tar.gz" ]] \
  && echo "PASS checkpoint: $(du -h "$WORK/${SAMPLE}_ckpt.tar.gz" | cut -f1)" \
  || { echo "FAIL: no checkpoint tarball -- torch pin still wrong" >&2; exit 1; }

echo "--- CHECK 2: filtered h5 written? ---"
[[ -s "$WORK/${SAMPLE}_filtered.h5" ]] \
  && echo "PASS filtered h5: $(du -h "$WORK/${SAMPLE}_filtered.h5" | cut -f1)" \
  || { echo "FAIL: no _filtered.h5" >&2; exit 1; }

echo "--- CHECK 3: can Seurat's Read10X_h5 parse it, and are the transgenes there? ---"
module purge; module load r-bundle-cran/2024.11-foss-2024a
Rscript --vanilla -e '
suppressPackageStartupMessages({library(Matrix); library(hdf5r)})
source("/scratch/user/uqschauq/Insular/Code/lib_read_cellbender.R")
f <- commandArgs(TRUE)[1]
h <- H5File$new(f, mode="r"); cat("h5 top-level groups:", paste(names(h), collapse=", "), "\n"); h$close_all()
m <- read_cellbender_h5(f)
cat("dims:", nrow(m), "features x", ncol(m), "cells\n")
tg <- c("EGFP","tdTomato","EBFP2","AAV_shared")
cat("transgenes present:", paste(intersect(tg, rownames(m)), collapse=", "), "\n")
stopifnot(all(tg %in% rownames(m)))
cat("counts are integer-valued:", all(m@x == round(m@x)), "\n")
cat("barcode example:", colnames(m)[1], "\n")
cat("PASS read\n")' "$WORK/${SAMPLE}_filtered.h5"

echo "ALL SMOKE CHECKS PASSED"
