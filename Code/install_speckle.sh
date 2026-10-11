#!/bin/bash
#SBATCH --job-name=inst_speckle
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/speckle_install_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/speckle_install_%j.err

# Housekeeping: install speckle (propeller) into the user R library for
# 8_composition. Bioconductor 3.20, so BiocManager rather than a GitHub tarball.
#
# No --constraint=epyc3 and no Makevars.portable, unlike install_presto.sh:
# speckle is pure R. limma, edgeR and statmod are already in the module bundle,
# so this resolves to speckle and any pure-R Bioconductor deps.

set -euo pipefail

module purge
module load r-bundle-cran/2024.11-foss-2024a

USERLIB=$HOME/R/x86_64-pc-linux-gnu-library/4.4

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
echo "R:        $(which Rscript)"
echo "user lib: $USERLIB"

mkdir -p "$USERLIB"

Rscript --vanilla - <<'RS'
lib <- path.expand("~/R/x86_64-pc-linux-gnu-library/4.4")
.libPaths(c(lib, .libPaths()))
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
BiocManager::install("speckle", lib = lib, update = FALSE, ask = FALSE)
RS

echo "[$(date)] verifying"
Rscript --vanilla - <<'RS'
suppressPackageStartupMessages(library(speckle))
cat("speckle:", as.character(packageVersion("speckle")), "\n")

# Smoke test at the ANIMAL level: 12 animals, 6 per diet, one cell type planted
# with a shift. propeller must call it and must not call the flat one.
set.seed(1)
animals <- data.frame(sample = paste0("a", 1:12),
                      diet = rep(c("SD", "HFD"), each = 6),
                      sex  = rep(rep(c("F", "M"), each = 3), 2))
shift <- ifelse(animals$diet == "HFD", 0.30, 0.15)
cl <- unlist(lapply(seq_len(12), function(i) {
  p <- c(shift[i], 0.30, 1 - shift[i] - 0.30)
  sample(c("moved", "flat1", "flat2"), 2000, TRUE, prob = p)
}))
smp <- rep(animals$sample, each = 2000)

props <- getTransformedProps(clusters = cl, sample = smp, transform = "asin")
des <- model.matrix(~ 0 + diet + sex,
                    data = animals[match(colnames(props$TransformedProps), animals$sample), ])
ctr <- setNames(rep(0, ncol(des)), colnames(des))
ctr["dietHFD"] <- 1; ctr["dietSD"] <- -1
res <- propeller.ttest(props, des, contrasts = ctr, robust = TRUE, trend = FALSE, sort = TRUE)
print(round(res[, c("PropRatio", "Tstatistic", "P.Value", "FDR")], 4))
stopifnot(res["moved", "FDR"] < 0.05, res["moved", "PropRatio"] > 1)
stopifnot(res["flat1", "FDR"] > 0.05)
cat("smoke test OK\n")
RS

echo "[$(date)] done"
