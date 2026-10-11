#!/bin/bash
#SBATCH --job-name=inst_scdbl
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --constraint=epyc3
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/scdbl_install_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/scdbl_install_%j.err

# Housekeeping: install scDblFinder into the user R library for 8_scdblfinder.
# Every dependency (scran, scuttle, SingleCellExperiment, bluster, xgboost,
# BiocNeighbors, BiocSingular, batchelor) is already present, so this resolves
# to scDblFinder alone. update=FALSE so nothing already in the library moves.
#
# Built on epyc3 with Code/Makevars.portable for the same reason as presto: the
# user library is shared across epyc3/4/5 and R's default -march=native would
# emit AVX-512 on a newer build node. See install_presto.sh.

set -euo pipefail

module purge
module load r-bundle-cran/2024.11-foss-2024a

PROJ=/scratch/user/uqschauq/Insular
USERLIB=$HOME/R/x86_64-pc-linux-gnu-library/4.4
export R_MAKEVARS_USER=$PROJ/Code/Makevars.portable

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
echo "arch: $(scontrol show node "$(hostname -s)" 2>/dev/null | grep -oP 'ActiveFeatures=\K\S+' || echo unknown)"
echo "user lib: $USERLIB"

Rscript --vanilla - <<'RS'
lib <- path.expand("~/R/x86_64-pc-linux-gnu-library/4.4")
.libPaths(c(lib, .libPaths()))
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")

# rtracklayer in the user library was built against OpenSSL 1.1, which this
# module does not provide (it ships OpenSSL 3), so its .so cannot dyn.load and
# ANY package importing it fails to load -- scDblFinder imports it for the ATAC
# helpers. Rebuild it against the current toolchain first. It is already broken,
# so reinstalling cannot regress anything. (rhdf5 has the same breakage and is
# NOT touched here -- reported to the user instead.)
ok <- tryCatch({ loadNamespace("rtracklayer"); TRUE }, error = function(e) FALSE)
cat("rtracklayer loadable before:", ok, "\n")
if (!ok) {
  cat("reinstalling rtracklayer against OpenSSL 3\n")
  # force=TRUE is required: the installed VERSION is current, so BiocManager
  # skips it by default -- but the binary is the thing that is broken.
  BiocManager::install("rtracklayer", lib = lib, update = FALSE, ask = FALSE,
                       force = TRUE)
  stopifnot(tryCatch({ loadNamespace("rtracklayer"); TRUE }, error = function(e) FALSE))
  cat("rtracklayer loadable after: TRUE\n")
}

if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install("scDblFinder", lib = lib, update = FALSE, ask = FALSE,
                       INSTALL_opts = "--no-test-load")
}
stopifnot(requireNamespace("scDblFinder", quietly = TRUE))
cat("scDblFinder version:", as.character(packageVersion("scDblFinder")), "\n")

# Smoke test on a small synthetic SCE, incl. the two-sample path used in step 8.
suppressPackageStartupMessages({
  library(scDblFinder); library(SingleCellExperiment); library(BiocParallel)
})
set.seed(1)
ng <- 400; nc <- 600
mu  <- rbind(matrix(rpois(ng*nc/2, 0.4), ng), matrix(rpois(ng*nc/2, 0.4), ng))
cnt <- matrix(rpois(ng*nc, 0.5), nrow = ng,
              dimnames = list(paste0("g", 1:ng), paste0("c", 1:nc)))
cnt[1:100, 1:300]   <- cnt[1:100, 1:300]   + rpois(100*300, 6)   # type A
cnt[101:200, 301:600] <- cnt[101:200, 301:600] + rpois(100*300, 6) # type B
sce <- SingleCellExperiment(list(counts = cnt))
sce$smp <- rep(c("s1", "s2"), each = nc/2)
sce <- scDblFinder(sce, samples = "smp", BPPARAM = SerialParam())
cat("classes:", paste(names(table(sce$scDblFinder.class)), table(sce$scDblFinder.class),
                      sep = "=", collapse = " "), "\n")
stopifnot(all(c("scDblFinder.score", "scDblFinder.class") %in% colnames(colData(sce))))
cat("smoke test OK\n")
RS

echo "[$(date)] done"
