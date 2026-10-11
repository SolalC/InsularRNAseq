#!/bin/bash
#SBATCH --job-name=inst_presto
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --constraint=epyc3
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/presto_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/presto_%j.err

# Housekeeping: install presto into the user R library so Seurat's
# FindAllMarkers uses the fast rank-sum path in 6_seurat_markers.
# presto is not on CRAN; installed from a pinned GitHub commit.
#
# Built with Code/Makevars.portable (-march=x86-64-v3) instead of R's default
# -march=native: `general` mixes epyc3/epyc4/epyc5 and the user library is
# shared across all of them, so a native build on a newer node emits AVX-512
# and SIGILLs on epyc3. Pinned to --constraint=epyc3 (the oldest ISA present)
# as a second guard, so the build node can never outrun the run nodes.

set -euo pipefail

module purge
module load r-bundle-cran/2024.11-foss-2024a

PROJ=/scratch/user/uqschauq/Insular
SOFT=$PROJ/Soft
USERLIB=$HOME/R/x86_64-pc-linux-gnu-library/4.4

PRESTO_SHA=3c97180d90d330524b5b8353a342ab2f4b990d21
TARBALL=$SOFT/presto-${PRESTO_SHA:0:7}.tar.gz
SRCDIR=$SOFT/presto-$PRESTO_SHA

export R_MAKEVARS_USER=$PROJ/Code/Makevars.portable

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
echo "arch:     $(scontrol show node "$(hostname -s)" 2>/dev/null | grep -oP 'ActiveFeatures=\K\S+' || echo unknown)"
echo "R:        $(which Rscript)"
echo "user lib: $USERLIB"
echo "presto:   $PRESTO_SHA"
echo "makevars: $R_MAKEVARS_USER"

mkdir -p "$SOFT" "$USERLIB"

if [ ! -f "$TARBALL" ]; then
  echo "[$(date)] downloading source"
  curl -fsSL -o "$TARBALL" \
    "https://github.com/immunogenomics/presto/archive/${PRESTO_SHA}.tar.gz"
else
  echo "[$(date)] source tarball already present, reusing"
fi

rm -rf "$SRCDIR"
tar -xzf "$TARBALL" -C "$SOFT"
test -f "$SRCDIR/DESCRIPTION"

echo "[$(date)] compiling and installing (portable ISA)"
rm -rf "$USERLIB/presto" "$USERLIB"/00LOCK-presto*
R CMD INSTALL --library="$USERLIB" "$SRCDIR"

echo "[$(date)] verifying"
Rscript --vanilla - <<'RS'
stopifnot(requireNamespace("presto", quietly = TRUE))
so <- system.file("libs", "presto.so", package = "presto")
cat("presto version:", as.character(packageVersion("presto")), "\n")
cat("shared object: ", so, "\n")

# The build must not contain AVX-512 (zmm registers): the user library is
# shared across epyc3/4/5 and epyc3 (Zen 3) cannot execute those instructions.
zmm <- system2("bash", c("-c", shQuote(paste0("objdump -d ", shQuote(so),
               " 2>/dev/null | grep -c '%zmm' || true"))), stdout = TRUE)
zmm <- suppressWarnings(as.integer(zmm[length(zmm)]))
cat("AVX-512 (zmm) instructions:", zmm, "\n")
if (is.na(zmm) || zmm > 0) stop("presto.so contains AVX-512; not portable to epyc3 nodes")

# Smoke test: presto's fast rank-sum on a small matrix.
set.seed(1)
m <- matrix(rpois(2000, 3), nrow = 100,
            dimnames = list(paste0("g", 1:100), paste0("c", 1:20)))
g <- rep(c("A", "B"), each = 10)
res <- presto::wilcoxauc(m, g)
cat("wilcoxauc rows:", nrow(res), "| cols:", paste(colnames(res), collapse = ","), "\n")
stopifnot(nrow(res) == 200, all(c("feature", "group", "auc", "pval") %in% colnames(res)))

# Confirm Seurat actually takes the presto path in FindMarkers.
suppressPackageStartupMessages(library(Seurat))
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
so2 <- CreateSeuratObject(counts = m)
so2 <- NormalizeData(so2, verbose = FALSE)
Idents(so2) <- factor(g)
mk <- FindMarkers(so2, ident.1 = "A", ident.2 = "B",
                  test.use = "wilcox", logfc.threshold = 0, min.pct = 0,
                  verbose = FALSE)
cat("FindMarkers returned", nrow(mk), "rows via presto path -- OK\n")
RS

echo "[$(date)] done"
