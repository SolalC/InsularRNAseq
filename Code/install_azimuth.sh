#!/bin/bash
#SBATCH --job-name=inst_azimuth
#SBATCH --partition=general
#SBATCH --account=a_imb_cpdg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --constraint=epyc3
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/inst_azimuth_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/inst_azimuth_%j.err

# Housekeeping: install Azimuth + the Allen mouse motor-cortex reference into
# the user R library, for reference-based annotation of the excitatory
# subclusters (10_2_excitatory_azimuth).
#
# AZIMUTH IS NOT ON CRAN -- it installs from GitHub, pinned here to a commit
# rather than "master" so a rerun installs the same code. Same for SeuratData.
#
# The dependency list is long but only 7 packages are actually missing from
# r-bundle-cran/2024.11 + the user library: glmGamPoi, Signac, TFBSTools,
# JASPAR2020, EnsDb.Hsapiens.v86, shinyBS, SeuratData. Several of those are
# HUMAN genome/ATAC packages (BSgenome.Hsapiens, EnsDb.Hsapiens, JASPAR2020,
# TFBSTools, Signac). They are irrelevant to a mouse RNA analysis, but Azimuth
# declares them under Imports rather than Suggests, so the package will not load
# without them. They are installed and then never used.
#
# Built with Code/Makevars.portable and --constraint=epyc3 for the same reason
# as install_presto.sh: the user library is shared across epyc3/4/5 and a
# -march=native build on a newer node SIGILLs on the older ones.
#
# TWO REPAIRS, both found by job 27867153 (FAILED 2026-09-02):
#
#   rhdf5 IN THE USER LIBRARY IS BROKEN AND HAS BEEN FOR A WHILE. Its
#   rhdf5.so links libcrypto.so.1.1, which this module set no longer provides
#   (OpenSSL 3 -> libcrypto.so.3), so `ldd` reports it "not found" and the
#   package cannot dyn.load. Nothing noticed because nothing used it:
#   lib_read_cellbender.R reads HDF5 through hdf5r, not rhdf5. Azimuth is the
#   first thing to need it -- glmGamPoi imports HDF5Array, which lazy-loads
#   rhdf5 at install time and dies. Rhdf5lib and rhdf5 are therefore rebuilt
#   against the current toolchain before glmGamPoi is attempted, and the script
#   refuses to continue unless rhdf5 actually loads.
#
#   TFMPvalue 1.0.0 (CRAN, 2026-01-19) DECLARES R >= 4.5.0. The cluster is on
#   R 4.4.2, so available.packages() hides it and TFBSTools fails with
#   "dependency 'TFMPvalue' is not available" -- which reads like a missing
#   package but is a version wall. Pinned here to the last R-4.4-compatible
#   release, 0.0.9, from the CRAN archive.

set -euo pipefail

module purge
module load r-bundle-cran/2024.11-foss-2024a

PROJ=/scratch/user/uqschauq/Insular
SOFT=$PROJ/Soft
USERLIB=$HOME/R/x86_64-pc-linux-gnu-library/4.4

# Pinned 2026-09-02.
AZIMUTH_SHA=ad5929686e005ae91d8452187a3a1d8c0563fdd5   # v0.5.1 master
SDATA_SHA=3e51f44303069b64f5dc4d68e6a3d4a343f55c39     # seurat-data main
REF_PKG=mousecortexref.SeuratData
REF_VER=1.0.0
REF_MD5=f8c82934994c744ee55a212811b5f1fb
REF_URL=https://seurat.nygenome.org/src/contrib/${REF_PKG}_${REF_VER}.tar.gz

export R_MAKEVARS_USER=$PROJ/Code/Makevars.portable
export R_LIBS_USER=$USERLIB

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
echo "arch:     $(scontrol show node "$(hostname -s)" 2>/dev/null | grep -oP 'ActiveFeatures=\K\S+' || echo unknown)"
echo "R:        $(which Rscript)"
echo "user lib: $USERLIB"
echo "azimuth:  $AZIMUTH_SHA"
echo "makevars: $R_MAKEVARS_USER"

mkdir -p "$SOFT" "$USERLIB"

## ---- 1. dependencies --------------------------------------------------------

echo "[$(date)] clearing stale install locks"
rm -rf "$USERLIB"/00LOCK-*

# Repair 1: the HDF5 stack that glmGamPoi needs. Each step is its OWN Rscript
# process on purpose -- once dyn.load of rhdf5.so has failed in a session, R
# will not retry it there, so a rebuild verified inside the same session always
# looks like it failed (job 27867427 died exactly that way, on an rhdf5 that had
# in fact been rebuilt correctly).
loads_rhdf5 () { Rscript --vanilla -e \
  'quit(status = if (requireNamespace("rhdf5", quietly = TRUE)) 0L else 1L)' >/dev/null 2>&1; }

if loads_rhdf5; then
  echo "[$(date)] rhdf5 loads -- no rebuild needed"
else
  echo "[$(date)] rhdf5 does not load -- rebuilding Rhdf5lib + rhdf5"
  Rscript --vanilla -e 'lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
    options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8L, timeout = 3600)
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib)
    for (p in c("Rhdf5lib", "rhdf5"))
      BiocManager::install(p, lib = lib, ask = FALSE, update = FALSE, force = TRUE)'
  if loads_rhdf5; then
    echo "[$(date)] rhdf5 rebuilt and loading"
  else
    echo "rhdf5 still will not load after rebuild:"
    ldd "$USERLIB/rhdf5/libs/rhdf5.so" | grep "not found" || true
    exit 1
  fi
fi

echo "[$(date)] installing missing dependencies"
Rscript --vanilla - <<'RS'
lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8L,
        timeout = 3600, warn = 1)
have <- function(p) requireNamespace(p, quietly = TRUE)
if (!have("BiocManager")) install.packages("BiocManager", lib = lib)

## repair 2: TFBSTools needs TFMPvalue, whose current CRAN release requires
## R >= 4.5. Take the last release that still supports R 4.4.
if (!have("TFMPvalue")) {
  url <- "https://cloud.r-project.org/src/contrib/Archive/TFMPvalue/TFMPvalue_0.0.9.tar.gz"
  cat("installing pinned TFMPvalue 0.0.9 (current CRAN 1.0.0 needs R >= 4.5)\n")
  install.packages(url, repos = NULL, type = "source", lib = lib)
  if (!have("TFMPvalue")) stop("TFMPvalue 0.0.9 failed to install")
}

cran <- c("Signac", "shinyBS")
bioc <- c("glmGamPoi", "TFBSTools", "JASPAR2020", "EnsDb.Hsapiens.v86")

need <- cran[!vapply(cran, have, logical(1))]
if (length(need)) { cat("CRAN:", paste(need, collapse = ", "), "\n")
                    install.packages(need, lib = lib) }

need <- bioc[!vapply(bioc, have, logical(1))]
if (length(need)) {
  if (!have("BiocManager")) install.packages("BiocManager", lib = lib)
  cat("Bioconductor", as.character(BiocManager::version()), ":",
      paste(need, collapse = ", "), "\n")
  BiocManager::install(need, lib = lib, ask = FALSE, update = FALSE, Ncpus = 8L)
}

all_deps <- c(cran, bioc, "TFMPvalue", "rhdf5", "HDF5Array")
missing <- all_deps[!vapply(all_deps, have, logical(1))]
if (length(missing)) stop("dependencies still missing: ", paste(missing, collapse = ", "))
cat("all dependencies present\n")
RS

## ---- 2. SeuratData + Azimuth from pinned GitHub commits ---------------------

install_gh () {   # $1 = owner/repo  $2 = sha  $3 = human name
  local tarball="$SOFT/$3-${2:0:7}.tar.gz" srcdir
  if [ ! -f "$tarball" ]; then
    echo "[$(date)] downloading $3 $2"
    curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors \
      -o "$tarball" "https://github.com/$1/archive/$2.tar.gz"
  else
    echo "[$(date)] $3 tarball already present, reusing"
  fi
  # awk, not `head -1`: head closes the pipe after one line, tar takes SIGPIPE,
  # and `set -o pipefail` turns that into a job failure (job 27868596).
  srcdir=$(tar -tzf "$tarball" | awk -F/ 'NR == 1 { print $1 }')
  rm -rf "${SOFT:?}/$srcdir"
  tar -xzf "$tarball" -C "$SOFT"
  test -f "$SOFT/$srcdir/DESCRIPTION"
  rm -rf "$USERLIB/$3" "$USERLIB"/00LOCK-*
  R CMD INSTALL --library="$USERLIB" "$SOFT/$srcdir"
}

install_gh satijalab/seurat-data "$SDATA_SHA" SeuratData
install_gh satijalab/azimuth    "$AZIMUTH_SHA" Azimuth

## ---- 3. the reference -------------------------------------------------------

REF_TARBALL=$SOFT/${REF_PKG}_${REF_VER}.tar.gz
if [ ! -f "$REF_TARBALL" ]; then
  echo "[$(date)] downloading $REF_PKG $REF_VER (57 MB)"
  curl -fsSL --retry 5 --retry-delay 5 --retry-all-errors -o "$REF_TARBALL" "$REF_URL"
fi
echo "[$(date)] verifying checksum"
echo "$REF_MD5  $REF_TARBALL" | md5sum -c -

rm -rf "$USERLIB/$REF_PKG" "$USERLIB"/00LOCK-*
R CMD INSTALL --library="$USERLIB" "$REF_TARBALL"

## ---- 4. verify --------------------------------------------------------------

echo "[$(date)] verifying"
Rscript --vanilla - <<'RS'
lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({ library(Seurat); library(SeuratData); library(Azimuth) })
cat("Seurat:     ", as.character(packageVersion("Seurat")), "\n")
cat("Azimuth:    ", as.character(packageVersion("Azimuth")), "\n")
cat("SeuratData: ", as.character(packageVersion("SeuratData")), "\n")

# Newly compiled shared objects must not carry AVX-512: the user library is
# shared across epyc3/4/5 and epyc3 (Zen 3) cannot execute those instructions.
for (p in c("Azimuth", "glmGamPoi", "TFBSTools", "Signac")) {
  d <- system.file("libs", package = p)
  if (!nzchar(d)) next
  for (so in list.files(d, pattern = "\\.so$", full.names = TRUE)) {
    n <- suppressWarnings(as.integer(system2("bash", c("-c",
      shQuote(sprintf("objdump -d %s 2>/dev/null | grep -c '%%zmm' || true", shQuote(so)))),
      stdout = TRUE)))
    cat(sprintf("  %-12s %-40s zmm=%s\n", p, basename(so), n))
    if (is.na(n) || n > 0) stop(basename(so), " contains AVX-512; not portable to epyc3")
  }
}

# The reference must be installed and readable, and must carry the annotation
# levels 10_2 will transfer.
inst <- InstalledData()
cat("installed references:", paste(rownames(inst), collapse = ", "), "\n")
stopifnot("mousecortexref.SeuratData" %in% rownames(inst))
refdir <- system.file("azimuth", package = "mousecortexref.SeuratData")
cat("reference dir:", refdir, "\n")
print(list.files(refdir))
stopifnot(all(c("ref.Rds", "idx.annoy") %in% list.files(refdir)))

ref <- readRDS(file.path(refdir, "ref.Rds"))
cat("reference: ", ncol(ref), " cells | assays: ", paste(Assays(ref), collapse = ", "),
    " | reductions: ", paste(Reductions(ref), collapse = ", "), "\n", sep = "")
cat("annotation columns available for transfer:\n")
print(names(slot(ref, "meta.data")))
for (a in intersect(c("class", "subclass", "cluster", "cross_species_cluster"),
                    names(slot(ref, "meta.data"))))
  cat("  ", a, ": ", length(unique(slot(ref, "meta.data")[[a]])), " levels -- ",
      paste(head(sort(unique(as.character(slot(ref, "meta.data")[[a]]))), 8),
            collapse = ", "), "\n", sep = "")
cat("OK\n")
RS

echo "[$(date)] done"
