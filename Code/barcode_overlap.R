#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# barcode_overlap.R  --  are excitatory-neuron barcodes unique to one animal?
# Run by barcode_overlap.sh. The computation is Barcode.r's, unchanged; this
# only makes it a job and writes the result. Housekeeping: reads the annotated
# object, changes nothing under Data/.
#
# The 10x barcode is the object's column name minus its `-<GEM group>` suffix.
# For each animal pair the matrix counts barcodes seen in both (off-diagonal);
# the diagonal counts barcodes seen in that animal and in NO other animal.
# `annot` is the object's manual call, as in Barcode.r (no Azimuth revision).
#
# Output: Results/reporting/barcode_overlap_excitatory.csv -- 12 x 12, rows and
# columns named by animal (GEM group -> mouse off the object's own metadata).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
OBJ_IN  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
OUT     <- file.path(PROJECT, "Results", "reporting", "barcode_overlap_excitatory.csv")

obj <- readRDS(OBJ_IN)
obj.n <- obj[, obj$annot == "Excitatory Neurons"]
vec.n <- colnames(obj.n)
bc.n  <- sub("-\\d+$", "", vec.n)
smp.n <- factor(as.integer(sub("^.*-", "", vec.n)))
inc.n <- 1 * (table(bc.n, smp.n) > 0)                     # barcodes x samples, 0/1
M.n   <- crossprod(inc.n)                                 # samples x samples
diag(M.n) <- colSums(inc.n[rowSums(inc.n) == 1, , drop = FALSE])

# GEM group -> animal, off the object, so the labels cannot drift
gem <- unique(data.table(gem = as.integer(obj.n$gem_group), mouse = as.character(obj.n$mouse)))
if (nrow(gem) != uniqueN(gem$gem)) stop("a GEM group maps to more than one animal")
lab <- gem$mouse[match(as.integer(colnames(M.n)), gem$gem)]
stopifnot(!anyNA(lab))
M.n <- unclass(M.n); dimnames(M.n) <- list(lab, lab)

cat("excitatory nuclei:", ncol(obj.n), "| distinct barcodes:", nrow(inc.n), "\n")
print(M.n)
cat("unique fraction per animal (Barcode.r's diag / colSums):\n")
print(round(diag(M.n) / colSums(M.n), 2))

fwrite(data.table(animal = lab, as.data.table(M.n)), OUT)
cat("wrote", OUT, "\n")
