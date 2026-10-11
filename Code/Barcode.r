library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)

PROJECT <- "/scratch/user/uqschauq/Insular"
OBJ_IN  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
obj <- readRDS(OBJ_IN)
vec <- colnames(obj)
bc  <- sub("-\\d+$", "", vec)
smp <- factor(as.integer(sub("^.*-", "", vec)))
inc <- 1 * (table(bc, smp) > 0)                     # barcodes x samples, 0/1
M   <- crossprod(inc)                               # samples x samples
diag(M) <- colSums(inc[rowSums(inc) == 1, , drop = FALSE])
round(diag(M)/colSums(M), 2)


# Same but only excitatory neuron:
obj.n <- obj[,obj$annot == 'Excitatory Neurons']
vec.n <- colnames(obj.n)
bc.n  <- sub("-\\d+$", "", vec.n)
smp.n <- factor(as.integer(sub("^.*-", "", vec.n)))
inc.n <- 1 * (table(bc.n, smp.n) > 0)                     # barcodes x samples, 0/1
M.n   <- crossprod(inc.n)                               # samples x samples
diag(M.n) <- colSums(inc.n[rowSums(inc.n) == 1, , drop = FALSE])

round(diag(M.n)/colSums(M.n), 2)
 