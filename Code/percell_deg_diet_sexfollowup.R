# Per-cell DEG, HFD vs SD, with a sex-specific follow-up. Seurat-native.
#
# Usage:   Rscript percell_deg_diet_sexfollowup.R [all | excitatory | inhibitory]
#   all          annotated object, groups = `annot`        -> .../all_celltypes
#   excitatory   excitatory subset, groups = Azimuth        -> .../excitatory_subclasses
#   inhibitory   inhibitory subset, groups = Azimuth        -> .../inhibitory_subclasses
#   For the subsets, `annot` is replaced by the per-nucleus Azimuth predicted.subclass,
#   keeping only subclasses whose nuclei are modally the lineage's own predicted.class
#   (as in pseudobulk_deg_seurat_diet.R and 10_3_deg.R).
#
# Step 1  Pooled:  FindMarkers(test.use = "wilcox") per `annot`, HFD vs SD nuclei, sexes pooled.
#                  Hit = Seurat's p_val_adj (Bonferroni over all genes in the object) < 0.05
#                  and |avg_log2FC| > LFC_CUT.
# Step 2  Per sex: same test within F and within M, on the step-1 hits of that cell type only.
# Joined:          pseudobulk DESeq2 (pseudobulk_deg_seurat_diet.R, sex ignored) and, if present,
#                  the edgeR diet x sex interaction (Code/9_pseudobulk_deg.R, model B).
# Added:           per-mouse means of the hits (mouse_means_hits_<ct>.csv, used by
#                  percell_deg_plots.R).
#
# The per-cell p-values treat nuclei as independent replicates. They are not: the unit of
# replication is the mouse (6 vs 6 pooled, 3 vs 3 per sex). Read them as a ranking, and use the
# pseudobulk p-value to judge whether an effect holds across animals.
#
# Untested. Wilcoxon is fast only with `presto` installed; without it the large neuronal types
# take hours. Peak memory is set by the largest type (Excitatory, ~137k nuclei); if it is
# still killed, request more or set MAX_CELLS (below).

suppressPackageStartupMessages({
  library(Seurat); library(data.table)
})

ROOT        <- "/QRISdata/Q9732/scratch"
PROJECT     <- "/scratch/user/uqschauq/Insular"   # where scripts 9 and 10_x write
args        <- commandArgs(trailingOnly = TRUE)
RUN         <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "all"
LINEAGE_CLASS <- c(excitatory = "Glutamatergic", inhibitory = "GABAergic")
if (!RUN %in% c("all", names(LINEAGE_CLASS)))
  stop("argument must be one of: all, ", paste(names(LINEAGE_CLASS), collapse = ", "))
AZIMUTH_RES <- 0.2    # resolution of the Azimuth predictions written by 10_2_azimuth.R
SUBDIR      <- if (RUN == "all") "all_celltypes" else paste0(RUN, "_subclasses")
first_existing <- function(x) c(Filter(file.exists, x), x[1])[1]

IN       <- if (RUN == "all") file.path(ROOT, "Data", "Seurat", "HFDseq_v2_annotated.rds") else
              file.path(ROOT, "Data", "Seurat", sprintf("HFDseq_v2_sub_%s.rds", RUN))
OUTDIR   <- file.path(ROOT, "Results", "percell_wilcox_diet", SUBDIR)
PB_DIR   <- file.path(ROOT, "Results", "pseudobulk_seurat_diet", SUBDIR)
# edgeR diet x sex interaction: 9_pseudobulk_deg.R (all) or 10_3_deg.R (subsets).
PB_B_REL <- if (RUN == "all") file.path("9_pseudobulk_deg", "pseudobulk_de_all.csv.gz") else
              file.path(sprintf("10_3_%s_deg", RUN), "pseudobulk_de_all.csv.gz")
PB_B_FILE <- first_existing(file.path(c(ROOT, PROJECT), "Results", PB_B_REL))
PRED_FILE <- if (RUN == "all") NA_character_ else
  first_existing(file.path(c(ROOT, PROJECT), "Results", sprintf("10_2_%s_azimuth", RUN),
                           sprintf("azimuth_predictions_res%s.rds", AZIMUTH_RES)))
TRANSGENES  <- c("EGFP", "tdTomato", "EBFP2", "AAV-shared")
MIN_CELLS   <- 20     # nuclei each diet needs (within the subset) for a test to run
MIN_PCT     <- 0.01   # Seurat v5 defaults for the pooled test
LOGFC_MIN   <- 0.1
MAX_CELLS   <- Inf    # max.cells.per.ident; set e.g. 2000 to cap the cell-number effect on p
BONF_CUT    <- 0.05
LFC_CUT     <- 0.5    # |avg_log2FC| a hit must exceed (log2, so 1.41-fold)
SEED        <- 42
DROP_TYPES  <- "Unidentified"   # `annot` levels left out of every analysis

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
msg  <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("presto", quietly = TRUE))
  msg("WARNING: presto not installed; FindMarkers falls back to the slow Wilcoxon")

## ---- object ------------------------------------------------------------------------
msg("run = ", RUN, "; reading ", basename(IN))
obj <- readRDS(IN)
stopifnot(all(c("sample", "diet", "sex") %in% names(obj@meta.data)))
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]], search = "counts")) > 1) obj <- JoinLayers(obj)
# Keep only what the tests use (RNA `data`): drops counts, scale.data, other assays,
# reductions and graphs, which otherwise sit in memory for the whole run.
obj <- DietSeurat(obj, assays = "RNA", layers = "data", dimreducs = NULL, graphs = NULL)
invisible(gc())
if (RUN != "all") {
  pred <- as.data.table(readRDS(PRED_FILE))
  i <- match(colnames(obj), pred$barcode)
  stopifnot(!anyNA(i))
  obj$annot <- as.character(pred$predicted.subclass[i])
  pc    <- as.character(pred$predicted.class[i])
  modal <- data.table(annot = obj$annot, pc = pc)[, .N, by = .(annot, pc)][order(-N)][, .SD[1], by = annot]
  keep  <- modal[pc == LINEAGE_CLASS[[RUN]], annot]
  msg("dropping ", sum(!obj$annot %in% keep), " nuclei outside ", LINEAGE_CLASS[[RUN]], " subclasses")
  obj <- subset(obj, cells = colnames(obj)[obj$annot %in% keep])
  rm(pred); invisible(gc())
}
stopifnot("annot" %in% names(obj@meta.data))
obj$annot  <- as.character(obj$annot)
obj$sample <- as.character(obj$sample)
obj$diet   <- factor(as.character(obj$diet), levels = c("SD", "HFD"))
obj$sex    <- factor(as.character(obj$sex),  levels = c("F", "M"))
FEATURES   <- setdiff(rownames(obj), TRANSGENES)
mice <- unique(as.data.table(obj[[]])[, .(sample, diet, sex)])
stopifnot(!anyDuplicated(mice$sample))
print(mice[order(diet, sex, sample)])

## ---- helpers -----------------------------------------------------------------------
# Wilcoxon HFD vs SD on `cells`. Returns NULL if either diet has < MIN_CELLS nuclei.
# Cells are passed as ident.1 / ident.2 (Seurat v5 accepts cell names there), so the
# object is never copied by subset().
wilcox_diet <- function(cells, features, min.pct, logfc) {
  diet <- obj$diet[match(cells, colnames(obj))]
  n <- table(diet)
  if (any(n < MIN_CELLS)) return(NULL)
  r <- FindMarkers(obj, ident.1 = cells[diet == "HFD"], ident.2 = cells[diet == "SD"],
                   test.use = "wilcox", features = features, min.pct = min.pct,
                   logfc.threshold = logfc, max.cells.per.ident = MAX_CELLS,
                   random.seed = SEED, verbose = FALSE)
  r <- data.table(gene = rownames(r), r)
  attr(r, "n") <- n
  invisible(gc())
  r
}

# Per-mouse mean expression on the linear scale, mean(expm1(data)), as AverageExpression does.
# Computed by hand to avoid AverageExpression renaming the groups.
mouse_means <- function(cells, genes) {
  x   <- expm1(LayerData(obj, assay = "RNA", layer = "data")[genes, cells, drop = FALSE])
  smp <- obj$sample[match(cells, colnames(obj))]
  g   <- Matrix::sparse.model.matrix(~ 0 + factor(smp))
  av  <- as.matrix(x %*% g) / matrix(Matrix::colSums(g), length(genes), ncol(g), byrow = TRUE)
  colnames(av) <- levels(factor(smp))
  d <- as.data.table(av, keep.rownames = "gene")
  d <- melt(d, id.vars = "gene", variable.name = "sample", value.name = "mean_expr")
  d[, sample := as.character(sample)]
  merge(d, mice, by = "sample")
}

read_pb <- function(ct) {
  f <- file.path(PB_DIR, paste0("DEG_HFDvsSD_", safe(ct), ".csv"))
  if (!file.exists(f)) return(NULL)
  fread(f)[, .(gene, pb_log2FC = avg_log2FC, pb_p_val = p_val, pb_FDR = FDR)]
}

pbB <- NULL
if (file.exists(PB_B_FILE)) {
  pbB <- fread(PB_B_FILE)
  if (!"annot" %in% names(pbB)) setnames(pbB, "subclass", "annot")   # 10_3_deg.R names it subclass
  pbB <- pbB[model == "B" & contrast == "dietHFD:sexM",
                          .(annot, gene, pbB_int_logFC = logFC, pbB_int_p = PValue, pbB_int_FDR = FDR)]
} else msg("no ", PB_B_FILE, "; edgeR interaction columns left out")

## ---- per cell type -----------------------------------------------------------------
summary <- list()
for (ct in setdiff(sort(unique(obj$annot)), DROP_TYPES)) {
  cells <- colnames(obj)[obj$annot == ct]
  msg(ct, ": ", length(cells), " nuclei")

  ## step 1: pooled
  r <- wilcox_diet(cells, FEATURES, MIN_PCT, LOGFC_MIN)
  if (is.null(r)) {
    msg("  skipped: < ", MIN_CELLS, " nuclei in a diet")
    summary[[ct]] <- data.table(annot = ct, n_nuclei = length(cells), tested = FALSE)
    next
  }
  n <- attr(r, "n")
  setnames(r, "p_val_adj", "p_bonferroni")   # Seurat's: p x nrow(obj), all genes in the object
  r[, FDR := p.adjust(p_val, "BH")]          # BH over the genes FindMarkers returned
  pb <- read_pb(ct)
  if (!is.null(pb)) r <- merge(r, pb, by = "gene", all.x = TRUE) else msg("  no pseudobulk table for ", ct)
  setorder(r, p_val)
  fwrite(r, file.path(OUTDIR, paste0("DEG_percell_HFDvsSD_", safe(ct), ".csv")))

  hits <- r[p_bonferroni < BONF_CUT & abs(avg_log2FC) > LFC_CUT, gene]
  msg("  pooled: ", nrow(r), " genes tested, ", length(hits), " at p_bonferroni < ", BONF_CUT,
      " and |log2FC| > ", LFC_CUT)
  s_row <- data.table(annot = ct, n_nuclei = length(cells), tested = TRUE,
                      nuc_SD = n[["SD"]], nuc_HFD = n[["HFD"]], n_genes = nrow(r),
                      n_hits = length(hits),
                      n_hits_up = r[gene %in% hits & avg_log2FC > 0, .N],
                      n_hits_pb_FDR05 = if (!is.null(pb)) r[gene %in% hits & pb_FDR < 0.05, .N] else NA_integer_)
  if (!length(hits)) { summary[[ct]] <- s_row; next }

  ## step 2: within each sex, hits only, no pct / logFC prefilter so every hit is returned
  hit_tab <- r[gene %in% hits]
  for (sx in c("F", "M")) {
    rs <- wilcox_diet(cells[obj$sex[match(cells, colnames(obj))] == sx], hits, 0, 0)
    if (is.null(rs)) { msg("  ", sx, ": skipped, < ", MIN_CELLS, " nuclei in a diet"); next }
    rs[, `:=`(p_bonferroni = p.adjust(p_val, "bonferroni"),   # over this cell type's hits
              FDR          = p.adjust(p_val, "BH"))]
    rs <- rs[, .(gene, avg_log2FC, pct.1, pct.2, p_val, p_bonferroni, FDR)]
    setnames(rs, names(rs)[-1], paste0(sx, "_", names(rs)[-1]))
    hit_tab <- merge(hit_tab, rs, by = "gene", all.x = TRUE)
  }

  ## per-mouse means of the hits, for the plots
  mm <- mouse_means(cells, hits)
  if (!is.null(pbB)) {
    hit_tab <- merge(hit_tab, pbB[annot == ct, !"annot"], by = "gene", all.x = TRUE)
    # Interaction BH re-applied over this cell type's hits only (those edgeR tested).
    # Valid because, in the balanced 3-per-group design, the pooled diet contrast used
    # for selection is orthogonal to the interaction contrast. pbB_int_FDR (over all
    # filterByExpr genes) is kept alongside.
    hit_tab[!is.na(pbB_int_p), pbB_int_FDR_hits := p.adjust(pbB_int_p, "BH")]
  }
  setorder(hit_tab, p_val)
  fwrite(hit_tab, file.path(OUTDIR, paste0("DEG_percell_sexsplit_", safe(ct), ".csv")))
  fwrite(mm[order(gene, diet, sex, sample)],
         file.path(OUTDIR, paste0("mouse_means_hits_", safe(ct), ".csv")))

  s_row[, `:=`(n_hits_F_bonf = if ("F_p_bonferroni" %in% names(hit_tab)) hit_tab[F_p_bonferroni < BONF_CUT & abs(F_avg_log2FC) > LFC_CUT, .N] else NA_integer_,
               n_hits_M_bonf = if ("M_p_bonferroni" %in% names(hit_tab)) hit_tab[M_p_bonferroni < BONF_CUT & abs(M_avg_log2FC) > LFC_CUT, .N] else NA_integer_)]
  summary[[ct]] <- s_row

}

summary <- rbindlist(summary, fill = TRUE)
fwrite(summary, file.path(OUTDIR, "DEG_percell_summary.csv"))
print(summary)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done")
