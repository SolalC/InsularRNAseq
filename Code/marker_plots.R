#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(data.table)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_umap.rds")
MARKERS <- file.path(PROJECT, "Results", "6_seurat_umap_markers", "markers_all.csv")
OUTDIR  <- file.path(PROJECT, "Results", "6_seurat_umap_markers")

args <- commandArgs(trailingOnly = TRUE)
N_GENES          <- if (length(args) >= 1) as.integer(args[1]) else 10
CELLS_PER_CLUST  <- if (length(args) >= 2) as.integer(args[2]) else 200
LOGFC_MIN        <- 1
CONCORDANT_FIRST <- TRUE
SEED             <- 42

set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
stopifnot(file.exists(IN), file.exists(MARKERS))

mk <- fread(MARKERS)
msg("markers: ", nrow(mk), " rows, ", uniqueN(mk$seurat_clusters), " clusters")

if (!"concordant" %in% names(mk)) {
  msg("WARNING: no `concordant` column, ranking on avg_log2FC alone")
  CONCORDANT_FIRST <- FALSE
}

sel <- mk[avg_log2FC > LOGFC_MIN]
msg("avg_log2FC > ", LOGFC_MIN, ": ", nrow(sel), " rows, ",
    uniqueN(sel$seurat_clusters), " clusters")

empty <- setdiff(unique(mk$seurat_clusters), unique(sel$seurat_clusters))
if (length(empty)) {
  msg("no gene above the floor in cluster(s) ", paste(empty, collapse = ", "),
      " -- taking their top rows unfiltered")
  sel <- rbind(sel, mk[seurat_clusters %in% empty])
}

if (CONCORDANT_FIRST) {
  sel[, conc := !is.na(concordant) & concordant]
  setorder(sel, seurat_clusters, -conc, -avg_log2FC)
} else {
  setorder(sel, seurat_clusters, -avg_log2FC)
}

top <- sel[, head(.SD, N_GENES), by = seurat_clusters]
top[, seurat_clusters := factor(
  as.character(seurat_clusters),
  levels = sort(unique(as.integer(as.character(seurat_clusters)))))]
setorder(top, seurat_clusters)
if (CONCORDANT_FIRST)
  msg("selected ", nrow(top), " pairs, ", sum(top$conc), " concordant (",
      round(100 * mean(top$conc), 1), "%)")

msg("reading ", basename(IN))
obj <- readRDS(IN)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "seurat_clusters"
msg("loaded ", ncol(obj), " nuclei | ", nlevels(obj$seurat_clusters), " clusters")

absent <- setdiff(unique(top$gene), rownames(obj[["RNA"]]))
if (length(absent)) {
  msg("WARNING: dropping ", length(absent), " gene(s) not in the RNA assay: ",
      paste(head(absent, 10), collapse = ", "))
  top <- top[!gene %in% absent]
}
fwrite(top, file.path(OUTDIR, sprintf("marker_plot_top%d_genes.csv", N_GENES)))

genes <- top$gene
ug    <- unique(genes)
msg(length(genes), " heatmap rows (", length(ug), " distinct genes)")

msg("downsampling to <=", CELLS_PER_CLUST, " cells per cluster")
sub <- subset(obj, downsample = CELLS_PER_CLUST)
msg("subset: ", ncol(sub), " cells")

msg("scaling ", length(ug), " genes on the subset")
sub <- ScaleData(sub, features = ug, verbose = FALSE)
drawn <- sum(ug %in% rownames(GetAssayData(sub, assay = "RNA", layer = "scale.data")))
msg("in scale.data: ", drawn, " / ", length(ug))
stopifnot(drawn > 0)

hm <- DoHeatmap(sub, features = genes, size = 3, angle = 90, raster = TRUE) +
  NoLegend() +
  theme(axis.text.y = element_text(size = 4))

h <- min(0.085 * length(genes) + 3, 200)
w <- min(0.010 * ncol(sub) + 8, 60)
msg("heatmap: ", round(w, 1), " x ", round(h, 1), " in")
ggsave(file.path(OUTDIR, sprintf("heatmap_top%d_markers.png", N_GENES)), hm,
       width = w, height = h, dpi = 150, limitsize = FALSE)
ggsave(file.path(OUTDIR, sprintf("heatmap_top%d_markers.pdf", N_GENES)), hm,
       width = w, height = h, limitsize = FALSE)

dp <- DotPlot(obj, features = ug) + RotatedAxis() +
  theme(axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 6)) +
  labs(x = NULL, y = "cluster")

dw <- min(0.11 * length(ug) + 4, 200)
dh <- min(0.16 * nlevels(obj$seurat_clusters) + 3, 60)
msg("dotplot: ", round(dw, 1), " x ", round(dh, 1), " in")
ggsave(file.path(OUTDIR, sprintf("dotplot_top%d_markers.png", N_GENES)), dp,
       width = dw, height = dh, dpi = 150, limitsize = FALSE)
ggsave(file.path(OUTDIR, sprintf("dotplot_top%d_markers.pdf", N_GENES)), dp,
       width = dw, height = dh, limitsize = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "sessionInfo_marker_plots.txt"))
msg("done: ", OUTDIR)
