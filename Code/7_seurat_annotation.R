#!/usr/bin/env Rscript

# Provisional cluster annotation for the insular-cortex snRNA-seq object.
# Labels are read from Results/annotation_markers/cluster_annotation.csv, which
# is hand-curated from 6_seurat_markers output and is meant to be edited: rerun
# this script after changing it rather than editing labels here.
#
# Adds two metadata columns:
#   cluster_label  fine label, one per seurat_clusters level
#   cell_class     merged/coarse class (several clusters collapse into one)
# and draws the marker dotplots.

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_umap.rds")
ANNO    <- file.path(PROJECT, "Results", "annotation_markers", "cluster_annotation.csv")
PANEL   <- file.path(PROJECT, "Results", "annotation_markers", "markers_insular_cortex.csv")
OUTDIR  <- file.path(PROJECT, "Results", "7_seurat_annotation")
OBJOUT  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")

DOT_COL_LO <- "#eef4fc"
DOT_COL_HI <- "#0d366b"
SEED       <- 42

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

obj <- readRDS(IN)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "seurat_clusters"
msg("loaded ", ncol(obj), " nuclei | ", nlevels(obj$seurat_clusters), " clusters")

anno <- fread(ANNO, colClasses = list(character = "cluster"))
stopifnot(!anyDuplicated(anno$cluster),
          setequal(anno$cluster, levels(obj$seurat_clusters)))

m <- match(as.character(obj$seurat_clusters), anno$cluster)
obj$cluster_label <- anno$label[m]
obj$cell_class    <- anno$cell_class[m]
obj$confidence    <- anno$confidence[m]

# Order clusters by class then label so the dotplot reads as blocks.
CLASS_ORDER <- c("Excitatory neuron", "Claustrum/endopiriform", "Striatal",
                 "Inhibitory neuron", "Astrocyte", "Oligodendrocyte", "OPC",
                 "Microglia", "Perivascular macrophage", "Endothelial",
                 "Pericyte", "Smooth muscle", "VLMC/fibroblast", "Meningeal",
                 "QC-flagged")
stopifnot(setequal(CLASS_ORDER, unique(anno$cell_class)))
anno[, cell_class := factor(cell_class, levels = CLASS_ORDER)]
anno[, cluster_int := as.integer(cluster)]   # setorder takes column names, not expressions
setorder(anno, cell_class, label, cluster_int)
cl_order <- anno$cluster

obj$seurat_clusters <- factor(as.character(obj$seurat_clusters), levels = cl_order)
obj$cluster_named   <- factor(paste0(obj$seurat_clusters, " ", obj$cluster_label),
                              levels = paste0(anno$cluster, " ", anno$label))
obj$cell_class      <- factor(obj$cell_class, levels = CLASS_ORDER)

fwrite(anno, file.path(OUTDIR, "cluster_annotation_ordered.csv"))
msg("classes: ", paste(sprintf("%s=%d", levels(anno$cell_class),
                               table(anno$cell_class)), collapse = " "))
### Temporary:
markersAll <- fread('/scratch/user/uqschauq/Insular/Results/6_seurat_umap_markers/markers_all.csv') %>% filter(pct.1 >0.15, pct.2 > 0.15) %>% group_by(seurat_clusters) %>%
       slice_min(order_by = pb_FDR, n= 10)
x <- DoHeatmap(obj, features = markersAll$gene) + NoLegend()
ggsave(x, filename ='/scratch/user/uqschauq/Insular/Results/temp/plotTemp.jpeg', width = 14, height = 14, dpi = 300)
###

# ---- marker panel, ordered by the class each gene marks ---------------------
pan <- fread(PANEL)
pan <- pan[gene %in% rownames(obj[["RNA"]])]
GENE_CLASS_ORDER <- c("Neuronal", "Excitatory neuron", "Inhibitory neuron",
                      "Glia", "Immune", "Vascular", "Ependymal", "Activity state")
pan[, cell_class := factor(cell_class, levels = GENE_CLASS_ORDER)]
setorder(pan, cell_class, cell_type, gene)
genes <- unique(pan$gene)
msg("panel: ", length(genes), " genes of interest present in the object")

dot <- function(group, file, w, h) {
  Idents(obj) <- group
  p <- DotPlot(obj, features = genes, cols = c(DOT_COL_LO, DOT_COL_HI),
               dot.scale = 4) +
    scale_y_discrete(limits = rev(levels(obj[[group]][, 1]))) +
    labs(x = NULL, y = NULL,
         title = "Canonical marker panel across provisional clusters",
         subtitle = "colour = scaled mean expression; size = % of nuclei expressing") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 7),
          panel.grid  = element_line(colour = "grey93"),
          plot.subtitle = element_text(colour = "grey35"))
  ggsave(file.path(OUTDIR, file), p, width = w, height = h, dpi = 150,
         limitsize = FALSE)
}

dot("cluster_named", "dotplot_markers_by_cluster.png", 20, 13)
dot("cell_class",    "dotplot_markers_by_class.png",   20, 6)
msg("dotplots written")

# ---- UMAPs -----------------------------------------------------------------
ggsave(file.path(OUTDIR, "umap_cluster_label.png"),
       DimPlot(obj, reduction = "umap", group.by = "cluster_label", label = TRUE,
               label.size = 2.6, repel = TRUE, raster = FALSE, pt.size = 0.1) +
         coord_equal() + NoLegend() + labs(title = "Provisional cluster labels"),
       width = 11, height = 10, dpi = 150)

ggsave(file.path(OUTDIR, "umap_cell_class.png"),
       DimPlot(obj, reduction = "umap", group.by = "cell_class", label = TRUE,
               label.size = 3, repel = TRUE, raster = FALSE, pt.size = 0.1) +
         coord_equal() + labs(title = "Merged cell classes"),
       width = 12, height = 10, dpi = 150)

ggsave(file.path(OUTDIR, "umap_confidence.png"),
       DimPlot(obj, reduction = "umap", group.by = "confidence", raster = FALSE,
               pt.size = 0.1, order = c("low", "medium", "high")) +
         coord_equal() + labs(title = "Annotation confidence"),
       width = 9, height = 8, dpi = 150)

# ---- composition tables ----------------------------------------------------
md <- as.data.table(obj@meta.data)
fwrite(md[, .N, by = .(cell_class, cluster_label, seurat_clusters)][
         order(cell_class, -N)], file.path(OUTDIR, "label_sizes.csv"))
fwrite(dcast(md[, .N, by = .(cell_class, mouse)], cell_class ~ mouse,
             value.var = "N", fill = 0), file.path(OUTDIR, "class_by_mouse.csv"))
fwrite(dcast(md[, .N, by = .(cell_class, diet)], cell_class ~ diet,
             value.var = "N", fill = 0), file.path(OUTDIR, "class_by_diet.csv"))

print(md[, .N, by = cell_class][order(-N)])

saveRDS(obj, OBJOUT, compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done")
