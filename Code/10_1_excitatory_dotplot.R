# Dotplot of the top marker genes of the EXCITATORY subclustering at the
# 16-cluster solution (resolution 0.2 of the step-10 sweep).
#
# WHY THIS IS A SEPARATE SCRIPT. Step 10 ran a resolution sweep and then did all
# of its marker work at the PRIMARY resolution (0.6, 28 subclusters), so
# Results/10_neuronal_subclustering/excitatory/markers_all.csv describes 28
# clusters and says nothing about the 16-cluster partition the user chose. The
# cluster labels are not nested in a way that lets the 28-cluster markers be
# collapsed, so the markers are recomputed from scratch on res 0.2 with the SAME
# two-layer test (Wilcoxon + mouse-paired pseudobulk edgeR) via
# lib_cluster_markers.R. Nothing about the embedding is recomputed: the object,
# the PCA, the UMAP and every RNA_snn_res.* column are read back from the step-10
# object exactly as saved.
#
# GENE SELECTION. N_GENES per cluster, ranked concordant-first then by
# avg_log2FC -- the same ranking step 10 used for its single-marker dotplot. A
# gene is only used once: if a cluster's next-best gene was already taken by an
# earlier cluster, the following candidate is used instead, so the panel has
# 16 x N_GENES distinct genes and the diagonal is readable. Genes that are NOT
# concordant (i.e. not supported by the animal-level test) are still allowed
# when a cluster has too few concordant ones, and are flagged in the CSV --
# never silently mixed in.
#
# Usage:
#   sbatch 10_1_excitatory_dotplot.sh            # res 0.2, 5 genes per cluster
#   sbatch 10_1_excitatory_dotplot.sh 0.2 5

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(data.table)
})

source("/scratch/user/uqschauq/Insular/Code/lib_cluster_markers.R")

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_sub_excitatory.rds")
OUTDIR  <- file.path(PROJECT, "Results", "10_1_excitatory_dotplot")

SEQ_LO <- "#eef4fc"
SEQ_HI <- "#0d366b"

args    <- commandArgs(trailingOnly = TRUE)
RES     <- if (length(args) >= 1 && nzchar(args[1])) as.numeric(args[1]) else 0.2
N_GENES <- if (length(args) >= 2 && nzchar(args[2])) as.integer(args[2]) else 5L
if (!is.finite(RES))     stop("first argument must be a resolution, e.g. 0.2")
if (!is.finite(N_GENES) || N_GENES < 1L) stop("second argument must be a positive integer")
EXPECT_K <- 16L   # what res 0.2 gave in the step-10 sweep; checked, not enforced

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
msg("resolution ", RES, " | ", N_GENES, " genes per cluster")

stopifnot(file.exists(IN))
msg("reading ", basename(IN), " (", round(file.size(IN) / 2^30, 1), " GB)")
sub <- readRDS(IN)
msg("loaded ", format(ncol(sub), big.mark = ","), " excitatory nuclei")

res_col <- paste0("RNA_snn_res.", RES)
if (!res_col %in% names(sub@meta.data))
  stop("no '", res_col, "' column -- resolutions present: ",
       paste(grep("^RNA_snn_res\\.", names(sub@meta.data), value = TRUE), collapse = ", "))

# Order the levels numerically (0,1,2,...,15) rather than as the strings Seurat
# leaves behind, so the dotplot rows run in cluster order.
cl <- sub@meta.data[[res_col]]
cl <- factor(as.character(cl), levels = as.character(sort(as.integer(unique(as.character(cl))))))
sub$dot_cluster <- cl
Idents(sub) <- "dot_cluster"
K <- nlevels(cl)
msg(res_col, ": ", K, " clusters")
if (K != EXPECT_K)
  msg("NOTE: expected ", EXPECT_K, " clusters at this resolution, got ", K,
      " -- proceeding with what is in the object")

sizes <- as.data.table(sub@meta.data)[, .N, by = dot_cluster][order(dot_cluster)]
fwrite(sizes, file.path(OUTDIR, sprintf("cluster_sizes_res%s.csv", RES)))
print(sizes)

## ---- markers: same two-layer test as step 10, on the 16-cluster partition ----

if (length(Layers(sub[["RNA"]], search = "counts")) > 1) sub <- JoinLayers(sub)
mk <- cluster_markers(sub, clusters = sub$dot_cluster, mouse = sub$mouse, log = msg)
setnames(mk, "cluster_id", "cluster")
fwrite(mk, file.path(OUTDIR, sprintf("markers_all_res%s.csv", RES)))
top20 <- mk[order(cluster, -concordant, -avg_log2FC)][, head(.SD, 20), by = cluster]
fwrite(top20, file.path(OUTDIR, sprintf("markers_top20_res%s.csv", RES)))
msg("markers: ", nrow(mk), " rows, ", sum(mk$concordant, na.rm = TRUE), " concordant")

conc <- mk[, .(n_wilcox = .N, n_concordant = sum(concordant, na.rm = TRUE)), by = cluster][
  order(cluster)]
conc[, pct_concordant := round(100 * n_concordant / n_wilcox, 1)]
fwrite(conc, file.path(OUTDIR, sprintf("concordance_by_cluster_res%s.csv", RES)))
print(conc)

## ---- pick N_GENES per cluster, no gene used twice ---------------------------

cand <- mk[order(cluster, -concordant, -avg_log2FC)]
used <- character(0)
picked <- rbindlist(lapply(levels(cl), function(k) {
  ck <- cand[cluster == k & !gene %in% used]
  take <- head(ck, N_GENES)
  used <<- c(used, take$gene)
  if (nrow(take) < N_GENES)
    msg("  cluster ", k, ": only ", nrow(take), " unused marker(s) available")
  if (nrow(take) && sum(!take$concordant, na.rm = TRUE))
    msg("  cluster ", k, ": ", sum(!take$concordant, na.rm = TRUE), " of ", nrow(take),
        " gene(s) NOT concordant with the animal-level test")
  take
}), fill = TRUE)
picked[, rank_in_cluster := seq_len(.N), by = cluster]
fwrite(picked, file.path(OUTDIR, sprintf("dotplot_genes_res%s.csv", RES)))
msg("selected ", nrow(picked), " genes across ", uniqueN(picked$cluster), " clusters")

genes <- picked$gene
stopifnot(!anyDuplicated(genes))

## ---- the dotplot ------------------------------------------------------------

# Rows = clusters in numeric order, top to bottom; columns = the selected genes
# in cluster order, so support runs down the diagonal.
p <- DotPlot(sub, features = genes, group.by = "dot_cluster",
             dot.scale = 4.5, scale = TRUE) +
  scale_colour_gradient(low = SEQ_LO, high = SEQ_HI, name = "scaled\nexpression") +
  scale_y_discrete(limits = rev(levels(cl))) +
  labs(x = NULL, y = "subcluster",
       title = sprintf("Excitatory neurons, %d subclusters (res %s)", K, RES),
       subtitle = sprintf("top %d markers per subcluster, concordant-first then avg_log2FC",
                          N_GENES)) +
  guides(size = guide_legend(title = "% nuclei\nexpressing")) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5,
                                   face = "italic"),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_line(linewidth = 0.2, colour = "grey92"),
        panel.grid.minor = element_blank())

W <- 0.16 * length(genes) + 3.2
H <- 0.30 * K + 3.0
ggsave(file.path(OUTDIR, sprintf("dotplot_top%d_res%s.png", N_GENES, RES)), p,
       width = W, height = H, dpi = 300, limitsize = FALSE)
ggsave(file.path(OUTDIR, sprintf("dotplot_top%d_res%s.pdf", N_GENES, RES)), p,
       width = W, height = H, limitsize = FALSE)
msg("wrote dotplot (", round(W, 1), " x ", round(H, 1), " in)")

# Same panel transposed -- genes down the page. Easier to read gene names when
# the figure is dropped into a portrait slide.
pt <- p + coord_flip() +
  scale_y_discrete(limits = levels(cl)) +
  labs(x = "subcluster", y = NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 8,
                                   face = "plain"),
        axis.text.y = element_text(size = 6.5, face = "italic"))
ggsave(file.path(OUTDIR, sprintf("dotplot_top%d_res%s_flipped.png", N_GENES, RES)), pt,
       width = 0.30 * K + 4.5, height = 0.15 * length(genes) + 2.5, dpi = 300,
       limitsize = FALSE)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)

## Check:
# Cluster 0:
p0 <- UMAPPlot(sub, raster = F, label = T)
gene <- c('Rorb', 'Rsoi1m', 'Zmat4', 'Pld5', 'Grik1', 'Slc17a6', 'Cux2', 'Slc17a7', 'Fezf2')
FeaturePlot(sub, features = gene, raster = F)
