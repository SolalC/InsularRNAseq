#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 6_seurat_umap_markers.R
#
# Embedding, clustering and cluster markers in one job. Merged from the former
# 5_seurat_umap.R and 6_seurat_markers.R on 2026-08-31, when the pipeline was
# reordered so that ambient removal and doublet calling (step 4) run before QC
# (step 5). Everything from normalisation onward is a function of the surviving
# nuclei, so it is all rebuilt here from HFDseq_v2_filtered.rds.
#
# The two halves share one 24 GB object in memory rather than round-tripping it
# through disk. The object IS still written between them, so a marker failure
# does not cost the embedding.
#
# Markers come in two layers written to one table:
#   layer 1  Seurat one-vs-rest Wilcoxon (presto-accelerated). Gives
#            avg_log2FC and the pct.1/pct.2 detection rates annotation keys on.
#   layer 2  Pseudobulk edgeR, cluster-vs-rest paired within mouse. Wilcoxon
#            treats every nucleus as an independent replicate, so its p-values
#            are uninterpretable at n ~ 20k and a "marker" can be driven by one
#            or two animals. Layer 2 re-tests each candidate with the mouse as
#            the unit of replication and reports how many animals support it.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(data.table)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_filtered.rds")
OUTDIR  <- file.path(PROJECT, "Results", "6_seurat_umap_markers")
OBJOUT  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_umap.rds")

# embedding
TRANSGENES <- c("EGFP", "tdTomato", "EBFP2", "AAV-shared")
N_HVG      <- 2000
N_PCS      <- 50
DIMS       <- 1:30
RESOLUTION <- 0.5
SEED       <- 42

MOUSE_BG <- "grey88"
SEQ_LO   <- "#eef4fc"
SEQ_HI   <- "#0d366b"

# markers, layer 1
MIN_PCT     <- 0.25
LOGFC_THR   <- 0.25
PREFILT_PCT <- 0.10   # keep a gene detected in >= this fraction of >=1 cluster
# markers, layer 2
MIN_CELLS   <- 10     # drop a (cluster, mouse) pseudobulk unit below this
MIN_MICE    <- 4      # clusters with fewer usable animals get no pseudobulk stats
TOP_N       <- 20

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

obj <- readRDS(IN)
msg("loaded ", ncol(obj), " nuclei x ", nrow(obj), " features")

counts <- LayerData(obj, assay = "RNA", layer = "counts")
tg <- intersect(TRANSGENES, rownames(counts))
mg <- setdiff(rownames(counts), tg)
stopifnot(length(tg) == length(TRANSGENES))

md  <- obj@meta.data
obj <- CreateSeuratObject(counts = counts[mg, ], meta.data = md, project = "HFDseq_v2")
obj[["TG"]] <- CreateAssay5Object(counts = counts[tg, , drop = FALSE])
DefaultAssay(obj) <- "RNA"
rm(counts, md); invisible(gc())
msg("RNA ", nrow(obj[["RNA"]]), " genes | TG ", nrow(obj[["TG"]]), " features")

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = N_HVG)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = N_PCS, verbose = FALSE)
msg("PCA done")

ggsave(file.path(OUTDIR, "elbow.png"),
       ElbowPlot(obj, ndims = N_PCS), width = 7, height = 5, dpi = 150)

obj <- FindNeighbors(obj, dims = DIMS)
obj <- FindClusters(obj, resolution = RESOLUTION)
msg("clusters: ", length(levels(obj$seurat_clusters)))

obj <- RunUMAP(obj, dims = DIMS, seed.use = SEED)
msg("UMAP done")

um <- function(group, file, label = FALSE, legend = TRUE) {
  p <- DimPlot(obj, reduction = "umap", group.by = group, label = label,
               label.size = 3, raster = FALSE, pt.size = 0.1, shuffle = TRUE) +
    coord_equal()
  if (!legend) p <- p + NoLegend()
  ggsave(file.path(OUTDIR, file), p, width = 9, height = 8, dpi = 150)
}

um("seurat_clusters", "umap_clusters.png", label = TRUE, legend = FALSE)
um("sample",          "umap_sample.png")
um("diet",            "umap_diet.png")
um("sex",             "umap_sex.png")

ggsave(file.path(OUTDIR, "umap_split_diet.png"),
       DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
               split.by = "diet", label = TRUE, label.size = 3,
               raster = FALSE, pt.size = 0.1) + coord_equal() + NoLegend(),
       width = 14, height = 7, dpi = 150)

obj$mouse <- factor(obj$mouse,
                    levels = unique(obj$mouse[order(as.integer(obj$sample))]))
msg("mice: ", paste(levels(obj$mouse), collapse = " "))

emb <- as.data.table(Embeddings(obj, "umap"))
setnames(emb, 1:2, c("UMAP_1", "UMAP_2"))
emb[, mouse := obj$mouse]
bg <- emb[, .(UMAP_1, UMAP_2)]

emb[, cluster := obj$seurat_clusters]
cl_pal <- scales::hue_pal()(nlevels(obj$seurat_clusters))
names(cl_pal) <- levels(obj$seurat_clusters)

ggsave(file.path(OUTDIR, "umap_mouse_facets.png"),
       ggplot(emb, aes(UMAP_1, UMAP_2)) +
         geom_point(data = bg, colour = MOUSE_BG, size = 0.3, stroke = 0) +
         geom_point(aes(colour = cluster), size = 0.35, stroke = 0) +
         facet_wrap(~ mouse, ncol = 4) +
         scale_colour_manual(values = cl_pal, name = "cluster") +
         guides(colour = guide_legend(override.aes = list(size = 2.5), ncol = 1)) +
         coord_equal() +
         labs(title = "Clusters by mouse of origin",
              subtitle = sprintf("each panel: that animal coloured by cluster, all %s nuclei in grey behind",
                                 format(ncol(obj), big.mark = ","))) +
         theme_classic(base_size = 10) +
         theme(strip.background = element_blank(),
               strip.text       = element_text(face = "bold"),
               axis.text        = element_text(colour = "grey40"),
               legend.key.height = unit(0.8, "lines"),
               plot.subtitle    = element_text(colour = "grey35")),
       width = 14.5, height = 10.5, dpi = 150)

cm     <- as.data.table(obj@meta.data)[, .N, by = .(seurat_clusters, mouse)]
wide_n <- dcast(cm, seurat_clusters ~ mouse, value.var = "N", fill = 0, drop = FALSE)
long_n <- melt(wide_n, id.vars = "seurat_clusters",
               variable.name = "mouse", value.name = "n")
long_n[, pct := 100 * n / sum(n), by = seurat_clusters]

even     <- 100 / nlevels(obj$mouse)
LABEL_AT <- round(2 * even)

ggsave(file.path(OUTDIR, "umap_mouse_cluster_heatmap.png"),
       ggplot(long_n, aes(mouse, seurat_clusters, fill = pct)) +
         geom_tile(colour = "white", linewidth = 0.5) +
         geom_text(aes(label = fifelse(pct >= LABEL_AT, sprintf("%.0f", pct), ""),
                       colour = pct > 55), size = 2.6) +
         scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#3d3d3a"),
                             guide = "none") +
         scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, limits = c(0, 100),
                             breaks = c(0, 25, 50, 75, 100), name = "% of\ncluster") +
         scale_y_discrete(limits = rev(levels(long_n$seurat_clusters))) +
         labs(x = NULL, y = "cluster",
              title = "Mouse contribution to each cluster",
              subtitle = sprintf("cells labelled above %.0f%%; even contribution would be %.1f%% per mouse",
                                 LABEL_AT, even)) +
         theme_minimal(base_size = 10) +
         theme(panel.grid     = element_blank(),
               axis.text.x    = element_text(angle = 45, hjust = 1),
               plot.subtitle  = element_text(colour = "grey35")),
       width = 9, height = 0.28 * nlevels(obj$seurat_clusters) + 2.5, dpi = 150)

long_n[, pct_animal := 100 * n / sum(n), by = mouse]

ggsave(file.path(OUTDIR, "umap_mouse_cluster_bars.png"),
       ggplot(long_n, aes(mouse, pct_animal, fill = seurat_clusters)) +
         geom_col(width = 0.8, colour = "white", linewidth = 0.3) +
         geom_text(aes(label = fifelse(pct_animal >= 4,
                                       as.character(seurat_clusters), "")),
                   position = position_stack(vjust = 0.5),
                   size = 2.4, colour = "grey15") +
         scale_fill_manual(values = cl_pal, name = "cluster") +
         scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
         guides(fill = guide_legend(ncol = 1)) +
         labs(x = NULL, y = "% of the animal's nuclei",
              title = "Cluster composition of each mouse",
              subtitle = "each bar sums to 100%; segments labelled above 4%") +
         theme_classic(base_size = 11) +
         theme(axis.text.x       = element_text(angle = 45, hjust = 1),
               legend.key.height = unit(0.8, "lines"),
               plot.subtitle     = element_text(colour = "grey35")),
       width = 11, height = 7, dpi = 150)

fwrite(dcast(long_n, seurat_clusters ~ mouse, value.var = "pct_animal"),
       file.path(OUTDIR, "mouse_cluster_percent.csv"))
fwrite(wide_n, file.path(OUTDIR, "cluster_mouse_counts.csv"))
fwrite(dcast(long_n, seurat_clusters ~ mouse, value.var = "pct"),
       file.path(OUTDIR, "cluster_mouse_percent.csv"))

dom <- long_n[, .(top_mouse = mouse[which.max(pct)], top_pct = round(max(pct), 1)),
              by = seurat_clusters][top_pct >= 50]
if (nrow(dom)) {
  msg("clusters >=50% from a single mouse:")
  print(dom)
} else {
  msg("no cluster is >=50% from a single mouse")
}

ggsave(file.path(OUTDIR, "umap_qc.png"),
       FeaturePlot(obj, reduction = "umap",
                   features = c("nFeature_mouse", "nCount_mouse", "percent_mt"),
                   raster = FALSE, pt.size = 0.1, ncol = 3),
       width = 16, height = 5.5, dpi = 150)

comp <- as.data.table(obj@meta.data)[, .N, by = .(seurat_clusters, sample, diet, sex)]
fwrite(comp[order(seurat_clusters, sample)],
       file.path(OUTDIR, "cluster_composition.csv"))
fwrite(as.data.table(table(obj$seurat_clusters)),
       file.path(OUTDIR, "cluster_sizes.csv"))
fwrite(data.table(feature = VariableFeatures(obj)),
       file.path(OUTDIR, "variable_features.csv"))

saveRDS(obj, OBJOUT, compress = FALSE)
msg("saved ", OBJOUT, " (", round(file.size(OBJOUT) / 2^30, 2), " GB)")

## ---- MARKERS --------------------------------------------------------------

DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]], search = "counts")) > 1) obj <- JoinLayers(obj)
Idents(obj) <- "seurat_clusters"

clust <- obj$seurat_clusters
mouse <- factor(obj$mouse)
msg("loaded ", ncol(obj), " nuclei x ", nrow(obj), " genes | ",
    nlevels(clust), " clusters | ", nlevels(mouse), " mice")

counts <- LayerData(obj, assay = "RNA", layer = "counts")

# ---- gene pre-filter -------------------------------------------------------
# Detection rate per gene per cluster, via a sparse indicator multiply.
# Binarise by overwriting the nonzero slot rather than `counts > 0`, which
# would materialise a second matrix of the same nnz.
bin <- counts
bin@x <- rep.int(1, length(bin@x))
ind_c <- sparseMatrix(i = seq_len(ncol(counts)), j = as.integer(clust),
                      x = 1, dims = c(ncol(counts), nlevels(clust)),
                      dimnames = list(NULL, levels(clust)))
det <- (bin %*% ind_c) %*% Diagonal(x = 1 / as.numeric(table(clust)))
max_det    <- apply(as.matrix(det), 1, max)
keep_genes <- rownames(counts)[max_det >= PREFILT_PCT]
rm(bin); invisible(gc())
msg("gene pre-filter: ", length(keep_genes), " / ", nrow(counts),
    " genes detected in >=", PREFILT_PCT * 100, "% of cells in at least one cluster")

# ---- layer 1: Seurat default Wilcoxon --------------------------------------
msg("layer 1: FindAllMarkers (wilcox, presto = ",
    requireNamespace("presto", quietly = TRUE), ")")
mk <- FindAllMarkers(obj, features = keep_genes, only.pos = TRUE,
                     min.pct = MIN_PCT, logfc.threshold = LOGFC_THR,
                     test.use = "wilcox", verbose = FALSE)
setDT(mk)
setnames(mk, "cluster", "seurat_clusters")
mk[, seurat_clusters := factor(as.character(seurat_clusters), levels = levels(clust))]
msg("layer 1: ", nrow(mk), " marker rows across ",
    uniqueN(mk$seurat_clusters), " clusters")

# ---- layer 2: pseudobulk, cluster vs rest, paired within mouse -------------
msg("layer 2: aggregating pseudobulk by (cluster x mouse)")
grp   <- interaction(clust, mouse, drop = FALSE, sep = "|")
ind_g <- sparseMatrix(i = seq_len(ncol(counts)), j = as.integer(grp),
                      x = 1, dims = c(ncol(counts), nlevels(grp)),
                      dimnames = list(NULL, levels(grp)))
pb    <- as.matrix(counts[keep_genes, ] %*% ind_g)          # genes x (cluster|mouse)
ncell <- as.integer(table(grp))
names(ncell) <- levels(grp)

key   <- tstrsplit(colnames(pb), "|", fixed = TRUE)
pb_cl <- factor(key[[1]], levels = levels(clust))
pb_ms <- factor(key[[2]], levels = levels(mouse))
# per-mouse totals, so "rest" is total minus the cluster being tested
tot_ms <- sapply(levels(mouse), function(m) rowSums(pb[, pb_ms == m, drop = FALSE]))
msg("layer 2: ", ncol(pb), " units, median ", median(ncell), " cells, ",
    sum(ncell < MIN_CELLS), " below ", MIN_CELLS, " cells")

pb_one <- function(k) {
  sel <- pb_cl == k & ncell >= MIN_CELLS
  mice_k <- as.character(pb_ms[sel])
  if (length(mice_k) < MIN_MICE) return(NULL)

  inmat  <- pb[, sel, drop = FALSE]
  restm  <- tot_ms[, mice_k, drop = FALSE] - inmat
  colnames(inmat) <- paste0("in_",   mice_k)
  colnames(restm) <- paste0("rest_", mice_k)

  mat <- cbind(inmat, restm)
  grp2 <- factor(rep(c("in", "rest"), each = length(mice_k)), levels = c("rest", "in"))
  ms2  <- factor(rep(mice_k, 2))

  y <- DGEList(counts = mat, group = grp2)
  design <- model.matrix(~ ms2 + grp2)
  y <- y[filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
  if (nrow(y) < 10) return(NULL)
  y <- normLibSizes(y)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  res <- topTags(glmQLFTest(fit, coef = ncol(design)), n = Inf, sort.by = "none")$table

  # per-mouse support: in how many animals is this cluster above that animal's rest
  cpm_all <- cpm(y, log = FALSE)
  sup <- rowSums(cpm_all[, grp2 == "in", drop = FALSE] >
                 cpm_all[, grp2 == "rest", drop = FALSE])

  data.table(seurat_clusters = k, gene = rownames(res),
             pb_log2FC = res$logFC, pb_PValue = res$PValue, pb_FDR = res$FDR,
             n_mice_supporting = as.integer(sup), n_mice_tested = length(mice_k))
}

pbres <- rbindlist(lapply(levels(clust), function(k) {
  r <- tryCatch(pb_one(k), error = function(e) {
    msg("  cluster ", k, ": pseudobulk failed -- ", conditionMessage(e)); NULL })
  if (!is.null(r)) msg("  cluster ", k, ": ", r$n_mice_tested[1], " mice, ",
                       nrow(r), " genes tested")
  r
}), fill = TRUE)

# Guard: if every cluster failed layer 2, keep the join well-typed so the
# Wilcoxon table still lands rather than dying at the merge.
if (nrow(pbres) == 0L) {
  pbres <- data.table(seurat_clusters = character(), gene = character(),
                      pb_log2FC = numeric(), pb_PValue = numeric(),
                      pb_FDR = numeric(),
                      n_mice_supporting = integer(), n_mice_tested = integer())
  msg("layer 2: WARNING -- no cluster produced pseudobulk stats")
}

skipped <- setdiff(levels(clust), unique(as.character(pbres$seurat_clusters)))
if (length(skipped))
  msg("layer 2: no pseudobulk stats for cluster(s) ", paste(skipped, collapse = ", "),
      " (< ", MIN_MICE, " usable animals) -- rely on layer 1 there")

# ---- join and write --------------------------------------------------------
mk[, seurat_clusters := as.character(seurat_clusters)]
pbres[, seurat_clusters := as.character(seurat_clusters)]
out <- merge(mk, pbres, by = c("seurat_clusters", "gene"), all.x = TRUE, sort = FALSE)

out[, concordant := !is.na(pb_FDR) & pb_FDR < 0.05 & pb_log2FC > 0 &
      n_mice_supporting >= pmax(MIN_MICE, ceiling(0.75 * n_mice_tested))]
out[, seurat_clusters := factor(seurat_clusters, levels = levels(clust))]
setorder(out, seurat_clusters, -avg_log2FC)
setcolorder(out, c("seurat_clusters", "gene", "avg_log2FC", "pct.1", "pct.2",
                   "p_val", "p_val_adj", "pb_log2FC", "pb_PValue", "pb_FDR",
                   "n_mice_supporting", "n_mice_tested", "concordant"))

fwrite(out, file.path(OUTDIR, "markers_all.csv"))
top <- out[order(seurat_clusters, -concordant, -avg_log2FC)][, head(.SD, TOP_N),
                                                             by = seurat_clusters]
fwrite(top, file.path(OUTDIR, "markers_top20.csv"))
msg("wrote ", nrow(out), " rows | ", sum(out$concordant, na.rm = TRUE),
    " concordant across both layers")

conc <- out[, .(n_wilcox = .N, n_concordant = sum(concordant, na.rm = TRUE)),
            by = seurat_clusters][order(seurat_clusters)]
conc[, pct_concordant := round(100 * n_concordant / n_wilcox, 1)]
fwrite(conc, file.path(OUTDIR, "concordance_by_cluster.csv"))
print(conc)

# top concordant marker per cluster, for a quick visual sanity check
best <- out[concordant == TRUE][order(seurat_clusters, -avg_log2FC)][
  , .(gene = gene[1]), by = seurat_clusters][!is.na(gene)]
if (nrow(best)) {
  ggsave(file.path(OUTDIR, "dotplot_top_marker.png"),
         DotPlot(obj, features = unique(best$gene)) +
           coord_flip() +
           labs(x = NULL, y = "cluster",
                title = "Top concordant marker per cluster") +
           theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                            size = 7),
                 axis.text.y = element_text(size = 7)),
         width = 11, height = 0.22 * length(unique(best$gene)) + 3,
         dpi = 150, limitsize = FALSE)
}


writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
