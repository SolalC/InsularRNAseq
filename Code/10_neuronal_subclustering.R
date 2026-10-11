#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10_neuronal_subclustering.R  --  run by 10_neuronal_subclustering.sh.
#
# The whole-dataset clustering resolves the neurons into two blocks, Excitatory
# and Inhibitory, and step 9 found essentially nothing between diets at that
# grain (user, 2026-09-01: "the grouping might be too high"). This takes each
# neuronal lineage on its own and re-runs the pipeline INSIDE it: new HVGs, new
# PCA, new neighbours, new clusters, new markers. The two lineages never see
# each other -- variable genes chosen on excitatory neurons would otherwise be
# dominated by the excitatory/inhibitory split itself, which is exactly the
# axis being conditioned away.
#
# Everything upstream is inherited: the nuclei are the QC-passing, singlet,
# CellBender-corrected set, and the lineage is the manual `annot` call.
#
# WHY A RESOLUTION SWEEP RATHER THAN ONE NUMBER. Subcluster count is the whole
# point here, so the script clusters at every resolution in RES_SWEEP, writes a
# table of what each one gives (cluster count, smallest cluster, how many
# clusters are single-animal), and draws a UMAP for each. PRIMARY_RES is only
# what the markers and the saved Idents use; it is a starting point, not a
# claim. Re-run with a different second argument to move it.
#
# THE SINGLE-ANIMAL CHECK IS THE ONE THAT MATTERS. Subclustering finds structure
# by construction: push the resolution and it will split anything. A subcluster
# drawn mostly from one mouse is a batch artefact, not a cell type, and at this
# grain that is the most likely failure. Per resolution the script reports the
# largest single-mouse share of any cluster, and at PRIMARY_RES it writes the
# full mouse-by-cluster table and names every cluster above SINGLE_MOUSE_WARN.
#
# Markers are the same two-layer test the parent clustering used (Wilcoxon plus
# pseudobulk paired within mouse) via lib_cluster_markers.R -- rank on
# `concordant`, not on the Wilcoxon p-value.
#
# Usage:
#   10_neuronal_subclustering.sh                    both lineages, res 0.6
#   10_neuronal_subclustering.sh excitatory         one lineage
#   10_neuronal_subclustering.sh inhibitory 0.8     one lineage, other resolution
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

source("/scratch/user/uqschauq/Insular/Code/lib_cluster_markers.R")

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
OUTROOT <- file.path(PROJECT, "Results", "10_neuronal_subclustering")
OBJDIR  <- file.path(PROJECT, "Data", "Seurat")

LINEAGES <- list(
  excitatory = "Excitatory Neurons",
  inhibitory = "Inhibitory Neurons")

# Same embedding parameters as 6_seurat_umap_markers.R, so a subcluster is not
# also a change of method.
N_HVG      <- 2000
N_PCS      <- 50
DIMS       <- 1:30
RES_SWEEP  <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
PRIMARY_RES_DEFAULT <- 0.6
MIN_CLUSTER_CELLS   <- 50    # reported, never merged automatically
SINGLE_MOUSE_WARN   <- 50    # % of a cluster from one animal before it is named
SEED       <- 42

SEQ_LO <- "#eef4fc"
SEQ_HI <- "#0d366b"

args <- commandArgs(trailingOnly = TRUE)
WHICH <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "both"
PRIMARY_RES <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(args[2]) else
                 PRIMARY_RES_DEFAULT
if (!WHICH %in% c("both", names(LINEAGES)))
  stop("first argument must be one of: both, ", paste(names(LINEAGES), collapse = ", "))
if (!is.finite(PRIMARY_RES)) stop("second argument must be a number (resolution)")
RES_SWEEP <- sort(unique(c(RES_SWEEP, PRIMARY_RES)))
todo <- if (WHICH == "both") names(LINEAGES) else WHICH

dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
msg("lineages: ", paste(todo, collapse = ", "), " | primary resolution ", PRIMARY_RES)

stopifnot(file.exists(IN))
msg("reading ", basename(IN))
full <- readRDS(IN)
if (!"annot" %in% names(full@meta.data))
  stop("no 'annot' column in the object -- the manual annotation was not saved")
full$annot_chr <- as.character(full$annot)
msg("loaded ", format(ncol(full), big.mark = ","), " nuclei")
print(sort(table(full$annot_chr), decreasing = TRUE))

res_col <- function(r) paste0("RNA_snn_res.", r)

for (lin in todo) {
  label  <- LINEAGES[[lin]]
  OUTDIR <- file.path(OUTROOT, lin)
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  msg("=== ", label, " -> ", OUTDIR)

  cells <- colnames(full)[full$annot_chr == label]
  if (!length(cells)) stop("no nuclei with annot == '", label, "'")
  msg(format(length(cells), big.mark = ","), " nuclei")

  sub <- subset(full, cells = cells)
  DefaultAssay(sub) <- "RNA"
  # The parent's clustering and HVGs are meaningless inside the subset; keep the
  # parent cluster id for traceability and drop everything else derived.
  sub$parent_cluster <- sub$seurat_clusters
  sub$parent_annot   <- sub$annot_chr

  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = N_HVG,
                              verbose = FALSE)
  sub <- ScaleData(sub, features = VariableFeatures(sub), verbose = FALSE)
  sub <- RunPCA(sub, features = VariableFeatures(sub), npcs = N_PCS, verbose = FALSE)
  msg("  PCA done on ", length(VariableFeatures(sub)), " HVGs")
  ggsave(file.path(OUTDIR, "elbow.png"), ElbowPlot(sub, ndims = N_PCS),
         width = 7, height = 5, dpi = 150)
  fwrite(data.table(feature = VariableFeatures(sub)),
         file.path(OUTDIR, "variable_features.csv"))

  sub <- FindNeighbors(sub, dims = DIMS, verbose = FALSE)
  for (r in RES_SWEEP) {
    sub <- FindClusters(sub, resolution = r, verbose = FALSE)
    msg("  resolution ", r, ": ", nlevels(sub@meta.data[[res_col(r)]]), " clusters")
  }
  sub <- RunUMAP(sub, dims = DIMS, seed.use = SEED, verbose = FALSE)
  msg("  UMAP done")

  # What each resolution buys, including the artefact check.
  mouse_f <- factor(sub$mouse)
  sweep <- rbindlist(lapply(RES_SWEEP, function(r) {
    cl <- droplevels(as.factor(sub@meta.data[[res_col(r)]]))
    tab <- table(cl, mouse_f)
    top_pct <- 100 * apply(tab, 1, max) / rowSums(tab)
    data.table(resolution = r, n_clusters = nlevels(cl),
               min_cluster = min(table(cl)), median_cluster = median(table(cl)),
               n_below_min = sum(table(cl) < MIN_CLUSTER_CELLS),
               max_single_mouse_pct = round(max(top_pct), 1),
               n_clusters_over_warn = sum(top_pct >= SINGLE_MOUSE_WARN))
  }))
  fwrite(sweep, file.path(OUTDIR, "resolution_sweep.csv"))
  msg("  resolution sweep:"); print(sweep)

  ggsave(file.path(OUTDIR, "umap_resolution_sweep.png"),
         wrap_plots(lapply(RES_SWEEP, function(r)
           DimPlot(sub, reduction = "umap", group.by = res_col(r), label = TRUE,
                   label.size = 3, raster = TRUE, raster.dpi = c(1200, 1200),
                   pt.size = 0.4) +
             coord_equal() + NoLegend() +
             labs(title = sprintf("res %.1f  |  %d clusters", r,
                                  nlevels(sub@meta.data[[res_col(r)]])))),
           ncol = 3),
         width = 15, height = 5 * ceiling(length(RES_SWEEP) / 3), dpi = 150)

  if (!res_col(PRIMARY_RES) %in% names(sub@meta.data))
    stop("primary resolution ", PRIMARY_RES, " was not computed")
  sub$subcluster <- droplevels(as.factor(sub@meta.data[[res_col(PRIMARY_RES)]]))
  Idents(sub) <- "subcluster"
  msg("  primary resolution ", PRIMARY_RES, ": ", nlevels(sub$subcluster), " subclusters")

  ggsave(file.path(OUTDIR, "umap_subclusters.png"),
         DimPlot(sub, reduction = "umap", group.by = "subcluster", label = TRUE,
                 label.size = 4, raster = TRUE, raster.dpi = c(1400, 1400),
                 pt.size = 0.4) + coord_equal() + NoLegend(),
         width = 9, height = 8.5, dpi = 200)
  for (g in c("sample", "diet", "sex", "parent_cluster"))
    ggsave(file.path(OUTDIR, paste0("umap_", g, ".png")),
           DimPlot(sub, reduction = "umap", group.by = g, raster = TRUE,
                   raster.dpi = c(1200, 1200), pt.size = 0.4, shuffle = TRUE) +
             coord_equal(),
           width = 9.5, height = 8, dpi = 150)

  ggsave(file.path(OUTDIR, "umap_qc.png"),
         FeaturePlot(sub, reduction = "umap",
                     features = c("nFeature_mouse", "nCount_mouse", "percent_mt"),
                     raster = TRUE, raster.dpi = c(1200, 1200), pt.size = 0.4, ncol = 3),
         width = 16, height = 5.5, dpi = 150)

  ## ---- composition: is any subcluster one animal? -------------------------

  cm <- as.data.table(sub@meta.data)[, .N, by = .(subcluster, mouse, sample, diet, sex)]
  fwrite(cm[order(subcluster, sample)], file.path(OUTDIR, "subcluster_composition.csv"))

  long_n <- cm[, .(n = sum(N)), by = .(subcluster, mouse)]
  long_n <- long_n[CJ(subcluster = subcluster, mouse = mouse, unique = TRUE),
                   on = .(subcluster, mouse)]
  long_n[is.na(n), n := 0L]
  long_n[, pct := 100 * n / sum(n), by = subcluster]
  fwrite(dcast(long_n, subcluster ~ mouse, value.var = "pct"),
         file.path(OUTDIR, "subcluster_mouse_percent.csv"))

  even <- 100 / uniqueN(long_n$mouse)
  ggsave(file.path(OUTDIR, "subcluster_mouse_heatmap.png"),
         ggplot(long_n, aes(mouse, subcluster, fill = pct)) +
           geom_tile(colour = "white", linewidth = 0.4) +
           geom_text(aes(label = fifelse(pct >= 2 * even, sprintf("%.0f", pct), ""),
                         colour = pct > 55), size = 2.6) +
           scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#3d3d3a"),
                               guide = "none") +
           scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, limits = c(0, 100),
                               name = "% of\ncluster") +
           scale_y_discrete(limits = rev(levels(long_n$subcluster))) +
           labs(x = NULL, y = "subcluster",
                subtitle = sprintf("even contribution would be %.1f%% per mouse", even)) +
           theme_minimal(base_size = 10) +
           theme(panel.grid = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1)),
         width = 9, height = 0.28 * nlevels(long_n$subcluster) + 2.5, dpi = 150)

  dom <- long_n[, .(top_mouse = mouse[which.max(pct)], top_pct = round(max(pct), 1)),
                by = subcluster][top_pct >= SINGLE_MOUSE_WARN]
  if (nrow(dom)) {
    msg("  WARNING: subcluster(s) >= ", SINGLE_MOUSE_WARN, "% from a single mouse -- ",
        "treat as artefact until shown otherwise:")
    print(dom)
  } else {
    msg("  no subcluster is >= ", SINGLE_MOUSE_WARN, "% from a single mouse")
  }

  sizes <- as.data.table(sub@meta.data)[, .N, by = subcluster][order(subcluster)]
  fwrite(sizes, file.path(OUTDIR, "subcluster_sizes.csv"))
  if (any(sizes$N < MIN_CLUSTER_CELLS))
    msg("  NOTE: ", sum(sizes$N < MIN_CLUSTER_CELLS), " subcluster(s) under ",
        MIN_CLUSTER_CELLS, " nuclei -- too small for the pseudobulk marker layer")

  ## ---- markers ------------------------------------------------------------

  if (length(Layers(sub[["RNA"]], search = "counts")) > 1) sub <- JoinLayers(sub)
  mk <- cluster_markers(sub, clusters = sub$subcluster, mouse = sub$mouse, log = msg)
  setnames(mk, "cluster_id", "subcluster")
  fwrite(mk, file.path(OUTDIR, "markers_all.csv"))
  top <- mk[order(subcluster, -concordant, -avg_log2FC)][, head(.SD, 20), by = subcluster]
  fwrite(top, file.path(OUTDIR, "markers_top20.csv"))

  conc <- mk[, .(n_wilcox = .N, n_concordant = sum(concordant, na.rm = TRUE)),
             by = subcluster][order(subcluster)]
  conc[, pct_concordant := round(100 * n_concordant / n_wilcox, 1)]
  fwrite(conc, file.path(OUTDIR, "concordance_by_subcluster.csv"))
  msg("  markers: ", nrow(mk), " rows, ", sum(mk$concordant, na.rm = TRUE),
      " concordant across both layers")
  print(conc)

  best <- mk[concordant == TRUE][order(subcluster, -avg_log2FC)][
    , .(gene = gene[1]), by = subcluster][!is.na(gene)]
  if (nrow(best)) {
    genes <- unique(best$gene)
    ggsave(file.path(OUTDIR, "dotplot_top_marker.png"),
           DotPlot(sub, features = genes, group.by = "subcluster") +
             coord_flip() + labs(x = NULL, y = "subcluster") +
             theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                   axis.text.y = element_text(size = 7)),
           width = 11, height = 0.22 * length(genes) + 3, dpi = 150, limitsize = FALSE)
  }

  objout <- file.path(OBJDIR, paste0("HFDseq_v2_sub_", lin, ".rds"))
  saveRDS(sub, objout, compress = FALSE)
  msg("  saved ", basename(objout), " (", round(file.size(objout) / 2^30, 2), " GB)")
  rm(sub); invisible(gc())
}

writeLines(capture.output(sessionInfo()), file.path(OUTROOT, "sessionInfo.txt"))
msg("done: ", OUTROOT)
