#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# figures_qc_to_annotation.R  --  every figure from QC through annotation, in
# one script, into one folder. Run by figures_qc_to_annotation.sh.
#
# Output goes to Results/reporting (asked for, 2026-09-07), alongside the
# doublet and decontamination figures reporting.R already writes there. The two
# scripts do not overlap: reporting.R covers the before/after CellBender
# comparison and the doublet violins, this one covers the pipeline figures from
# QC to annotation. Filenames are prefixed fig_NN_ so they sort in that order.
#
# HOUSEKEEPING, NOT A PIPELINE STEP, so it carries no number: it consumes
# finished objects and tables and changes nothing under Data/.
#
# NO PLOT CARRIES A TITLE OR SUBTITLE (asked for, 2026-09-07) -- every figure
# goes through nt(), which strips both, including the title Seurat's DimPlot
# and FeaturePlot set from group.by. Axis labels, legends and facet strips are
# kept: they are labels, not titles, and without them a faceted panel cannot be
# read at all. Where a multi-panel plot would have identified its panels by
# title alone, the panel's name is moved to the y-axis label instead.
#
# WHAT IT DRAWS, in pipeline order
# `annot` is shown as revised by 10_4_reclassify_from_azimuth.R when that
# script's output is present in this folder -- see section C.
#
#   A  ambient + doublets   step 4, from CSV     figs 01-02
#   B  QC                   step 5, from CSV     figs 03-09
#   C  embedding + clusters step 6, object+CSV   figs 10-18
#   D  markers              step 6, object+CSV   figs 19-20
#   E  annotation           step 7, object       figs 21-24
#
# Sections A and B read only the summary tables written by steps 4 and 5, so
# they cost nothing. Sections C-E need HFDseq_v2_annotated.rds (18 GB), which is
# read once and reused. Pass "tables" as the first argument to draw only A and B
# and skip the object entirely.
#
# TRAPS BAKED IN, all of them previously paid for in this project:
#   - step 6 scaled ONLY the 2,000 HVGs, and DoHeatmap silently drops features
#     missing from scale.data -- so the heatmap genes are re-scaled first or the
#     panel renders mostly blank with no error.
#   - DoHeatmap draws one column per cell; 247k columns is not a figure, so
#     cells are downsampled per cluster.
#   - geom_text(aes(y = -Inf)) on a scale_y_log10() panel gives log10(-Inf) =
#     NaN and silently drops every label. Top-anchored labels use y = Inf.
#   - the cluster x mouse heatmap fill is FIXED at c(0, 100). Auto-ranging makes
#     every cell mid-blue, because with even contribution the whole grid sits
#     near 8.3%.
#   - transgene features are `AAV-shared` etc inside a Seurat object (dash, not
#     underscore) -- CreateSeuratObject rewrites `_` to `-`.
#
# Usage (via sbatch figures_qc_to_annotation.sh):
#   figures_qc_to_annotation.sh            -> everything
#   figures_qc_to_annotation.sh tables     -> sections A and B only, no object
#
# Every output lands in Results/reporting.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(speckle)
  library(limma)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
OBJ_IN  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
D4      <- file.path(PROJECT, "Results", "4_ambient_doublets")
D5      <- file.path(PROJECT, "Results", "5_seurat_qc")
D6      <- file.path(PROJECT, "Results", "6_seurat_umap_markers")
DANN    <- file.path(PROJECT, "Results", "annotation_markers")
OUTDIR  <- file.path(PROJECT, "Results", "reporting")

SEQ_LO <- "#eef4fc"
SEQ_HI <- "#0d366b"
DIET_COLS <- c(SD = "#0f7ba2", HFD = "#dd5129")
SEX_COLS  <- c(F = "#bb292c", M = "#376795")
TOP_N     <- 5      # marker genes per cluster in the dotplot / heatmap
HEAT_CELLS <- 150   # cells per cluster drawn in the heatmap
SEED <- 42

args   <- commandArgs(trailingOnly = TRUE)
TABLES_ONLY <- length(args) >= 1 && identical(tolower(args[1]), "tables")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

source(file.path(PROJECT, "Code", "lib_celltype_palette.R"))

MOUSE_LEVELS <- paste0("Ech", 1:12)

# strip titles from a ggplot or a patchwork/Seurat multi-panel object
nt <- function(p) {
  if (inherits(p, "patchwork"))
    p & labs(title = NULL, subtitle = NULL, caption = NULL) &
      theme(plot.title = element_blank(), plot.subtitle = element_blank())
  else
    p + labs(title = NULL, subtitle = NULL, caption = NULL) +
      theme(plot.title = element_blank(), plot.subtitle = element_blank())
}

save_fig <- function(p, name, w, h, dpi = 300, pdf = TRUE) {
  p <- nt(p)
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h,
         dpi = dpi, limitsize = FALSE)
  if (pdf)
    ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h,
           limitsize = FALSE)
  msg("  ", name)
}

as_mouse <- function(x) factor(as.character(x), levels = MOUSE_LEVELS)

## ===========================================================================
## A. ambient correction and doublets  (step 4, tables only)
## ===========================================================================

msg("section A: ambient + doublets")

amb <- fread(file.path(D4, "cellbender_ambient_summary.csv"))
amb[, mouse := as_mouse(mouse)]
p01 <- ggplot(amb, aes(mouse, pct_removed, fill = diet)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", pct_removed)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = DIET_COLS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "mouse-gene UMIs removed by CellBender (%)") +
  facet_grid(~ sex, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA))
save_fig(p01, "fig_01_ambient_pct_removed", 9, 4.5)

dbl <- fread(file.path(D4, "doublet_by_mouse.csv"))
dbl[, mouse := as_mouse(mouse)]
p02 <- ggplot(dbl, aes(mouse, pct_doublet, fill = diet)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", pct_doublet)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = DIET_COLS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "nuclei called doublet by scDblFinder (%)") +
  facet_grid(~ sex, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA))
save_fig(p02, "fig_02_doublet_rate_per_animal", 9, 4.5)

## ===========================================================================
## B. QC  (step 5, per-cell table + cutoffs)
## ===========================================================================

msg("section B: QC")

qc <- fread(file.path(D5, "qc_metrics_per_cell.csv.gz"))
qc[, mouse := as_mouse(mouse)]
cut <- fread(file.path(D5, "qc_cutoffs_per_sample.csv"))
cut[, mouse := as_mouse(mouse)]

violin_qc <- function(metric, lo, hi, ylab) {
  b <- cut[, .(mouse, lo = get(lo), hi = get(hi))]
  ggplot(qc[get(metric) > 0], aes(mouse, get(metric))) +
    geom_violin(aes(fill = diet), scale = "width", linewidth = 0.25,
                colour = "#3d3d3a") +
    geom_errorbar(data = b, aes(x = mouse, ymin = lo, ymax = hi),
                  inherit.aes = FALSE, width = 0.55, colour = "#bb292c",
                  linewidth = 0.45) +
    scale_fill_manual(values = DIET_COLS, name = NULL) +
    scale_y_log10() +
    labs(x = NULL, y = ylab) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
}
save_fig(violin_qc("nFeature_mouse", "nFeature_lo", "nFeature_hi",
                   "detected mouse genes per nucleus (log10)"),
         "fig_03_qc_violin_nFeature", 10, 4.8)
save_fig(violin_qc("nCount_mouse", "nCount_lo", "nCount_hi",
                   "mouse-gene UMIs per nucleus (log10)"),
         "fig_04_qc_violin_nCount", 10, 4.8)

# This file has an embedded newline inside its quoted `lab` field. fread
# mis-detects the separator on it ("Detected 1 column names but the data has 6
# columns") and returns garbage; read.csv parses the quoting correctly.
fits <- as.data.table(read.csv(file.path(D5, "qc_count_vs_feature_fits.csv"),
                               check.names = FALSE))
fits <- merge(fits, unique(qc[, .(sample, mouse)]), by = "sample")
fits[, lab := sprintf("slope %.3f   R2 %.3f", slope, r2)]
p05 <- ggplot(qc[nCount_mouse > 0 & nFeature_mouse > 0],
              aes(nCount_mouse, nFeature_mouse)) +
  geom_point(aes(colour = qc_pass), size = 0.15, alpha = 0.25, shape = 16) +
  geom_abline(data = fits, aes(intercept = intercept, slope = slope),
              colour = "#dd5129", linewidth = 0.4) +
  geom_text(data = fits, aes(x = Inf, y = Inf, label = lab), hjust = 1.05,
            vjust = 1.6, size = 2.7, inherit.aes = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "#376795", `FALSE` = "grey70"),
                      name = NULL, labels = c("removed", "kept"),
                      guide = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ mouse, ncol = 4) +
  labs(x = "mouse-gene UMIs per nucleus (log10)",
       y = "detected mouse genes per nucleus (log10)") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p05, "fig_05_qc_scatter_count_vs_feature", 12, 9, pdf = FALSE)

mt_cap <- quantile(qc$percent_mt, 0.99, na.rm = TRUE)
p06 <- ggplot(qc[percent_mt <= mt_cap], aes(percent_mt)) +
  geom_histogram(aes(fill = diet), bins = 80) +
  geom_vline(xintercept = unique(cut$percent_mt_hi), colour = "#bb292c",
             linewidth = 0.4, linetype = "22") +
  scale_fill_manual(values = DIET_COLS, name = NULL) +
  scale_y_log10() +
  facet_wrap(~ mouse, ncol = 4, scales = "free_y") +
  labs(x = sprintf("percent mitochondrial UMIs (x capped at p99 = %.2f)", mt_cap),
       y = "nuclei (log10)") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p06, "fig_06_qc_hist_percent_mt", 12, 8, pdf = FALSE)

p07 <- ggplot(qc, aes(percent_mt, colour = mouse)) +
  stat_ecdf(linewidth = 0.4) +
  geom_vline(xintercept = unique(cut$percent_mt_hi), colour = "#bb292c",
             linewidth = 0.4, linetype = "22") +
  scale_x_continuous(limits = c(0, 20)) +
  scale_colour_viridis_d(option = "turbo", name = NULL) +
  labs(x = "percent mitochondrial UMIs", y = "cumulative fraction of nuclei") +
  theme_bw(base_size = 11)
save_fig(p07, "fig_07_qc_ecdf_percent_mt", 8, 5)

ret <- fread(file.path(D5, "qc_retention_per_sample.csv"))
ret[, mouse := as_mouse(mouse)]
p08 <- ggplot(ret, aes(mouse, pct_kept, fill = diet)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.1f", pct_kept)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = DIET_COLS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = NULL, y = "nuclei passing QC (%)") +
  facet_grid(~ sex, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA))
save_fig(p08, "fig_08_qc_retention_per_animal", 9, 4.5)

rem <- melt(fread(file.path(D5, "qc_removal_by_criterion.csv")),
            measure.vars = patterns("."), variable.name = "criterion",
            value.name = "n")
tot <- rem[criterion == "total", n]
rem <- rem[criterion != "total"]
rem[, criterion := factor(criterion, levels = rem[order(-n), as.character(criterion)])]
p09 <- ggplot(rem, aes(criterion, n)) +
  geom_col(fill = "#376795", width = 0.7) +
  geom_text(aes(label = sprintf("%s  (%.2f%%)", format(n, big.mark = ","),
                                100 * n / tot)), vjust = -0.4, size = 3.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = sprintf("nuclei removed (of %s entering QC)",
                             format(tot, big.mark = ","))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank())
save_fig(p09, "fig_09_qc_removal_by_criterion", 7.5, 4.5)

if (TABLES_ONLY) {
  writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
  msg("tables-only run finished: ", OUTDIR); quit(save = "no")
}

## ===========================================================================
## C. embedding and clusters  (step 6)
## ===========================================================================

msg("section C: embedding + clusters")
msg("reading ", basename(OBJ_IN), " (", round(file.size(OBJ_IN) / 2^30, 1), " GB)")
obj <- readRDS(OBJ_IN)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]], search = "counts")) > 1) obj <- JoinLayers(obj)
msg("loaded ", format(ncol(obj), big.mark = ","), " nuclei")

obj$mouse <- as_mouse(obj$mouse)
if (!"annot" %in% names(obj@meta.data))
  stop("no 'annot' column -- this object predates the manual annotation")

# 10_4 revises `annot` where the Azimuth mapping disagrees with the manual call:
# nuclei inside a neuronal blob that map onto Astro/Oligo/Micro-PVM/Endo/VLMC/
# OPC/Peri/SMC, and nuclei that map onto the other neuronal lineage. Applied
# here rather than baked into the object, so the hand-made annotation on disk
# stays intact. Absent the file, the figures show the manual call unchanged.
REVISED <- file.path(OUTDIR, "reclass_annot_revised.csv.gz")
if (file.exists(REVISED)) {
  rv <- fread(REVISED)[changed == TRUE]
  md <- data.table(barcode = colnames(obj), annot = as.character(obj$annot))
  md[rv, on = "barcode", annot := i.annot_azimuth]
  n_hit <- sum(md$barcode %in% rv$barcode)
  if (n_hit != nrow(rv))
    stop(nrow(rv) - n_hit, " revised barcodes are absent from the object")
  obj$annot <- md$annot
  msg("applied Azimuth reclassification to ", format(n_hit, big.mark = ","),
      " nuclei (", sprintf("%.2f%%", 100 * n_hit / ncol(obj)), ")")
} else {
  msg("no reclassification file -- figures show the manual annot unchanged")
}

RAST <- list(raster = TRUE, raster.dpi = c(1400, 1400), pt.size = 0.4)
dim_plot <- function(group, cols = NULL, label = FALSE, legend = TRUE) {
  p <- do.call(DimPlot, c(list(obj, reduction = "umap", group.by = group,
                               label = label, label.size = 3.4, repel = TRUE),
                          RAST))
  if (!is.null(cols)) p <- p + scale_colour_manual(values = cols, name = NULL)
  p <- p + coord_equal() + theme_bw(base_size = 11)
  if (!legend) p <- p + NoLegend()
  p
}

save_fig(dim_plot("seurat_clusters", label = TRUE, legend = FALSE),
         "fig_10_umap_clusters", 8, 8, pdf = FALSE)
save_fig(dim_plot("diet", cols = DIET_COLS), "fig_11_umap_diet", 8.6, 8, pdf = FALSE)
save_fig(dim_plot("sex",  cols = SEX_COLS),  "fig_12_umap_sex",  8.6, 8, pdf = FALSE)
save_fig(dim_plot("sample"), "fig_13_umap_sample", 9.6, 8, pdf = FALSE)

# Panel titles would normally carry the metric name; titles are stripped, so
# the name is moved onto each panel's y axis instead.
qc_feats <- intersect(c("nFeature_mouse", "nCount_mouse", "percent_mt"),
                      names(obj@meta.data))
fp <- lapply(qc_feats, function(f) {
  hi <- quantile(obj@meta.data[[f]], 0.99, na.rm = TRUE)
  do.call(FeaturePlot, c(list(obj, features = f, reduction = "umap",
                              max.cutoff = hi), RAST)) +
    scale_colour_gradient(low = SEQ_LO, high = SEQ_HI, name = NULL) +
    coord_equal() + labs(y = f) + theme_bw(base_size = 10)
})
save_fig(wrap_plots(fp, nrow = 1), "fig_14_umap_qc_metrics",
         6 * length(fp), 5.6, pdf = FALSE)

um <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "barcode")
setnames(um, 2:3, c("UMAP_1", "UMAP_2"))
um[, `:=`(mouse = obj$mouse, annot = as.character(obj$annot))]
bg <- um[, .(UMAP_1, UMAP_2)]
acols <- celltype_cols(sort(unique(um$annot)))
p15 <- ggplot(um, aes(UMAP_1, UMAP_2)) +
  geom_point(data = bg, colour = "grey88", size = 0.30, shape = 16) +
  geom_point(aes(colour = annot), size = 0.35, shape = 16) +
  scale_colour_manual(values = acols, name = NULL,
                      guide = guide_legend(override.aes = list(size = 2.6))) +
  facet_wrap(~ mouse, ncol = 4) + coord_equal() +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p15, "fig_15_umap_facet_by_animal", 14, 11, pdf = FALSE)

heat_long <- function(f, idcol) {
  d <- fread(file.path(D6, f))
  m <- melt(d, id.vars = idcol, variable.name = "mouse", value.name = "pct")
  setnames(m, idcol, "cluster")
  m[, mouse := as_mouse(mouse)]
  m[, cluster := factor(as.character(cluster),
                        levels = as.character(sort(as.integer(unique(cluster)))))]
  m
}
cm <- heat_long("cluster_mouse_percent.csv", "seurat_clusters")
p16 <- ggplot(cm, aes(mouse, cluster, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = fifelse(pct >= 17, sprintf("%.0f", pct), ""),
                colour = pct > 55), size = 2.5) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#3d3d3a"), guide = "none") +
  scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, limits = c(0, 100),
                      name = "% of\ncluster") +
  scale_y_discrete(limits = rev(levels(cm$cluster))) +
  labs(x = NULL, y = "cluster") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank())
save_fig(p16, "fig_16_cluster_by_mouse_heatmap", 7.5, 0.26 * uniqueN(cm$cluster) + 2.5)

mc <- heat_long("mouse_cluster_percent.csv", "seurat_clusters")
ccols <- setNames(scales::hue_pal()(uniqueN(mc$cluster)), levels(mc$cluster))
p17 <- ggplot(mc, aes(mouse, pct, fill = cluster)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = fifelse(pct >= 4, as.character(cluster), "")),
            position = position_stack(vjust = 0.5), size = 2.3) +
  scale_fill_manual(values = ccols, name = "cluster") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% of that animal's nuclei") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank())
save_fig(p17, "fig_17_mouse_cluster_bars", 11, 6)

if ("pca" %in% Reductions(obj)) {
  sd <- Stdev(obj, reduction = "pca")
  p18 <- ggplot(data.table(PC = seq_along(sd), sd = sd), aes(PC, sd)) +
    geom_point(size = 1.4, colour = "#376795") +
    geom_vline(xintercept = 30, colour = "#bb292c", linetype = "22", linewidth = 0.4) +
    labs(x = "principal component", y = "standard deviation") +
    theme_bw(base_size = 11)
  save_fig(p18, "fig_18_pca_elbow", 7, 4.5)
}

## ===========================================================================
## D. cluster markers  (step 6)
## ===========================================================================

msg("section D: markers")

mk <- fread(file.path(D6, "markers_all.csv"))
mk <- mk[avg_log2FC > 0]
setorder(mk, seurat_clusters, -concordant, -avg_log2FC)
top <- mk[, head(.SD, TOP_N), by = seurat_clusters]
genes <- unique(top$gene)
genes <- intersect(genes, rownames(obj))
msg("  ", length(genes), " marker genes over ", uniqueN(top$seurat_clusters), " clusters")

Idents(obj) <- "seurat_clusters"
p19 <- DotPlot(obj, features = genes, group.by = "seurat_clusters") +
  scale_colour_gradient(low = SEQ_LO, high = SEQ_HI) +
  labs(x = NULL, y = "cluster") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5))
save_fig(p19, "fig_19_dotplot_top_markers_by_cluster",
         0.16 * length(genes) + 3, 0.22 * uniqueN(top$seurat_clusters) + 3)

# step 6 scaled only the 2,000 HVGs; DoHeatmap drops anything missing from
# scale.data without erroring, so the marker genes are re-scaled here.
need <- setdiff(genes, rownames(GetAssayData(obj, layer = "scale.data")))
if (length(need)) {
  msg("  re-scaling ", length(need), " genes absent from scale.data")
  obj <- ScaleData(obj, features = union(genes, rownames(
    GetAssayData(obj, layer = "scale.data"))), verbose = FALSE)
}
present <- intersect(genes, rownames(GetAssayData(obj, layer = "scale.data")))
msg("  ", length(present), "/", length(genes), " genes present in scale.data")
set.seed(SEED)
small <- subset(obj, downsample = HEAT_CELLS)
p20 <- DoHeatmap(small, features = present, group.by = "seurat_clusters",
                 size = 2.6, angle = 90, raster = TRUE) +
  scale_fill_gradient2(low = "#0f7ba2", mid = "white", high = "#bb292c",
                       midpoint = 0, name = NULL) +
  theme(axis.text.y = element_text(size = 4.5))
save_fig(p20, "fig_20_heatmap_top_markers_by_cluster",
         0.30 * uniqueN(small$seurat_clusters) + 6, 0.09 * length(present) + 4,
         pdf = FALSE)
rm(small); invisible(gc())

## ===========================================================================
## E. annotation  (step 7)
## ===========================================================================

msg("section E: annotation")

# The same marker panel and the same banded layout the manual annotation was
# read off, so the handover figure and the call that produced `annot` agree --
# with ONE deliberate exception: Eno2 is dropped from the excitatory band (asked
# for, 2026-09-08). It is a pan-neuronal enolase, not an excitatory marker, and
# it fired across both neuronal bands. It stays in the panel 7_annotation_
# featureplots.R plots, which is the record of what the manual call actually
# looked at; only the figure changes.
markers <- list(
  "Excitatory Neurons" = c("Slc17a7","Slc17a6","Slc17a8","Camk2a","Grin2b","Tubb3","Syp"),
  "Inhibitory Neurons" = c("Gad1","Gad2","Slc32a1","Pvalb","Sst"),
  "Microglia"          = c("Aif1","Itgam","Ptprc","Cx3cr1","Tmem119"),
  "Astrocytes"         = c("Gfap","Aldh1l1","S100b","Sox9","Slc1a2"),
  "Oligodendrocytes"   = c("Olig2","Sox10","Mbp","Plp1","Myrf"),
  "OPC"                = c("Pdgfra","Cspg4","Olig2","Sox10"),
  "VLCM"               = c("Dcn","Bmp6","Cped1","Cfh","Slc6a13"),
  "Endothelial Cell"   = c("Flt1","Ly6c1","Pecam1","Cdh5","Cldn5"),
  "Meningeal Cell"     = c("Bnc2","Slc47a1","Slc6a13"),
  "Pericytes"          = c("Kcnj8","Vtn","Abcc9","Higd1b","Cd248"),
  "Smooth Muscle"      = c("Acta2","Myh11","Tagln","Myocd","Lmod1"),
  "Macrophage"         = c("Mrc1","Pf4","Lyz2","Ms4a7","F13a1","Lyve1"))
markers <- setNames(lapply(seq_along(markers), function(i)
  setdiff(markers[[i]], unlist(markers[seq_len(i - 1)]))), names(markers))
markers <- lapply(markers, function(g) intersect(g, rownames(obj)))
celltypes <- names(markers)
agenes <- unlist(markers, use.names = FALSE)

ct_levels <- rev(c(celltypes, setdiff(unique(as.character(obj$annot)), celltypes)))
obj$annot <- factor(as.character(obj$annot), levels = ct_levels)
acols <- celltype_cols(levels(obj$annot))

save_fig(dim_plot("annot", cols = acols, label = TRUE),
         "fig_21_umap_annotation", 10, 8, pdf = FALSE)

nmk   <- lengths(markers)
bands <- data.frame(xmin = cumsum(c(0, head(nmk, -1))) + 0.5, xmax = cumsum(nmk) + 0.5)
pd <- DotPlot(obj, features = agenes, group.by = "annot")
pd$layers <- c(annotate("rect", xmin = bands$xmin, xmax = bands$xmax,
                        ymin = -Inf, ymax = Inf,
                        fill = rep(c("grey90", "white"), length.out = nrow(bands)),
                        alpha = 0.5),
               pd$layers)
p22 <- pd +
  annotate("text", x = (bands$xmin + bands$xmax) / 2, y = Inf, label = celltypes,
           vjust = rep(c(-2.4, -0.7), length.out = length(celltypes)), size = 3) +
  scale_colour_gradient(low = SEQ_LO, high = SEQ_HI) +
  coord_cartesian(clip = "off") + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.margin = margin(t = 48, r = 6, b = 6, l = 6))
save_fig(p22, "fig_22_dotplot_annotation_markers", 16, 8, pdf = FALSE)

comp <- as.data.table(obj@meta.data)[, .N, by = .(mouse, diet, sex, annot)]
comp[, pct := 100 * N / sum(N), by = mouse]
fwrite(dcast(comp, mouse + diet + sex ~ annot, value.var = "pct", fill = 0),
       file.path(OUTDIR, "annotation_composition_per_animal.csv"))
p23 <- ggplot(comp, aes(mouse, pct, fill = annot)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = acols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% of that animal's nuclei") +
  facet_grid(~ sex, scales = "free_x", space = "free_x") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA))
save_fig(p23, "fig_23_annotation_composition_by_animal", 11, 6)

# The composition test for THIS figure, run here rather than read from step 8:
# step 8 tested the pre-revision `annot`, and the labels plotted above are the
# revised ones, so its p-values would not describe these points. Same test as
# step 8 and as 10_5 -- propeller, asin-sqrt, ~ 0 + diet + sex, HFD - SD,
# robust empirical-Bayes shrinkage, n = 12 animals.
md24 <- as.data.table(obj@meta.data)[, .(sample = as.character(sample),
                                         diet, sex, annot = as.character(annot))]
md24[, `:=`(diet = factor(diet, levels = c("SD", "HFD")),
            sex  = factor(sex,  levels = c("F", "M")))]
props   <- getTransformedProps(clusters = md24$annot, sample = md24$sample,
                               transform = "asin")
animals <- unique(md24[, .(sample, diet, sex)])
animals <- animals[match(colnames(props$TransformedProps), sample)]
stopifnot(identical(animals$sample, colnames(props$TransformedProps)))
des <- model.matrix(~ 0 + diet + sex, data = animals)
ctr <- setNames(rep(0, ncol(des)), colnames(des))
ctr["dietHFD"] <- 1; ctr["dietSD"] <- -1
pr <- as.data.table(propeller.ttest(props, des, contrasts = ctr, robust = TRUE,
                                    trend = FALSE, sort = FALSE),
                    keep.rownames = "annot")
fwrite(pr[order(P.Value)], file.path(OUTDIR, "annotation_composition_propeller.csv"))
msg("  propeller on the revised annot: ", pr[FDR < 0.05, .N],
    " cell type(s) at FDR < 0.05")
pr[, `:=`(lab = sprintf("p = %s\nFDR = %s", signif(P.Value, 2), signif(FDR, 2)),
          sig = FDR < 0.05)]
pr[, annot := factor(annot, levels = levels(comp$annot))]

p24 <- ggplot(comp, aes(diet, pct, colour = annot)) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.8,
             show.legend = FALSE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.45, linewidth = 0.3,
               colour = "#3d3d3a", show.legend = FALSE) +
  geom_text(data = pr, aes(x = 1.5, y = Inf, label = lab, colour = sig),
            inherit.aes = FALSE, vjust = 1.25, size = 2.8, lineheight = 0.95,
            show.legend = FALSE) +
  facet_wrap(~ annot, scales = "free_y", ncol = 5) +
  scale_colour_manual(values = c(acols, `TRUE` = "#bb292c", `FALSE` = "#3d3d3a")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  labs(x = NULL, y = "% of that animal's nuclei") +
  guides(colour = "none") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
save_fig(p24, "fig_24_annotation_composition_by_diet", 13, 7)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
figs <- list.files(OUTDIR, pattern = "^fig_.*\\.png$")
msg("done: ", length(figs), " figures in ", OUTDIR)
