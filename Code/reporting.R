#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# reporting.R  --  figures for the write-up, run by reporting.sh.
#
# HOUSEKEEPING, NOT A PIPELINE STEP, so it carries no number: it consumes
# finished objects and changes nothing. Run it whenever a figure needs
# regenerating; nothing downstream depends on its output.
#
# Two independent sections, chosen by the first argument (default "all"):
#
#  A  "doublet"  Violins of UMI and detected-gene counts in the nuclei
#                scDblFinder called doublet vs singlet. Read from
#                HFDseq_v2_ambient_doublets.rds because that is the LAST
#                object that still contains the doublets -- 5_seurat_qc.R
#                subsets to singlets before it writes anything, so
#                qc_metrics_per_cell.csv.gz cannot answer this question.
#                The comparison is therefore made on ALL CellBender nuclei,
#                i.e. at the moment the call was made and before any QC
#                cutoff, which is the honest denominator: filtering first
#                would remove the high-count tail the plot is about.
#
#                Counts are MOUSE-GENE counts (the 4 transgene rows removed),
#                the same nCount_mouse / nFeature_mouse that 5_seurat_qc.R
#                cuts on, so the numbers in this figure and in the QC tables
#                are the same numbers.
#
#  B  "decon"    UMAP feature plots of a few genes BEFORE and AFTER CellBender.
#                "Before" is not stored anywhere, so it is rebuilt here: the
#                12 raw cellranger matrices are re-read and restricted to the
#                barcodes that survived into the plotted object. Both sides
#                are log-normalised HERE, each on its OWN mouse-gene library
#                size (pre on the raw total, post on the corrected total),
#                because that is what the two objects would look like if you
#                had normalised each one on its own -- normalising the pre
#                matrix on post totals would manufacture a difference.
#                Both panels of a gene share one colour scale, cut at the p99
#                of the PRE non-zero values; anything above is clamped.
#
# Neither section writes an object or touches Data/.
#
# Usage (via sbatch reporting.sh):
#   reporting.sh                       -> both sections, default genes
#   reporting.sh doublet
#   reporting.sh decon Plp1,Mbp,Snap25
#   reporting.sh decon Plp1,Mbp /path/to/other_object.rds
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(scattermore)
})

source("/scratch/user/uqschauq/Insular/Code/lib_celltype_palette.R")

PROJECT <- "/scratch/user/uqschauq/Insular"
OUTDIR  <- file.path(PROJECT, "Results", "reporting")
AGGCSV  <- file.path(PROJECT, "Data", "Aggr_v2", "HFDseq_aggr_v2", "outs", "aggregation.csv")
IN_DBL  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_ambient_doublets.rds")

# Default object for section B: the manually annotated object written at the end
# of 7_annotation_featureplots.R, whose `annot` column is the working cell-type
# assignment. Pass another object as the 3rd argument to override. Either way it
# must be derived from the current HFDseq_v2_umap.rds, and the script warns if
# the object it is given is older than that.
IN_OBJ_DEFAULT <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
IN_UMAP        <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_umap.rds")

TRANSGENES <- c("EGFP", "tdTomato", "EBFP2", "AAV_shared")

# Seurat rewrote AAV_shared to AAV-shared when the object was built (see
# Results/4_ambient_doublets/feature_name_changes.csv), so the object and the
# raw h5 files need different spellings of the same four features. 5_seurat_qc.R
# does the same substitution for the same reason.
TRANSGENES_OBJ <- gsub("_", "-", TRANSGENES, fixed = TRUE)

# Genes for the pre/post panels. Chosen as high-expressing, cell-type-restricted
# transcripts that are the usual ambient sources in brain nuclei preps -- myelin
# (Plp1, Mbp), neuronal (Snap25), astrocytic (Slc1a2), choroid plexus (Ttr) and
# blood (Hbb-bs). Ambient shows up as a low, even wash of these over every
# cluster; decontamination should erase the wash and leave the island.
GENES_DEFAULT <- c("Plp1", "Mbp", "Snap25", "Slc1a2", "Ttr", "Hbb-bs")

SCALE_FACTOR <- 1e4      # matches NormalizeData() in 6_seurat_umap_markers.R
Q_HI         <- 0.99     # colour-scale cut, as in 7_annotation_featureplots.R
CUT_FLOOR    <- 1
COL_LO       <- "grey92"
COL_HI       <- "#08306b"
COL_SINGLET  <- "#0d366b"
COL_DOUBLET  <- "#d1341f"
PT_PX        <- 2.2      # scattermore point size, in pixels
RAST_PX      <- c(1400, 1400)
SEED         <- 42

args    <- commandArgs(trailingOnly = TRUE)
SECTION <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "all"
GENES   <- if (length(args) >= 2 && nzchar(args[2]))
             trimws(strsplit(args[2], ",", fixed = TRUE)[[1]]) else GENES_DEFAULT
IN_OBJ  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else IN_OBJ_DEFAULT

if (!SECTION %in% c("all", "doublet", "decon"))
  stop("section must be one of: all, doublet, decon (got '", SECTION, "')")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
msg("section = ", SECTION, " | out = ", OUTDIR)

theme_set(theme_bw(base_size = 11))

## ===========================================================================
## A. doublets carry more UMIs and more genes
## ===========================================================================

if (SECTION %in% c("all", "doublet")) {

  stopifnot(file.exists(IN_DBL))
  msg("reading ", basename(IN_DBL), " (all CellBender nuclei, doublets flagged)")
  obj <- readRDS(IN_DBL)
  stopifnot(all(c("sample", "mouse", "diet", "sex", "dbl_class") %in% names(obj@meta.data)))

  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
  tg  <- intersect(TRANSGENES_OBJ, rownames(cnt))
  if (length(tg) != length(TRANSGENES_OBJ))
    msg("WARNING: transgene rows absent from the object: ",
        paste(setdiff(TRANSGENES_OBJ, tg), collapse = ", "))

  # Same construction as 5_seurat_qc.R: subtract the transgene block from the
  # totals Seurat already computed, rather than re-summing 33,700 rows.
  tg_umi      <- if (length(tg)) Matrix::colSums(cnt[tg, , drop = FALSE]) else 0
  tg_detected <- if (length(tg)) Matrix::colSums(cnt[tg, , drop = FALSE] > 0) else 0
  rm(cnt); invisible(gc())

  d <- as.data.table(obj@meta.data, keep.rownames = "barcode")
  d[, `:=`(nCount_mouse   = as.numeric(nCount_RNA)  - as.numeric(tg_umi),
           nFeature_mouse = as.integer(nFeature_RNA) - as.integer(tg_detected))]
  d[, dbl_class := factor(dbl_class, levels = c("singlet", "doublet"))]
  rm(obj); invisible(gc())

  msg(nrow(d), " nuclei | ", d[dbl_class == "doublet", .N], " doublets (",
      round(100 * mean(d$dbl_class == "doublet"), 2), "%)")

  # Kept short: the longer form ran off the right edge of the panel, and
  # "mouse genes" is already in the figure subtitle.
  METRICS <- c(nCount_mouse = "UMIs per nucleus",
               nFeature_mouse = "genes detected per nucleus")

  # Summary table: overall and per animal, medians + a rank-based effect size.
  # With n in the hundreds of thousands the Wilcoxon p is uninformative (it is
  # ~0 for any real difference), so auc -- P(doublet > singlet) -- is the number
  # to quote alongside the medians.
  stat_block <- function(dd, metric) {
    x <- dd[dbl_class == "doublet",  get(metric)]
    y <- dd[dbl_class == "singlet",  get(metric)]
    if (!length(x) || !length(y))
      return(data.table(metric = metric, n_singlet = length(y), n_doublet = length(x)))
    w <- suppressWarnings(wilcox.test(x, y))
    num <- function(v) unname(as.numeric(v))
    data.table(
      metric        = metric,
      n_singlet     = length(y),               n_doublet     = length(x),
      median_singlet = num(median(y)),         median_doublet = num(median(x)),
      q25_singlet   = num(quantile(y, .25)),   q25_doublet   = num(quantile(x, .25)),
      q75_singlet   = num(quantile(y, .75)),   q75_doublet   = num(quantile(x, .75)),
      mean_singlet  = round(mean(y), 1),       mean_doublet  = round(mean(x), 1),
      ratio_median  = round(num(median(x)) / num(median(y)), 3),
      auc           = round(as.numeric(w$statistic) /
                            (as.numeric(length(x)) * as.numeric(length(y))), 4),
      wilcox_p      = w$p.value)
  }

  stats <- rbindlist(c(
    lapply(names(METRICS), function(m) cbind(scope = "all animals", sample = "ALL",
                                             stat_block(d, m))),
    lapply(names(METRICS), function(m)
      d[, cbind(scope = "per animal", stat_block(.SD, m)), by = sample])),
    fill = TRUE)
  fwrite(stats, file.path(OUTDIR, "doublet_umi_gene_stats.csv"))
  msg("pooled medians:")
  print(stats[scope == "all animals",
              .(metric, median_singlet, median_doublet, ratio_median, auc)])

  # Long form for plotting. Zeros would drop off a log axis with a warning;
  # clamp them at 1 and say how many were clamped.
  long <- melt(d[, .(barcode, sample, mouse, diet, dbl_class,
                     nCount_mouse = as.numeric(nCount_mouse),
                     nFeature_mouse = as.numeric(nFeature_mouse))],
               id.vars = c("barcode", "sample", "mouse", "diet", "dbl_class"),
               variable.name = "metric", value.name = "value")
  nz <- long[value < 1, .N]
  if (nz) msg("clamping ", nz, " value(s) < 1 to 1 for the log axis")
  long[, value := pmax(value, 1)]
  long[, metric := factor(metric, levels = names(METRICS), labels = METRICS)]

  meds <- long[, .(value = median(value), n = .N), by = .(metric, dbl_class)]
  meds_s <- long[, .(value = median(value)), by = .(metric, sample, dbl_class)]

  fills <- c(singlet = COL_SINGLET, doublet = COL_DOUBLET)
  base_violin <- function(dd) {
    ggplot(dd, aes(dbl_class, value, fill = dbl_class)) +
      geom_violin(scale = "width", linewidth = 0.2, colour = "grey25") +
      geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2, fill = "white") +
      scale_fill_manual(values = fills, guide = "none") +
      scale_y_log10(labels = function(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)) +
      labs(x = NULL, y = NULL)
  }

  sub_pooled <- sprintf(
    "all %s CellBender nuclei, before any QC cutoff | %s doublets (%.2f%%) | transgene features excluded",
    format(nrow(d), big.mark = ","), format(d[dbl_class == "doublet", .N], big.mark = ","),
    100 * mean(d$dbl_class == "doublet"))

  panels <- lapply(levels(long$metric), function(m) {
    dd <- long[metric == m]
    mm <- meds[metric == m]
    st <- stats[scope == "all animals" & metric == names(METRICS)[METRICS == m]]
    # No median label drawn here: the subtitle already carries both medians, and
    # a label placed anywhere inside the panel lands on the violin or the box.
    # The per-animal facets below keep theirs -- they have no subtitle.
    base_violin(dd) +
      labs(title = m,
           subtitle = sprintf("median %s vs %s  (x%.2f) | P(doublet > singlet) = %.3f",
                              format(round(mm[dbl_class == "doublet", value]), big.mark = ","),
                              format(round(mm[dbl_class == "singlet", value]), big.mark = ","),
                              st$ratio_median, st$auc)) +
      theme(plot.title = element_text(size = 12, face = "bold"),
            plot.subtitle = element_text(size = 9, colour = "grey35"))
  })

  p <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title = "Doublet nuclei carry more UMIs and more genes",
      subtitle = sub_pooled,
      caption = "scDblFinder, called per animal on CellBender-corrected mouse genes; box = IQR, label = median",
      theme = theme(plot.title = element_text(size = 16, face = "bold"),
                    plot.subtitle = element_text(size = 10, colour = "grey35")))
  ggsave(file.path(OUTDIR, "doublet_violin_pooled.png"), p, width = 9, height = 5.5, dpi = 300)
  ggsave(file.path(OUTDIR, "doublet_violin_pooled.pdf"), p, width = 9, height = 5.5)
  msg("written doublet_violin_pooled.png/.pdf")

  for (i in seq_along(METRICS)) {
    m  <- METRICS[i]
    dd <- long[metric == m]
    ps <- base_violin(dd) +
      geom_text(data = meds_s[metric == m], inherit.aes = FALSE,
                aes(dbl_class, value, label = format(round(value), big.mark = ",")),
                nudge_x = 0.32, vjust = -0.4, size = 2.6, colour = "grey15") +
      facet_wrap(~ sample, ncol = 4) +
      labs(title = paste0(m, ", per animal"),
           subtitle = sub_pooled,
           caption = "labels are medians") +
      theme(strip.background = element_blank(),
            plot.title = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9, colour = "grey35"))
    f <- file.path(OUTDIR, sprintf("doublet_violin_per_animal_%s.png", names(METRICS)[i]))
    ggsave(f, ps, width = 12, height = 8, dpi = 300)
    ggsave(sub("\\.png$", ".pdf", f), ps, width = 12, height = 8)
    msg("written ", basename(f), "/.pdf")
  }

  rm(long, d); invisible(gc())
}

## ===========================================================================
## B. before / after CellBender, on the UMAP
## ===========================================================================

if (SECTION %in% c("all", "decon")) {

  stopifnot(file.exists(IN_OBJ), file.exists(AGGCSV))
  if (file.exists(IN_UMAP) && !identical(IN_OBJ, IN_UMAP) &&
      file.mtime(IN_OBJ) < file.mtime(IN_UMAP))
    msg("WARNING: ", basename(IN_OBJ), " is OLDER than ", basename(IN_UMAP),
        " -- check it was rebuilt on the current clustering")

  msg("reading ", basename(IN_OBJ))
  obj <- readRDS(IN_OBJ)
  DefaultAssay(obj) <- "RNA"
  stopifnot("umap" %in% Reductions(obj))
  bc <- colnames(obj)
  msg("loaded ", format(length(bc), big.mark = ","), " nuclei x ", nrow(obj), " features")

  # Group column for the reference UMAP: whatever annotation exists, else the
  # cluster id. `annot` is the manual pass in 7_annotation_featureplots.R,
  # `cluster_label` the one 7_seurat_annotation.R writes.
  GROUP <- c("annot", "cluster_label", "cell_class", "seurat_clusters")
  GROUP <- GROUP[GROUP %in% names(obj@meta.data)][1]
  if (is.na(GROUP)) stop("no grouping column found in the object metadata")
  msg("grouping by '", GROUP, "' (", length(unique(obj@meta.data[[GROUP]])), " levels)")

  post_cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
  genes <- GENES[GENES %in% rownames(post_cnt)]
  if (length(setdiff(GENES, genes)))
    msg("WARNING: dropping gene(s) absent from the object: ",
        paste(setdiff(GENES, genes), collapse = ", "))
  if (!length(genes)) stop("none of the requested genes are in the object")
  msg("genes: ", paste(genes, collapse = ", "))

  tg <- intersect(TRANSGENES_OBJ, rownames(post_cnt))
  post <- as.matrix(post_cnt[genes, , drop = FALSE])
  post_total <- if ("nCount_mouse" %in% names(obj@meta.data)) {
    as.numeric(obj$nCount_mouse)
  } else {
    as.numeric(obj$nCount_RNA) -
      (if (length(tg)) Matrix::colSums(post_cnt[tg, , drop = FALSE]) else 0)
  }
  rm(post_cnt); invisible(gc())

  ## ---- rebuild the pre-CellBender counts for these barcodes ---------------

  agg <- fread(AGGCSV)
  setnames(agg, gsub("\r", "", names(agg), fixed = TRUE))
  agg[, gem_group := .I]

  gem <- as.integer(sub("^.*-", "", bc))
  if (anyNA(gem)) stop("barcodes without a -N GEM suffix")
  if (!all(unique(gem) %in% agg$gem_group)) stop("GEM groups outside aggregation.csv")

  pre       <- matrix(0, nrow = length(genes), ncol = length(bc),
                      dimnames = list(genes, bc))
  pre_total <- setNames(rep(NA_real_, length(bc)), bc)

  for (i in sort(unique(gem))) {
    s     <- agg$sample_id[i]
    rawh5 <- file.path(PROJECT, "Data", "Counts_v2", s, "outs", "raw_feature_bc_matrix.h5")
    if (!file.exists(rawh5)) stop("missing raw matrix: ", rawh5)

    cells <- bc[gem == i]
    raw <- Read10X_h5(rawh5)
    if (is.list(raw)) raw <- raw[["Gene Expression"]]

    # The object's barcodes are <16bp>-<gem>; the raw matrix uses cellranger's
    # per-sample <16bp>-1. Rebuild the raw name and refuse to guess if it misses.
    want <- paste0(sub("-[0-9]+$", "", cells), "-1")
    if (!all(want %in% colnames(raw)))
      stop(sum(!want %in% colnames(raw)), " barcode(s) of ", s, " absent from ", rawh5)
    missing_g <- setdiff(genes, rownames(raw))
    if (length(missing_g))
      stop("gene(s) absent from the raw matrix of ", s, ": ",
           paste(missing_g, collapse = ", "))

    tg_raw <- intersect(TRANSGENES, rownames(raw))
    rk <- raw[, want, drop = FALSE]
    pre[, cells] <- as.matrix(rk[genes, , drop = FALSE])
    pre_total[cells] <- Matrix::colSums(rk) -
      (if (length(tg_raw)) Matrix::colSums(rk[tg_raw, , drop = FALSE]) else 0)
    rm(raw, rk); invisible(gc())
    msg("  ", s, " (gem ", i, "): ", format(length(cells), big.mark = ","),
        " nuclei re-read from the raw matrix")
  }
  if (anyNA(pre_total)) stop("barcodes left without a raw library size")

  ## ---- normalise each side on its own library size ------------------------

  lognorm <- function(m, total) log1p(sweep(m, 2, pmax(total, 1), "/") * SCALE_FACTOR)
  pre_n  <- lognorm(pre,  pre_total)
  post_n <- lognorm(post, post_total)

  emb <- as.data.table(Embeddings(obj, "umap")[, 1:2], keep.rownames = "barcode")
  setnames(emb, 2:3, c("UMAP_1", "UMAP_2"))
  emb[, group := as.character(obj@meta.data[[GROUP]])]
  stopifnot(identical(emb$barcode, bc))

  ## ---- per-gene numbers ---------------------------------------------------

  grp <- emb$group
  gstats <- rbindlist(lapply(genes, function(g) {
    a <- pre[g, ]; b <- post[g, ]
    an <- pre_n[g, ]; bn <- post_n[g, ]
    top <- names(which.max(tapply(bn, grp, mean)))   # where the gene really lives
    off <- grp != top
    data.table(
      gene = g, n_nuclei = length(a),
      umi_pre = sum(a), umi_post = sum(b),
      pct_umi_removed = round(100 * (sum(a) - sum(b)) / max(sum(a), 1), 2),
      pct_detected_pre  = round(100 * mean(a > 0), 2),
      pct_detected_post = round(100 * mean(b > 0), 2),
      mean_lognorm_pre  = round(mean(an), 4),
      mean_lognorm_post = round(mean(bn), 4),
      top_group = top,
      pct_detected_offtarget_pre  = round(100 * mean(a[off] > 0), 2),
      pct_detected_offtarget_post = round(100 * mean(b[off] > 0), 2),
      pct_detected_ontarget_pre   = round(100 * mean(a[!off] > 0), 2),
      pct_detected_ontarget_post  = round(100 * mean(b[!off] > 0), 2))
  }))
  gstats[, offtarget_drop_pct := round(100 * (pct_detected_offtarget_pre -
                                              pct_detected_offtarget_post) /
                                       pmax(pct_detected_offtarget_pre, 1e-9), 1)]
  fwrite(gstats, file.path(OUTDIR, "decon_gene_stats.csv"))
  msg("per-gene ambient removal:")
  print(gstats[, .(gene, umi_pre, umi_post, pct_umi_removed,
                   pct_detected_pre, pct_detected_post,
                   top_group, pct_detected_offtarget_pre, pct_detected_offtarget_post)])

  ## ---- plots ---------------------------------------------------------------

  # Cluster ids get Seurat's default colours; a cell-type column gets the shared
  # palette, so this UMAP matches the composition figure and the annotated UMAP.
  ref_cols <- if (GROUP == "seurat_clusters") NULL
              else celltype_cols(levels(factor(obj@meta.data[[GROUP]])))
  ref <- DimPlot(obj, reduction = "umap", group.by = GROUP, label = TRUE,
                 label.size = 3, repel = TRUE, raster = TRUE, cols = ref_cols,
                 raster.dpi = RAST_PX, pt.size = 0.35) +
    coord_equal() + NoLegend() +
    labs(title = sprintf("%s | %s nuclei", GROUP, format(length(bc), big.mark = ",")))
  ggsave(file.path(OUTDIR, "decon_umap_reference.png"), ref,
         width = 11, height = 10, dpi = 200)
  rm(obj); invisible(gc())

  panel_theme <- theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey35"),
          axis.title = element_blank(), axis.text = element_text(size = 6),
          legend.key.width = grid::unit(0.25, "cm"))

  # One colour scale per gene, cut at the p99 of the PRE non-zero values, so the
  # two panels are comparable and the pre panel is not blown out by its own tail.
  fp <- function(v, ttl, sub, cut) {
    d <- data.table(x = emb$UMAP_1, y = emb$UMAP_2, v = pmin(v, cut))
    setorder(d, v)                                   # expressing nuclei on top
    ggplot(d, aes(x, y, colour = v)) +
      geom_scattermore(pointsize = PT_PX, pixels = RAST_PX) +
      scale_colour_gradient(low = COL_LO, high = COL_HI, limits = c(0, cut),
                            name = "log-norm") +
      coord_equal() + labs(title = ttl, subtitle = sub) + panel_theme
  }

  cuts <- vapply(genes, function(g) {
    e <- pre_n[g, ]; e <- e[e > 0]
    if (!length(e)) CUT_FLOOR else max(as.numeric(quantile(e, Q_HI)), CUT_FLOOR)
  }, numeric(1))
  fwrite(data.table(gene = genes, max_cutoff = round(cuts, 3), quantile = Q_HI),
         file.path(OUTDIR, "decon_scale_cutoffs.csv"))

  rows <- list()
  for (g in genes) {
    st <- gstats[gene == g]
    a <- fp(pre_n[g, ], paste0(g, "  -  pre-CellBender"),
            sprintf("%s%% of nuclei detected | %s UMIs", st$pct_detected_pre,
                    format(st$umi_pre, big.mark = ",")), cuts[[g]])
    b <- fp(post_n[g, ], paste0(g, "  -  post-CellBender"),
            sprintf("%s%% of nuclei detected | %s UMIs (-%s%%)", st$pct_detected_post,
                    format(st$umi_post, big.mark = ","), st$pct_umi_removed), cuts[[g]])
    row <- a | b
    rows[[g]] <- row
    ggsave(file.path(OUTDIR, sprintf("decon_featureplot_%s.png", gsub("[^A-Za-z0-9]+", "_", g))),
           row + plot_annotation(
             title = g,
             subtitle = sprintf("shared colour scale, cut at the p99 of the pre non-zero values (%.2f); off-%s detection %s%% -> %s%%",
                                cuts[[g]], st$top_group,
                                st$pct_detected_offtarget_pre, st$pct_detected_offtarget_post),
             theme = theme(plot.title = element_text(size = 16, face = "bold"),
                           plot.subtitle = element_text(size = 9, colour = "grey35"))),
           width = 11, height = 5.6, dpi = 300)
    msg("written decon_featureplot_", g, ".png")
  }

  all_p <- wrap_plots(rows, ncol = 1) +
    plot_annotation(
      title = "Ambient RNA removal (CellBender), before and after",
      subtitle = paste0("log-normalised counts on the same UMAP, each side normalised on its own ",
                        "mouse-gene library size; one colour scale per gene"),
      caption = "cellbender remove-background, fpr = 0; transgene features excluded",
      theme = theme(plot.title = element_text(size = 18, face = "bold"),
                    plot.subtitle = element_text(size = 10, colour = "grey35")))
  ggsave(file.path(OUTDIR, "decon_featureplots_all.png"), all_p,
         width = 11, height = 5.2 * length(genes) + 1, dpi = 200, limitsize = FALSE)
  msg("written decon_featureplots_all.png (", length(genes), " genes)")
}

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, paste0("sessionInfo_reporting_", SECTION, ".txt")))
msg("done: ", OUTDIR)
