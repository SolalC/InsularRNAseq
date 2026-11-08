# Plots and pathway enrichment of the per-cell DEG (percell_deg_diet_sexfollowup.R),
# HFD vs SD, per `annot`.
#
# Hit:     p_bonferroni < 0.05 and |avg_log2FC| > 0.5 in DEG_percell_HFDvsSD_<celltype>.csv,
#          used for every plot and for the enrichment.
# Genes:   heatmaps show the HEAT_MAX hits with the lowest p_val.
# Plot 1   animals: as in pseudobulk_deg_seurat_diet.R. Counts summed per mouse
#          (AggregateExpression), log2(CPM + 1), z-scored per gene across the mice; an
#          extra column at the right of each diet shows that diet's mean z.
# Plot 2   cells:   all nuclei, Seurat DoHeatmap; ScaleData z-scores of RNA `data` across the
#          cell type's nuclei, colour clipped at +/- Z_CLIP; columns grouped by mouse.
# Both:    columns grouped by diet then mouse; rows in the same order (clustered on the
#          animal matrix), so the two plots can be read side by side.
# Plot 3   volcano: every tested gene, avg_log2FC vs -log10 per-cell p_val. Hits
#          coloured; ringed if also pseudobulk FDR < 0.05.
#          Wilcoxon p-values underflow to 0 for large cell types; these are drawn at
#          the smallest non-zero p as triangles.
# Output:  in DEG_DIR, one page per cell type in
#          heatmap_percell_animals.pdf and heatmap_percell_cells.pdf, and a single
#          volcano_percell_all.jpeg faceted by cell type.
# Plot 4   sex: F vs M log2FC of each pooled hit (sexsplit tables), one facet per cell type,
#          y = x drawn; points off the diagonal respond differently by sex.
# Plot 5   sex: interaction plots of the SEX_TOP hits per cell type with the largest
#          |F - M| log2FC: per-mouse means by diet, one line per sex.
#          Output: sex_log2FC_F_vs_M.jpeg, sex_interaction_tophits.pdf (one page per type).
# Then     enrichment: gprofiler2::gost, KEGG and WikiPathways (section at the end).
# Usage:   Rscript percell_deg_plots.R [all | excitatory | inhibitory]
#          Same argument as percell_deg_diet_sexfollowup.R; reads and writes that run's
#          subfolder (all_celltypes, excitatory_subclasses, inhibitory_subclasses).
# Untested. Needs the memory of the DEG run; submit with sbatch, not on the login node.

suppressPackageStartupMessages({
  library(Seurat); library(data.table); library(ggplot2); library(gprofiler2)
})

ROOT            <- "/QRISdata/Q9732/scratch"
PROJECT         <- "/scratch/user/uqschauq/Insular"   # where scripts 9 and 10_x write
args            <- commandArgs(trailingOnly = TRUE)
RUN             <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "all"
LINEAGE_CLASS   <- c(excitatory = "Glutamatergic", inhibitory = "GABAergic")
if (!RUN %in% c("all", names(LINEAGE_CLASS)))
  stop("argument must be one of: all, ", paste(names(LINEAGE_CLASS), collapse = ", "))
AZIMUTH_RES     <- 0.2
SUBDIR          <- if (RUN == "all") "all_celltypes" else paste0(RUN, "_subclasses")
first_existing  <- function(x) c(Filter(file.exists, x), x[1])[1]
IN              <- if (RUN == "all") file.path(ROOT, "Data", "Seurat", "HFDseq_v2_annotated.rds") else
                     file.path(ROOT, "Data", "Seurat", sprintf("HFDseq_v2_sub_%s.rds", RUN))
PRED_FILE       <- if (RUN == "all") NA_character_ else
  first_existing(file.path(c(ROOT, PROJECT), "Results", sprintf("10_2_%s_azimuth", RUN),
                           sprintf("azimuth_predictions_res%s.rds", AZIMUTH_RES)))
DEG_DIR         <- file.path(ROOT, "Results", "percell_wilcox_diet", SUBDIR)
TRANSGENES      <- c("EGFP", "tdTomato", "EBFP2", "AAV-shared")
BONF_CUT        <- 0.05
LFC_CUT         <- 0.5   # |avg_log2FC| a hit must exceed (log2, so 1.41-fold)
HEAT_MAX        <- 50    # at most this many genes, so labels stay readable
Z_CLIP          <- 2.5   # colour limit for the cell heatmap (DoHeatmap's default); single-nucleus z is heavy-tailed
VOLC_LABEL      <- 8     # hits labelled per cell type facet (lowest p_val)
SEX_TOP         <- 9     # genes per cell type in the sex interaction plots
SEX_LABEL       <- 6     # genes labelled per facet in the F vs M scatter
SEED            <- 42
DROP_TYPES  <- "Unidentified"   # `annot` levels left out of every analysis

safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
msg  <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
is_hit <- function(d) d$p_bonferroni < BONF_CUT & abs(d$avg_log2FC) > LFC_CUT
set.seed(SEED)

## ---- object ------------------------------------------------------------------------
msg("run = ", RUN, "; reading ", basename(IN))
obj <- readRDS(IN)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]], search = "counts")) > 1) obj <- JoinLayers(obj)
obj <- DietSeurat(obj, assays = "RNA", layers = c("counts", "data"), dimreducs = NULL, graphs = NULL)
invisible(gc())
if (RUN != "all") {   # same subclass labels and filter as the DEG script
  pred <- as.data.table(readRDS(PRED_FILE))
  i <- match(colnames(obj), pred$barcode)
  stopifnot(!anyNA(i))
  obj$annot <- as.character(pred$predicted.subclass[i])
  pc    <- as.character(pred$predicted.class[i])
  modal <- data.table(annot = obj$annot, pc = pc)[, .N, by = .(annot, pc)][order(-N)][, .SD[1], by = annot]
  obj   <- subset(obj, cells = colnames(obj)[obj$annot %in% modal[pc == LINEAGE_CLASS[[RUN]], annot]])
  rm(pred); invisible(gc())
}
obj$annot  <- as.character(obj$annot)
ANNOT_NAMES <- setNames(unique(obj$annot), safe(unique(obj$annot)))   # file-name form -> label
obj$sample <- as.character(obj$sample)
mice <- unique(as.data.table(obj[[]])[, .(sample, diet = factor(as.character(diet), levels = c("SD", "HFD")),
                                          sex = factor(as.character(sex), levels = c("F", "M")))])
setorder(mice, diet, sex, sample)
MOUSE_ORDER <- mice$sample
N_BONF      <- nrow(obj)   # Seurat's Bonferroni multiplier (all genes in the object)

# Pseudobulk once for all cell types. AggregateExpression turns "_" into "-" in group names.
pb <- AggregateExpression(obj, assays = "RNA", group.by = c("annot", "sample"), return.seurat = TRUE)
pb_cnt <- LayerData(pb, assay = "RNA", layer = "counts")
pb_cnt <- pb_cnt[setdiff(rownames(pb_cnt), TRANSGENES), ]
pb_md  <- data.table(library = colnames(pb), annot = pb$annot, sample = pb$sample)

## ---- shared plot pieces --------------------------------------------------------------
heat_theme <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), panel.spacing.x = unit(1, "pt"),
        strip.text.x = element_text(angle = 90, hjust = 0))

# ggrepel if installed, plain text otherwise
label_layer <- function(data) {
  if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(data = data, aes(label = gene), size = 2.3, max.overlaps = Inf,
                             min.segment.length = 0, segment.size = 0.2)
  else geom_text(data = data, aes(label = gene), size = 2.3, vjust = -0.6)
}

# Volcano data for one cell type; plotted together after the loop.
# p = 0 (Wilcoxon underflow) is floored at that cell type's smallest non-zero p.
volcano_data <- function(all, ct) {
  v <- all[, .(gene, avg_log2FC, p_val, p_bonferroni,
               pb_FDR = if ("pb_FDR" %in% names(all)) pb_FDR else NA_real_)]
  pmin_nz <- min(v[p_val > 0, p_val], 1e-300)
  v[, `:=`(annot = ct, capped = p_val == 0, y = -log10(pmax(p_val, pmin_nz)))]
  v[, cls := fifelse(is_hit(v) & avg_log2FC > 0, "up in HFD",
             fifelse(is_hit(v) & avg_log2FC < 0, "down in HFD", "n.s."))]
  v
}

row_labels <- function(r) setNames(sprintf("%s (p_bonf = %.1e)", r$gene, r$p_bonferroni), r$gene)

## ---- per cell type -----------------------------------------------------------------
volc <- list(); heat_animals <- list(); heat_cells <- list()
for (ct in setdiff(sort(unique(obj$annot)), DROP_TYPES)) {
  f <- file.path(DEG_DIR, paste0("DEG_percell_HFDvsSD_", safe(ct), ".csv"))
  if (!file.exists(f)) { msg(ct, ": no DEG table, skipped"); next }
  all <- fread(f)
  volc[[ct]] <- volcano_data(all, ct)
  r <- all[is_hit(all)][order(p_val)]
  r <- head(r[gene %in% rownames(pb_cnt)], HEAT_MAX)
  if (nrow(r) < 2) { msg(ct, ": ", nrow(r), " hit(s), skipped"); next }
  g <- r$gene
  msg(ct, ": ", length(g), " genes")

  ## animals: log2 CPM of the per-mouse pseudobulk, z-scored per gene
  lib <- pb_md[annot == gsub("_", "-", ct)]
  lib[, sample := gsub("-", "_", sample)]
  cnt <- as.matrix(pb_cnt[, lib$library, drop = FALSE])
  za  <- log2(t(t(cnt) / colSums(cnt)) * 1e6 + 1)[g, , drop = FALSE]
  colnames(za) <- lib$sample
  za  <- t(scale(t(za)))
  za[is.na(za)] <- 0                                   # gene constant across mice
  gene_order <- rownames(za)[hclust(dist(za))$order]

  d <- as.data.table(as.table(za))
  setnames(d, c("gene", "sample", "z"))
  d <- merge(d, mice, by = "sample")
  # Group average: mean z of the diet's mice, one extra column at the right of each diet
  avg <- d[, .(z = mean(z)), by = .(gene, diet)][, sample := paste(diet, "average")]
  d <- rbind(d[, .(gene, sample, z, diet)], avg[, .(gene, sample, z, diet)])
  col_levels <- unlist(lapply(levels(mice$diet), function(dt)
    c(mice[diet == dt, sample], paste(dt, "average"))))
  d[, `:=`(sample = factor(sample, levels = col_levels), gene = factor(gene, levels = gene_order))]
  sep <- mice[, .(x = .N + 0.5), by = diet]             # white line before each average column
  lim <- max(abs(d$z))
  gp <- ggplot(d, aes(sample, gene, fill = z)) +
    geom_tile() +
    geom_vline(data = sep, aes(xintercept = x), colour = "white", linewidth = 1.2) +
    facet_grid(~ diet, scales = "free_x", space = "free_x") +
    scale_y_discrete(labels = row_labels(r)) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-lim, lim), name = "z-score\nlog2 CPM") +
    labs(x = NULL, y = NULL,
         title = paste0(ct, ": HFD vs SD, per mouse, ", length(g), " per-cell hits")) +
    heat_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                       strip.text.x = element_text(angle = 0))
  heat_animals[[ct]] <- gp

  ## cells: every nucleus of the cell type, Seurat's DoHeatmap on ScaleData z-scores
  # (scaled across this cell type's nuclei). Groups = mice, ordered by diet then sex.
  sc <- subset(obj, cells = colnames(obj)[obj$annot == ct], features = g)
  sc$mouse <- factor(sc$sample, levels = MOUSE_ORDER)
  sc <- ScaleData(sc, features = g, verbose = FALSE)
  gp <- DoHeatmap(sc, features = rev(gene_order), group.by = "mouse", slot = "scale.data",
                  disp.min = -Z_CLIP, disp.max = Z_CLIP, size = 2.5, angle = 90,
                  raster = TRUE) +
    scale_y_discrete(labels = row_labels(r)) +
    labs(title = paste0(ct, ": HFD vs SD, all ", ncol(sc), " nuclei, ", length(g), " per-cell hits"),
         subtitle = "columns grouped by mouse: SD (F, M) then HFD (F, M)") +
    theme(axis.text.y = element_text(size = 6), legend.position = "right")
  rm(sc)
  heat_cells[[ct]] <- gp
}

## ---- heatmaps: one multi-page PDF per type, one page per cell type -------------------
# Page height fits HEAT_MAX genes; cell types with fewer genes leave white space below.
write_pdf <- function(plots, file, width) {
  if (!length(plots)) { msg("no heatmap for ", basename(file)); return(invisible()) }
  pdf(file, width = width, height = 1.5 + 0.14 * HEAT_MAX + 1, onefile = TRUE)
  for (gp in plots) print(gp)
  dev.off()
}
write_pdf(heat_animals, file.path(DEG_DIR, "heatmap_percell_animals.pdf"), width = 3 + 0.3 * length(MOUSE_ORDER) + 2)
write_pdf(heat_cells,   file.path(DEG_DIR, "heatmap_percell_cells.pdf"),   width = 12)

## ---- volcano: all cell types, one facet each ------------------------------------------
if (!length(volc)) stop("no DEG_percell_HFDvsSD_*.csv in ", DEG_DIR, " -- run percell_deg_diet_sexfollowup.sh ", RUN, " first")
v <- rbindlist(volc)
v[, cls := factor(cls, levels = c("up in HFD", "down in HFD", "n.s."))]
setorder(v, -cls, -p_val)                                # n.s. drawn first, hits on top
n_hit <- v[, .(lab = sprintf("%s\n%d significant genes", annot[1], sum(cls != "n.s."))), by = annot]
v[, facet := factor(n_hit$lab[match(annot, n_hit$annot)], levels = n_hit[order(annot), lab])]
lab <- v[cls != "n.s."][order(p_val), head(.SD, VOLC_LABEL), by = annot]
gp <- ggplot(v, aes(avg_log2FC, y, colour = cls)) +
  geom_point(aes(shape = capped), size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(BONF_CUT / N_BONF), linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ facet, scales = "free", ncol = 4) +
  scale_colour_manual(values = c("up in HFD" = "#b2182b", "down in HFD" = "#2166ac", "n.s." = "grey75"),
                      drop = FALSE, name = paste0("p_bonferroni < ", BONF_CUT, ", |log2FC| > ", LFC_CUT)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none")
if (v[pb_FDR < 0.05, .N])
  gp <- gp + geom_point(data = v[pb_FDR < 0.05], shape = 1, size = 1.5, colour = "black", stroke = 0.3)
gp <- gp + label_layer(lab) +
  labs(x = "avg_log2FC (HFD / SD)", y = "-log10 per-cell p (Wilcoxon)",
       title = "HFD vs SD, per nucleus, by cell type",
       subtitle = paste0("dashed: Bonferroni 0.05 and |log2FC| = ", LFC_CUT, "; ringed: pseudobulk FDR < 0.05; ",
                         "triangles: p = 0, drawn at the cell type's smallest non-zero p")) +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(DEG_DIR, "volcano_percell_all.jpeg"), gp,
       width = 14, height = 1 + 3.2 * ceiling(uniqueN(v$annot) / 4), dpi = 300, bg = "white",
       limitsize = FALSE)

rm(obj, pb, pb_cnt); invisible(gc())

## ---- sex: F vs M diet effect -----------------------------------------------------------
# Reads DEG_percell_sexsplit_<ct>.csv and mouse_means_hits_<ct>.csv from the DEG script.
# These compare effect sizes; they are not a test of the sex difference. The only test
# here is the pseudobulk interaction, present if 9_pseudobulk_deg.R was run: its BH is
# re-applied over the cell type's hits only (pbB_int_FDR_hits, computed by the DEG script).
sx <- list(); mm_all <- list()
for (ct in setdiff(sort(unique(names(volc))), DROP_TYPES)) {
  f  <- file.path(DEG_DIR, paste0("DEG_percell_sexsplit_", safe(ct), ".csv"))
  fm <- file.path(DEG_DIR, paste0("mouse_means_hits_", safe(ct), ".csv"))
  if (!file.exists(f)) next
  d <- fread(f)
  if (!all(c("F_avg_log2FC", "M_avg_log2FC") %in% names(d))) next
  d <- d[is_hit(d) & !is.na(F_avg_log2FC) & !is.na(M_avg_log2FC)]
  if (!nrow(d)) next
  d[, `:=`(annot = ct, sex_diff = M_avg_log2FC - F_avg_log2FC,
           int = if ("pbB_int_FDR_hits" %in% names(d)) pbB_int_FDR_hits else NA_real_)]
  sx[[ct]] <- d[, .(annot, gene, F_avg_log2FC, M_avg_log2FC, sex_diff, int)]
  if (file.exists(fm)) mm_all[[ct]] <- fread(fm)[gene %in% d$gene][, annot := ct]
}

if (length(sx)) {
  sx <- rbindlist(sx)
  strip <- sx[, .(lab = sprintf("%s\nn = %d, r = %.2f", annot[1], .N,
                                if (.N > 2) cor(F_avg_log2FC, M_avg_log2FC) else NA_real_)), by = annot]
  sx[, facet := factor(strip$lab[match(annot, strip$annot)], levels = strip[order(annot), lab])]
  lab <- sx[order(-abs(sex_diff)), head(.SD, SEX_LABEL), by = annot]
  lim <- max(abs(c(sx$F_avg_log2FC, sx$M_avg_log2FC)))
  gp <- ggplot(sx, aes(F_avg_log2FC, M_avg_log2FC)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(size = 0.8, alpha = 0.7, colour = "black") +
    label_layer(lab) +
    facet_wrap(~ facet, ncol = 4) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = "log2FC HFD / SD, females", y = "log2FC HFD / SD, males",
         title = "Diet effect by sex, per-cell hits (sexes pooled for selection)",
         subtitle = paste0("dashed: same effect in both sexes; labelled: largest |M - F|; ",
                           "r = Pearson correlation of F and M log2FC")) +
    theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank())
  ggsave(file.path(DEG_DIR, "sex_log2FC_F_vs_M.jpeg"), gp,
         width = 13, height = 1.2 + 3.4 * ceiling(uniqueN(sx$annot) / 4), dpi = 300, bg = "white",
         limitsize = FALSE)

  ## interaction plots, one PDF page per cell type
  if (length(mm_all)) {
    mm_all <- rbindlist(mm_all, fill = TRUE)
    mm_all[, `:=`(diet = factor(diet, levels = c("SD", "HFD")), sex = factor(sex, levels = c("F", "M")))]
    pdf(file.path(DEG_DIR, "sex_interaction_tophits.pdf"), width = 9, height = 7.5, onefile = TRUE)
    for (ct in sort(unique(mm_all$annot))) {
      top <- sx[annot == ct][order(-abs(sex_diff)), head(gene, SEX_TOP)]
      d   <- mm_all[annot == ct & gene %in% top]
      if (!nrow(d)) next
      info <- sx[annot == ct & gene %in% top]
      info[, lab := sprintf("%s\nF %+.2f, M %+.2f%s", gene, F_avg_log2FC, M_avg_log2FC,
                            fifelse(is.na(int), "", sprintf(", int FDR %.2g", int)))]
      d[, gene := factor(info$lab[match(gene, info$gene)], levels = info$lab[match(top, info$gene)])]
      means <- d[, .(mean_expr = mean(mean_expr)), by = .(gene, diet, sex)]
      gp <- ggplot(d, aes(diet, mean_expr, colour = sex, group = sex)) +
        geom_point(position = position_dodge(width = 0.3), size = 1.8) +
        geom_line(data = means, position = position_dodge(width = 0.3), linewidth = 0.6) +
        facet_wrap(~ gene, scales = "free_y", ncol = 3) +
        scale_colour_manual(values = c(F = "#d95f02", M = "#1b9e77")) +
        labs(x = NULL, y = "mean expression per mouse (linear, normalised)",
             title = paste0(ct, ": top ", length(top), " hits by |M - F| log2FC"),
             subtitle = "points: mice (n = 3 per sex x diet); lines: group means; F/M: per-cell log2FC within sex") +
        theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank())
      print(gp)
    }
    dev.off()
  }
} else msg("no sexsplit tables with hits; sex plots skipped")

## ---- pathway enrichment -----------------------------------------------------------
# Query:      per cell type, all hits (up and down together).
# Background: the genes FindMarkers tested in that cell type (domain_scope = "custom").
#             The default, all annotated mouse genes, would make any pathway of
#             brain-expressed genes look enriched. These genes already passed
#             FindMarkers' |log2FC| >= LOGFC_MIN prefilter in the DEG script.
# p-value:    gost's p_value, already corrected (g:SCS, gprofiler's default).
# Plots:      one bar plot per source, x = -log10 p, y = pathway, facet by cell type;
#             TOP_N terms per cell type.
# Output:     in DEG_DIR/enrichment: enrichment_percell_all.csv,
#             enrichment_percell_KEGG.jpeg, enrichment_percell_WP.jpeg.
# gost queries the g:Profiler web server, so the node needs internet access. A failed
# query is logged and skipped; the plots above are already saved by then.

OUTDIR   <- file.path(DEG_DIR, "enrichment")
SOURCES  <- c(KEGG = "KEGG", WP = "WikiPathways")
MIN_HITS <- 5     # a cell type needs this many hits to be queried
TOP_N    <- 15    # terms shown per cell type
MAX_TERM <- 200   # drop very large terms (often generic: "Metabolic pathways")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(DEG_DIR, pattern = "^DEG_percell_HFDvsSD_.*\\.csv$", full.names = TRUE)
files <- files[!basename(files) %in% paste0("DEG_percell_HFDvsSD_", safe(DROP_TYPES), ".csv")]
if (!length(files)) stop("no DEG_percell_HFDvsSD_*.csv in ", DEG_DIR)

res <- list()
for (f in files) {
  ct   <- sub("^DEG_percell_HFDvsSD_(.*)\\.csv$", "\\1", basename(f))   # safe() form of annot
  d    <- fread(f)
  hits <- d[is_hit(d), gene]
  msg(ct, ": ", nrow(d), " tested, ", length(hits), " hits")
  if (length(hits) < MIN_HITS) next
  g <- tryCatch(gost(query = hits, organism = "mmusculus", sources = names(SOURCES),
                      custom_bg = d$gene, domain_scope = "custom", user_threshold = 0.05,
                      correction_method = "g_SCS", significant = TRUE, evcodes = TRUE),
                 error = function(err) { msg("  gost failed: ", conditionMessage(err)); NULL })
  if (is.null(g)) { msg("  no term at p < 0.05"); next }
  res[[ct]] <- as.data.table(g$result)[, .(annot = ct, source, term_id, term_name,
                                           p_value, term_size, query_size, intersection_size,
                                           precision, recall, intersection)]
}
if (!length(res)) { msg("no enriched term in any cell type"); quit(save = "no") }
res <- rbindlist(res)
res[, annot := fifelse(annot %in% names(ANNOT_NAMES), ANNOT_NAMES[annot], gsub("_", " ", annot))]
setorder(res, source, annot, p_value)
fwrite(res, file.path(OUTDIR, "enrichment_percell_all.csv"))

## ---- bar plots, one per source ------------------------------------------------------
for (src in names(SOURCES)) {
  d <- res[source == src & term_size <= MAX_TERM]
  if (!nrow(d)) { msg(src, ": nothing to plot"); next }
  d <- d[order(p_value), head(.SD, TOP_N), by = annot]
  d[, mlp := -log10(p_value)]
  # y key unique per facet, ordered by p within the facet; label strips the key
  d[, key := paste(term_name, annot, sep = "___")]
  d[, key := factor(key, levels = d[order(annot, mlp), key])]
  n_rows <- d[, .N, by = annot][, max(N)]
  n_ct   <- uniqueN(d$annot)
  gp <- ggplot(d, aes(mlp, key)) +
    geom_col(fill = "grey40", width = 0.8) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~ annot, scales = "free_y", ncol = 3) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    labs(x = "-log10 adjusted p (g:SCS)", y = NULL,
         title = paste0(SOURCES[[src]], " enrichment of per-cell HFD vs SD hits (p_bonferroni < ",
                        BONF_CUT, ", |log2FC| > ", LFC_CUT, ")"),
         subtitle = paste0("background: genes tested per cell type; top ", TOP_N,
                           " terms per cell type; terms > ", MAX_TERM, " genes dropped")) +
    theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
  ggsave(file.path(OUTDIR, paste0("enrichment_percell_", src, ".jpeg")), gp,
         width = 15, height = 1.5 + ceiling(n_ct / 3) * (0.6 + 0.16 * n_rows),
         dpi = 300, bg = "white", limitsize = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done")
