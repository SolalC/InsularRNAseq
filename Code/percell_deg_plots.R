# Plots and pathway enrichment of the per-cell DEG (percell_deg_diet_sexfollowup.R),
# HFD vs SD, per `annot`.
#
# Hit:     p_bonferroni < 0.05 and |avg_log2FC| > 0.5 in DEG_percell_HFDvsSD_<celltype>.csv,
#          used for every plot and for the enrichment.
# Genes:   heatmaps show the HEAT_MAX hits with the lowest p_val.
# Plot 1   animals: counts summed per mouse (AggregateExpression), log2(CPM + 1), then
#          z-scored per gene against the SD (standard diet) mice:
#            z = (x - mean of SD mice) / sd of SD mice
#          so SD mice scatter around 0 with unit spread and an HFD mouse at z = 2 lies two
#          SD-group standard deviations from the SD mean. Colour clipped at +/- ZSD_CLIP.
#          An extra column at the right of each diet shows that diet's mean.
# Plot 1b  animals, log2FC copy (sanity check of plot 1): same matrix, numerator only,
#            log2FC = x - mean of SD mice   (log2 CPM, so a log2 fold change vs the SD mean)
#          Same rows, same order. Plots 1 and 1b differ per gene only by the factor
#          1 / sd(SD mice), so every row has the same sign pattern in both. Row labels give
#          the per-cell avg_log2FC; the HFD average column should agree with it in sign and
#          roughly in size (not exactly: avg_log2FC is a log of means over nuclei, this is a
#          mean of logs over mice).
# Plot 2   cells:   all nuclei, Seurat DoHeatmap; RNA `data` z-scored per gene against the
#          cell type's SD nuclei (mean and sd of SD nuclei), colour clipped at +/- Z_CLIP;
#          columns grouped by mouse.
# Both:    columns grouped by diet then mouse; rows in the same order (clustered on the
#          animal matrix), so the two plots can be read side by side.
# Plot 3   volcano: every tested gene, avg_log2FC vs -log10 per-cell p_val. Hits
#          coloured; ringed if also pseudobulk FDR < 0.05.
#          Wilcoxon p-values underflow to 0 for large cell types; these are drawn at
#          the smallest non-zero p as triangles.
# Output:  in DEG_DIR, one page per cell type in heatmap_percell_animals_zSD.pdf,
#          heatmap_percell_animals_log2FC.pdf and heatmap_percell_cells_zSD.pdf, and a single
#          volcano_percell_all.pdf faceted by cell type.
# Plot 4   sex: F vs M log2FC of each pooled hit (sexsplit tables), one facet per cell type,
#          y = x drawn; points off the diagonal respond differently by sex.
# Plot 5   sex: interaction plots of the SEX_TOP hits per cell type with the largest
#          |F - M| log2FC: per-mouse means by diet, one line per sex.
#          Output: sex_log2FC_F_vs_M.pdf, sex_interaction_tophits.pdf (one page per type).
# Plot 6   sex: Venn diagram per cell type of female vs male DEG. Sex DEG = within-sex
#          F_/M_p_bonferroni < BONF_CUT and |F_/M_avg_log2FC| > LFC_CUT in the sexsplit
#          table. Caveats, all from the DEG script's design:
#          (a) only pooled hits were tested within sex, so the universe is the pooled hits;
#              a gene changed in one sex only but missed by the pooled test is absent;
#          (b) within-sex Bonferroni is over the cell type's pooled hits, not all genes;
#          (c) each sex has half the nuclei, so lower power: a gene in one circle only is
#              not evidence of a sex-specific effect (see plot 4 and the interaction FDR).
#          Output: sex_venn_F_vs_M.pdf (one facet per cell type), sex_venn_membership.csv.
# Plot 7   main heatmap per cell type: SD mice | per-cell log2FC | HFD mice (section before
#          the enrichment; standalone copy in percell_heatmap_main.R).
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
ZSD_CLIP        <- 4     # colour limit for the animal z vs SD heatmap; sd of 6 SD mice can be small, so HFD z can be large
VOLC_LABEL      <- 8     # hits labelled per cell type facet (lowest p_val)
SEX_TOP         <- 9     # genes per cell type in the sex interaction plots
SEX_LABEL       <- 6     # genes labelled per facet in the F vs M scatter
SEED            <- 42
MAKE_VOLCANO    <- FALSE # TRUE to write volcano_percell_all.pdf
DROP_TYPES  <- "Unidentified"   # `annot` levels left out of every analysis

# All single-page figures are written as PDF: on Bunya compute nodes the JPEG device hung
# (no CPU use) at ggsave, while pdf() works.
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

# z-score of each row (gene) against the reference columns `ref` (SD mice or SD nuclei):
# (x - mean(x[ref])) / sd(x[ref]). A gene with sd = 0 in the reference (e.g. undetected in
# every SD sample) has no defined z; it is set to 0 where x equals the reference mean and to
# +/- clip elsewhere, so it shows as saturated rather than disappearing.
ref_z <- function(m, ref, clip) {
  mu <- rowMeans(m[, ref, drop = FALSE])
  s  <- apply(m[, ref, drop = FALSE], 1, sd)
  dv <- m - mu
  z  <- dv / s
  bad <- s == 0 | is.na(s)
  if (any(bad)) z[bad, ] <- sign(dv[bad, , drop = FALSE]) * clip
  z
}

# Animal heatmap from a gene x mouse matrix: diet facets, one "average" column per diet
# (mean of that diet's mice), rows in gene_order, symmetric colour scale clipped at +/- lim.
animal_heatmap <- function(m, gene_order, labels, lim, fill_name, title) {
  d <- as.data.table(as.table(m))
  setnames(d, c("gene", "sample", "value"))
  d <- merge(d, mice, by = "sample")
  avg <- d[, .(value = mean(value)), by = .(gene, diet)][, sample := paste(diet, "average")]
  d <- rbind(d[, .(gene, sample, value, diet)], avg[, .(gene, sample, value, diet)])
  col_levels <- unlist(lapply(levels(mice$diet), function(dt)
    c(mice[diet == dt, sample], paste(dt, "average"))))
  d[, `:=`(sample = factor(sample, levels = col_levels), gene = factor(gene, levels = gene_order))]
  sep <- mice[, .(x = .N + 0.5), by = diet]             # white line before each average column
  ggplot(d, aes(sample, gene, fill = value)) +
    geom_tile() +
    geom_vline(data = sep, aes(xintercept = x), colour = "white", linewidth = 1.2) +
    facet_grid(~ diet, scales = "free_x", space = "free_x") +
    scale_y_discrete(labels = labels) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-lim, lim), oob = scales::squish, name = fill_name) +
    labs(x = NULL, y = NULL, title = title) +
    heat_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
                       strip.text.x = element_text(angle = 0))
}

## ---- per cell type -----------------------------------------------------------------
volc <- list(); heat_animals <- list(); heat_animals_lfc <- list(); heat_cells <- list()
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

  ## animals: log2 CPM of the per-mouse pseudobulk, referenced to the SD mice
  lib <- pb_md[annot == gsub("_", "-", ct)]
  lib[, sample := gsub("-", "_", sample)]
  cnt <- as.matrix(pb_cnt[, lib$library, drop = FALSE])
  lc  <- log2(t(t(cnt) / colSums(cnt)) * 1e6 + 1)[g, , drop = FALSE]
  colnames(lc) <- lib$sample
  is_sd <- colnames(lc) %in% mice[diet == "SD", sample]
  stopifnot(sum(is_sd) >= 2, sum(!is_sd) >= 1)
  zsd <- ref_z(lc, is_sd, ZSD_CLIP)                    # z vs SD mice (plot 1)
  lfc <- lc - rowMeans(lc[, is_sd, drop = FALSE])      # log2FC vs SD mean (plot 1b)
  gene_order <- rownames(zsd)[hclust(dist(pmin(pmax(zsd, -ZSD_CLIP), ZSD_CLIP)))$order]

  heat_animals[[ct]] <- animal_heatmap(
    zsd, gene_order, row_labels(r), lim = ZSD_CLIP,
    fill_name = paste0("z vs SD mice\n(x - mean SD) / sd SD\nlog2 CPM,\nclipped at +/-", ZSD_CLIP),
    title = paste0(ct, ": HFD vs SD, per mouse, ", length(g), " per-cell hits, z-score vs SD"))
  lfc_lim <- max(abs(lfc), 0.1)
  heat_animals_lfc[[ct]] <- animal_heatmap(
    lfc, gene_order,
    setNames(sprintf("%s (per-cell log2FC = %+.2f)", r$gene, r$avg_log2FC), r$gene), lim = lfc_lim,
    fill_name = "log2FC vs SD mean\nlog2 CPM",
    title = paste0(ct, ": HFD vs SD, per mouse, ", length(g), " per-cell hits, log2FC vs SD (check of z plot)"))

  ## cells: every nucleus of the cell type, DoHeatmap on RNA `data` z-scored against the
  # cell type's SD nuclei (written to scale.data). Groups = mice, ordered by diet then sex.
  sc <- subset(obj, cells = colnames(obj)[obj$annot == ct], features = g)
  sc$mouse <- factor(sc$sample, levels = MOUSE_ORDER)
  dat <- as.matrix(LayerData(sc, assay = "RNA", layer = "data")[g, , drop = FALSE])
  LayerData(sc, assay = "RNA", layer = "scale.data") <- ref_z(dat, as.character(sc$diet) == "SD", Z_CLIP)
  rm(dat)
  gp <- DoHeatmap(sc, features = rev(gene_order), group.by = "mouse", slot = "scale.data",
                  disp.min = -Z_CLIP, disp.max = Z_CLIP, size = 2.5, angle = 90,
                  raster = TRUE) +
    scale_y_discrete(labels = row_labels(r)) +
    labs(title = paste0(ct, ": HFD vs SD, all ", ncol(sc), " nuclei, ", length(g), " per-cell hits"),
         subtitle = paste0("z = (x - mean of SD nuclei) / sd of SD nuclei, RNA data, clipped at +/-", Z_CLIP,
                           "; columns grouped by mouse: SD (F, M) then HFD (F, M)")) +
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
aw <- 3 + 0.3 * length(MOUSE_ORDER) + 2.5
write_pdf(heat_animals,     file.path(DEG_DIR, "heatmap_percell_animals_zSD.pdf"),    width = aw)
write_pdf(heat_animals_lfc, file.path(DEG_DIR, "heatmap_percell_animals_log2FC.pdf"), width = aw)
write_pdf(heat_cells,       file.path(DEG_DIR, "heatmap_percell_cells_zSD.pdf"),      width = 12)

## ---- volcano: all cell types, one facet each ------------------------------------------
if (!length(volc)) stop("no DEG_percell_HFDvsSD_*.csv in ", DEG_DIR, " -- run percell_deg_diet_sexfollowup.sh ", RUN, " first")
if (MAKE_VOLCANO) {
msg("volcano: building")
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
msg("volcano: writing pdf")
ggsave(file.path(DEG_DIR, "volcano_percell_all.pdf"), gp,
       width = 14, height = 1 + 3.2 * ceiling(uniqueN(v$annot) / 4), bg = "white",
       limitsize = FALSE)
} else msg("volcano: skipped (MAKE_VOLCANO = FALSE)")

rm(obj, pb, pb_cnt); invisible(gc())

## ---- sex: F vs M diet effect -----------------------------------------------------------
# Reads DEG_percell_sexsplit_<ct>.csv and mouse_means_hits_<ct>.csv from the DEG script.
# These compare effect sizes; they are not a test of the sex difference. The only test
# here is the pseudobulk interaction, present if 9_pseudobulk_deg.R was run: its BH is
# re-applied over the cell type's hits only (pbB_int_FDR_hits, computed by the DEG script).
msg("sex: reading sexsplit tables")
sx <- list(); mm_all <- list(); venn <- list()
for (ct in setdiff(sort(unique(names(volc))), DROP_TYPES)) {
  f  <- file.path(DEG_DIR, paste0("DEG_percell_sexsplit_", safe(ct), ".csv"))
  fm <- file.path(DEG_DIR, paste0("mouse_means_hits_", safe(ct), ".csv"))
  if (!file.exists(f)) next
  d <- fread(f)
  if (!all(c("F_avg_log2FC", "M_avg_log2FC") %in% names(d))) next
  if (all(c("F_p_bonferroni", "M_p_bonferroni") %in% names(d)))   # every pooled hit, for the Venn
    venn[[ct]] <- d[is_hit(d), .(annot = ct, gene, F_avg_log2FC, M_avg_log2FC, F_p_bonferroni, M_p_bonferroni)]
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
  msg("sex: writing F vs M scatter")
  ggsave(file.path(DEG_DIR, "sex_log2FC_F_vs_M.pdf"), gp,
         width = 13, height = 1.2 + 3.4 * ceiling(uniqueN(sx$annot) / 4), bg = "white",
         limitsize = FALSE)

  ## interaction plots, one PDF page per cell type
  if (length(mm_all)) {
    mm_all <- rbindlist(mm_all, fill = TRUE)
    mm_all[, `:=`(diet = factor(diet, levels = c("SD", "HFD")), sex = factor(sex, levels = c("F", "M")))]
    msg("sex: writing interaction plots")
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

## ---- sex: Venn of female vs male DEG, one facet per cell type ------------------------
# Universe = the cell type's pooled hits (the only genes tested within sex). Sex DEG uses the
# same cut-offs as the pooled hits: within-sex p_bonferroni < BONF_CUT (Bonferroni over the
# pooled hits) and |within-sex avg_log2FC| > LFC_CUT. A gene not tested in a sex (NA) counts
# as not DEG in that sex. Counts under each region: up / down in HFD; in the overlap, genes
# whose F and M log2FC have opposite signs are counted separately.
if (length(venn)) {
  vn <- rbindlist(venn)
  vn[, `:=`(F_deg = !is.na(F_p_bonferroni) & F_p_bonferroni < BONF_CUT & abs(F_avg_log2FC) > LFC_CUT,
            M_deg = !is.na(M_p_bonferroni) & M_p_bonferroni < BONF_CUT & abs(M_avg_log2FC) > LFC_CUT)]
  vn[, region := fifelse(F_deg & M_deg, "both", fifelse(F_deg, "F only", fifelse(M_deg, "M only", "neither")))]
  fwrite(vn[order(annot, region, gene)], file.path(DEG_DIR, "sex_venn_membership.csv"))

  updown <- function(lfc) sprintf("%d up / %d down", sum(lfc > 0), sum(lfc < 0))
  cnt <- vn[, {
    fo <- region == "F only"; mo <- region == "M only"; bo <- region == "both"
    opp <- bo & sign(F_avg_log2FC) != sign(M_avg_log2FC)
    data.table(x = c(-1.05, 0, 1.05),
               n = c(sum(fo), sum(bo), sum(mo)),
               sub = c(updown(F_avg_log2FC[fo]),
                       paste0(updown(F_avg_log2FC[bo & !opp]), if (any(opp)) sprintf("\n%d opposite sign", sum(opp)) else ""),
                       updown(M_avg_log2FC[mo])),
               n_hit = .N, n_neither = sum(region == "neither"))
  }, by = annot]
  strip <- cnt[, .(lab = sprintf("%s\n%d pooled hits; %d DEG in neither sex", annot[1], n_hit[1], n_neither[1])), by = annot]
  cnt[, facet := factor(strip$lab[match(annot, strip$annot)], levels = strip[order(annot), lab])]

  th <- seq(0, 2 * pi, length.out = 200)
  circ <- rbind(data.table(sex = "F", x = -0.6 + cos(th), y = sin(th)),
                data.table(sex = "M", x =  0.6 + cos(th), y = sin(th)))
  gp <- ggplot() +
    geom_polygon(data = circ, aes(x, y, fill = sex, group = sex), alpha = 0.25, colour = "grey30", linewidth = 0.3) +
    geom_text(data = data.table(x = c(-1.2, 1.2), y = 1.2, lab = c("Female", "Male")),   # repeated in every facet
              aes(x, y, label = lab), size = 3, fontface = "bold") +
    geom_text(data = cnt, aes(x, 0.12, label = n), size = 4.5) +
    geom_text(data = cnt, aes(x, -0.25, label = sub), size = 1.9, lineheight = 0.9) +
    facet_wrap(~ facet, ncol = 4) +
    scale_fill_manual(values = c(F = "#d95f02", M = "#1b9e77"), guide = "none") +
    coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.1, 1.4)) +
    labs(title = "HFD vs SD DEG within each sex, per cell type",
         subtitle = paste0("Sex DEG: within-sex Wilcoxon p_bonferroni < ", BONF_CUT, " (Bonferroni over the cell type's pooled hits) ",
                           "and |log2FC| > ", LFC_CUT, ".\nOnly pooled hits (p_bonferroni < ", BONF_CUT, " over all genes, |log2FC| > ",
                           LFC_CUT, ") were tested within sex. One circle only does not imply a sex-specific effect: ",
                           "each sex has half the nuclei.")) +
    theme_void(base_size = 8) +
    theme(strip.text = element_text(size = 7, margin = margin(2, 0, 2, 0)),
          plot.subtitle = element_text(size = 6.5), plot.margin = margin(6, 6, 6, 6))
  msg("sex: writing Venn")
  ggsave(file.path(DEG_DIR, "sex_venn_F_vs_M.pdf"), gp,
         width = 12, height = 1.2 + 2.6 * ceiling(uniqueN(vn$annot) / 4), bg = "white",
         limitsize = FALSE)
} else msg("no sexsplit tables with F_/M_p_bonferroni; Venn skipped")

## ---- main heatmap: per-mouse expression (SD | log2FC | HFD) ------------------------
# One figure per cell type: left, the SD mice; right, the HFD mice; middle, the per-cell
# DEG avg_log2FC (HFD / SD) of each gene. Rows: the HEAT_MAX hits with the lowest p_val,
# sorted by avg_log2FC (most up at the top).
# Expression: mouse_means_hits_<ct>.csv (DEG script), the mean of expm1(RNA data) over the
#   mouse's nuclei of that cell type, i.e. the same per-mouse average Seurat's avg_log2FC is
#   built from. Shown as log2(mean + 1), then per gene as set by MAIN_SCALE:
#     "zSD"  : (x - mean of SD mice) / sd of SD mice   (SD mice centre on 0)
#     "zall" : (x - mean of all mice) / sd of all mice (conventional row scaling)
#     "log2" : log2(mean + 1) unscaled (colour then mostly shows gene abundance)
# Output: DEG_DIR/heatmap_main_all.pdf, one page per cell type
MAIN_SCALE <- "zSD"
MAIN_CLIP  <- 3      # colour limit for the z panels
main_pages <- list()

main_files <- list.files(DEG_DIR, pattern = "^mouse_means_hits_.*\\.csv$", full.names = TRUE)
main_files <- main_files[!basename(main_files) %in% paste0("mouse_means_hits_", safe(DROP_TYPES), ".csv")]
for (fm in main_files) {
  ct_safe <- sub("^mouse_means_hits_(.*)\\.csv$", "\\1", basename(fm))
  fd <- file.path(DEG_DIR, paste0("DEG_percell_HFDvsSD_", ct_safe, ".csv"))
  if (!file.exists(fd)) next
  deg <- fread(fd)
  hits <- deg[p_bonferroni < BONF_CUT & abs(avg_log2FC) > LFC_CUT][order(p_val)]
  mm   <- fread(fm)
  hits <- head(hits[gene %in% mm$gene], HEAT_MAX)
  if (nrow(hits) < 2) next
  gene_lv <- hits[order(avg_log2FC), gene]               # bottom to top: most up at the top
  mm <- mm[gene %in% hits$gene]
  mm[, `:=`(diet = factor(diet, levels = c("SD", "HFD")), sex = factor(sex, levels = c("F", "M")))]
  mm[, x := log2(mean_expr + 1)]
  mm[, value := switch(MAIN_SCALE,
    zSD  = { s <- sd(x[diet == "SD"]);  if (is.na(s) || s == 0) sign(x - mean(x[diet == "SD"])) * MAIN_CLIP else (x - mean(x[diet == "SD"])) / s },
    zall = { s <- sd(x);                if (is.na(s) || s == 0) 0 else (x - mean(x)) / s },
    log2 = x), by = gene]
  mice_ord <- unique(mm[, .(sample, diet, sex)])[order(diet, sex, sample)]
  mm[, `:=`(gene = factor(gene, levels = gene_lv),
            col  = factor(sample, levels = mice_ord$sample))]
  xlab_of <- setNames(paste0(mice_ord$sample, " (", mice_ord$sex, ")"), mice_ord$sample)

  expr_fill <- if (MAIN_SCALE == "log2")
    scale_fill_viridis_c(name = "log2(mean + 1)") else
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-MAIN_CLIP, MAIN_CLIP),
                         oob = scales::squish,
                         name = paste0(if (MAIN_SCALE == "zSD") "z vs SD mice" else "z across mice",
                                       "\nlog2 mean expr\n(clipped +/-", MAIN_CLIP, ")"))
  base_th <- theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top",
          legend.title = element_text(size = 7), legend.key.height = unit(6, "pt"),
          plot.margin = margin(2, 1, 2, 1))
  panel_expr <- function(dt, title, show_y) {
    ggplot(mm[diet == dt], aes(col, gene, fill = value)) + geom_tile() +
      scale_x_discrete(labels = xlab_of) + expr_fill + labs(title = title) + base_th +
      theme(axis.text.y = if (show_y) element_text(size = 6) else element_blank())
  }
  lfc_lim <- max(abs(hits$avg_log2FC))
  p_lfc <- ggplot(hits[, .(gene = factor(gene, levels = gene_lv), avg_log2FC, col = "log2FC")],
                  aes(col, gene, fill = avg_log2FC)) +
    geom_tile() + geom_text(aes(label = sprintf("%.2f", avg_log2FC)), size = 1.8) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-lfc_lim, lfc_lim), name = "log2FC\nHFD / SD") +
    labs(title = "log2FC") + base_th + theme(axis.text.y = element_blank())
  pg <- lapply(list(panel_expr("SD", "SD", TRUE), p_lfc, panel_expr("HFD", "HFD", FALSE)),
               ggplotGrob)
  n_col <- c(sum(mice_ord$diet == "SD"), 1.6, sum(mice_ord$diet == "HFD"))   # panel widths, in tiles
  for (i in 1:3) {
    pc <- unique(pg[[i]]$layout[grepl("^panel", pg[[i]]$layout$name), "l"])
    pg[[i]]$widths[pc] <- unit(n_col[i], "null")
  }
  g <- cbind(pg[[1]], pg[[2]], pg[[3]], size = "max")
  ttl <- sprintf("%s: HFD vs SD, %d per-cell hits (p_bonferroni < %g, |log2FC| > %g)",
                 gsub("_", " ", ct_safe), nrow(hits), BONF_CUT, LFC_CUT)
  g <- gtable::gtable_add_rows(g, unit(14, "pt"), pos = 0)
  g <- gtable::gtable_add_grob(g, grid::textGrob(ttl, gp = grid::gpar(fontsize = 9, fontface = "bold")),
                               t = 1, l = 1, r = ncol(g))
  main_pages[[ct_safe]] <- list(g = g, n = nrow(hits), m = nrow(mice_ord))
  msg("main heatmap: ", ct_safe, " (", nrow(hits), " genes)")
}
# One PDF, one page per cell type; page size fits the largest cell type.
if (length(main_pages)) {
  pdf(file.path(DEG_DIR, "heatmap_main_all.pdf"),
      width  = 3 + 0.3 * (max(sapply(main_pages, `[[`, "m")) + 1.6) + 1.5,
      height = 2 + 0.13 * max(sapply(main_pages, `[[`, "n")), onefile = TRUE)
  for (pgi in main_pages) { grid::grid.newpage(); grid::grid.draw(pgi$g) }
  invisible(dev.off())
  msg("main heatmap: wrote heatmap_main_all.pdf (", length(main_pages), " pages)")
} else msg("main heatmap: no cell type with >= 2 hits")

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
#             enrichment_percell_KEGG.pdf, enrichment_percell_WP.pdf.
# gost queries the g:Profiler web server, so the node needs internet access. A failed
# query is logged and skipped; the plots above are already saved by then.

OUTDIR   <- file.path(DEG_DIR, "enrichment")
SOURCES  <- c(KEGG = "KEGG", WP = "WikiPathways")
MIN_HITS <- 5     # a cell type needs this many hits to be queried
TOP_N    <- 15    # terms shown per cell type
MAX_TERM <- 200   # drop very large terms (often generic: "Metabolic pathways")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
msg("enrichment: querying g:Profiler")

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
  ggsave(file.path(OUTDIR, paste0("enrichment_percell_", src, ".pdf")), gp,
         width = 15, height = 1.5 + ceiling(n_ct / 3) * (0.6 + 0.16 * n_rows),
         bg = "white", limitsize = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done")
