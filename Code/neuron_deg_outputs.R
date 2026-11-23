# Outputs of the per-cell DEG for the neuronal runs, from the CSVs written by
# percell_deg_diet_sexfollowup.R (no Seurat object needed).
#
# Usage:   Rscript neuron_deg_outputs.R [excitatory_layer | inhibitory | excitatory]
#   excitatory_layer  excitatory nuclei pooled by layer (L2/3, L4/5, L5, L6; LAYER_MAP in the
#                     DEG script)                       -> .../percell_wilcox_diet/excitatory_layers
#   inhibitory        Azimuth inhibitory subclasses     -> .../inhibitory_subclasses
#
# Hit (DEG): pooled-sex Wilcoxon p_bonferroni < BONF_CUT (Bonferroni over all genes in the
#            object, Seurat's p_val_adj) and |avg_log2FC| > LFC_CUT.
# 1  volcano_<run>.<ext>       every tested gene, one facet per group.
# 2  heatmap_allDEG.pdf        every hit, one page per group: SD mean | per-cell log2FC | HFD mean.
# 2b heatmap_log2FC_sex.pdf    every hit, one page per group: log2FC (HFD / SD) pooled, in
#                              females and in males; red = higher in HFD.
# 3  sex_venn_F_vs_M.<ext>     female vs male within-sex DEG per group (2 circles), and
#    sex_venn_F_vs_M_updown.<ext>  the same split by direction (4 ellipses).
#    Sex DEG: within-sex p_bonferroni < BONF_CUT (Bonferroni over the group's pooled hits) and
#    |within-sex log2FC| > LFC_CUT. Only pooled hits were tested within sex.
# 4  enrichment/enrichment_<KEGG|WP>.<ext>  gprofiler2::gost per group, hits vs the genes tested
#    in that group, g:SCS-adjusted p < 0.05; table in enrichment/enrichment_all.csv.
# No plot titles or subtitles (thresholds are stated here); multi-page heatmaps carry only
# the group name so the pages can be told apart.
# <ext> = FIG_EXT. PDF by default: the JPEG device hung in batch jobs on Bunya. Set
# FIG_EXT <- "jpeg" when running interactively if you want JPEGs.
# The enrichment step needs internet (g:Profiler web server); run it on a node that has it.

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

args     <- commandArgs(trailingOnly = TRUE)
RUN      <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "excitatory_layer"
SUBDIR   <- c(excitatory_layer = "excitatory_layers", inhibitory = "inhibitory_subclasses",
              excitatory = "excitatory_subclasses")[RUN]
if (is.na(SUBDIR)) stop("argument must be one of: excitatory_layer, inhibitory, excitatory")
# Same root as the DEG script: env PERCELL_OUT if set, else QRISdata.
DEG_DIR  <- file.path(Sys.getenv("PERCELL_OUT", "/QRISdata/Q9732/scratch/Results/percell_wilcox_diet"), SUBDIR)
FIG_EXT  <- "pdf"
BONF_CUT <- 0.05
LFC_CUT  <- 0.5
HEAT_MAX <- 50       # used only if MAIN_FILE gets a "top" entry
VOLC_LABEL <- 8      # hits labelled per volcano facet (lowest p_val)
DROP_TYPES <- "Unidentified"
DO_ENRICH  <- TRUE

safe   <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
pretty <- function(x) gsub("_", " ", sub("^L(\\d)_(\\d)", "L\\1/\\2", x))   # L2_3 -> L2/3
msg    <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
is_hit <- function(d) d$p_bonferroni < BONF_CUT & abs(d$avg_log2FC) > LFC_CUT

deg_files <- list.files(DEG_DIR, pattern = "^DEG_percell_HFDvsSD_.*\\.csv$", full.names = TRUE)
deg_files <- deg_files[!basename(deg_files) %in% paste0("DEG_percell_HFDvsSD_", safe(DROP_TYPES), ".csv")]
if (!length(deg_files)) stop("no DEG_percell_HFDvsSD_*.csv in ", DEG_DIR,
                             " -- run percell_deg_diet_sexfollowup.R ", RUN, " first")
group_of <- function(f, prefix) sub(paste0("^", prefix, "(.*)\\.csv$"), "\\1", basename(f))
msg("run = ", RUN, "; ", length(deg_files), " groups in ", DEG_DIR)

## ---- 1. volcano ---------------------------------------------------------------------
# p = 0 (Wilcoxon underflow) is drawn at the group's smallest non-zero p, as a triangle.
# Bonferroni line: Seurat multiplies by the number of genes in the object, recovered here
# from p_bonferroni / p_val.
v <- rbindlist(lapply(deg_files, function(f) {
  d <- fread(f)
  pmin_nz <- min(d[p_val > 0, p_val], 1e-300)
  d[, .(annot = pretty(group_of(f, "DEG_percell_HFDvsSD_")), gene, avg_log2FC, p_val, p_bonferroni,
        capped = p_val == 0, y = -log10(pmax(p_val, pmin_nz)),
        cls = fifelse(is_hit(d) & avg_log2FC > 0, "up in HFD",
              fifelse(is_hit(d) & avg_log2FC < 0, "down in HFD", "n.s.")))]
}))
n_bonf <- v[p_val > 0 & p_bonferroni < 1, median(p_bonferroni / p_val)]
v[, cls := factor(cls, levels = c("up in HFD", "down in HFD", "n.s."))]
setorder(v, -cls, -p_val)                                # n.s. drawn first, hits on top
n_hit <- v[, .(lab = sprintf("%s\n%d up, %d down", annot[1], sum(cls == "up in HFD"),
                             sum(cls == "down in HFD"))), by = annot]
v[, facet := factor(n_hit$lab[match(annot, n_hit$annot)], levels = n_hit[order(annot), lab])]
lab <- v[cls != "n.s."][order(p_val), head(.SD, VOLC_LABEL), by = annot]
gp <- ggplot(v, aes(avg_log2FC, y, colour = cls)) +
  geom_point(aes(shape = capped), size = 0.5, alpha = 0.7) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ facet, scales = "free", ncol = 4) +
  scale_colour_manual(values = c("up in HFD" = "#b2182b", "down in HFD" = "#2166ac", "n.s." = "grey75"),
                      drop = FALSE, name = paste0("p_bonferroni < ", BONF_CUT, ", |log2FC| > ", LFC_CUT)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none")
if (is.finite(n_bonf))
  gp <- gp + geom_hline(yintercept = -log10(BONF_CUT / n_bonf), linetype = "dashed", linewidth = 0.3)
gp <- gp +
  (if (requireNamespace("ggrepel", quietly = TRUE))
     ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.3, max.overlaps = Inf,
                              min.segment.length = 0, segment.size = 0.2, show.legend = FALSE)
   else geom_text(data = lab, aes(label = gene), size = 2.3, vjust = -0.6, show.legend = FALSE)) +
  labs(x = "avg_log2FC (HFD / SD)", y = "-log10 per-cell p (Wilcoxon)") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(DEG_DIR, paste0("volcano_", RUN, ".", FIG_EXT)), gp, width = 14, dpi = 300,
       height = 1.2 + 3.2 * ceiling(uniqueN(v$annot) / 4), bg = "white", limitsize = FALSE)
msg("1. volcano written")
rm(v, lab)

## ---- 2. heatmap of every DEG: SD mean | log2FC | HFD mean ----------------------------
# Rows: every pooled hit, sorted by avg_log2FC (most up at the top).
# Left / right: mean expression of the diet's 6 mice. Per mouse, mouse_means_hits_<group>.csv
#   holds the mean of expm1(RNA data) over that mouse's nuclei (the average Seurat's
#   avg_log2FC is built from); the column is the mean of those 6 values, shown as
#   log2(mean + 1) on one colour scale shared by SD and HFD, value printed in the tile.
# Middle: per-cell DEG avg_log2FC (HFD / SD); red = higher in HFD.
# With more than MAIN_LABEL_MAX genes, gene names and values are left out.
MAIN_LABEL_MAX <- 150
heat_pages <- list()
mm_files <- list.files(DEG_DIR, pattern = "^mouse_means_hits_.*\\.csv$", full.names = TRUE)
mm_files <- mm_files[!basename(mm_files) %in% paste0("mouse_means_hits_", safe(DROP_TYPES), ".csv")]
for (fm in mm_files) {
  ct <- group_of(fm, "mouse_means_hits_")
  fd <- file.path(DEG_DIR, paste0("DEG_percell_HFDvsSD_", ct, ".csv"))
  if (!file.exists(fd)) next
  deg  <- fread(fd)
  mm   <- fread(fm)
  hits <- deg[is_hit(deg) & gene %in% mm$gene]
  if (nrow(hits) < 2) next
  lv     <- hits[order(avg_log2FC), gene]
  lab_ok <- nrow(hits) <= MAIN_LABEL_MAX
  gm <- mm[gene %in% hits$gene, .(n_mice = .N, value = log2(mean(mean_expr) + 1)), by = .(gene, diet)]
  gm[, `:=`(gene = factor(gene, levels = lv), diet = factor(diet, levels = c("SD", "HFD")))]
  n_mice <- gm[, max(n_mice), by = diet][, setNames(V1, diet)]
  lf <- hits[, .(gene = factor(gene, levels = lv), lfc = avg_log2FC)]
  th <- theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 8), legend.position = "top",
          legend.title = element_text(size = 6.5), legend.key.height = unit(5, "pt"),
          legend.key.width = unit(14, "pt"), plot.margin = margin(2, 1, 2, 1))
  expr_scale <- scale_fill_viridis_c(limits = range(gm$value), name = "log2(mean expr + 1)")
  side <- function(dt, show_y) {
    ggplot(gm[diet == dt], aes(diet, gene, fill = value)) + geom_tile() +
      (if (lab_ok) geom_text(aes(label = sprintf("%.2f", value),
                                 colour = value > mean(range(gm$value))), size = 1.8, show.legend = FALSE)) +
      scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
      expr_scale + labs(title = dt) + th +
      theme(axis.text.y = if (show_y && lab_ok) element_text(size = 6) else element_blank())
  }
  lfc_lim <- max(abs(lf$lfc))
  p_mid <- ggplot(lf, aes("log2FC", gene, fill = lfc)) + geom_tile() +
    (if (lab_ok) geom_text(aes(label = sprintf("%.2f", lfc)), size = 1.8)) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-lfc_lim, lfc_lim), name = "log2FC HFD / SD") +
    labs(title = "log2FC") + th + theme(axis.text.y = element_blank())
  pg <- lapply(list(side("SD", TRUE), p_mid, side("HFD", FALSE)), ggplotGrob)
  g  <- cbind(pg[[1]], pg[[2]], pg[[3]], size = "max")
  ttl <- pretty(ct)   # group name only, to tell the pages apart
  g <- gtable::gtable_add_rows(g, unit(14, "pt"), pos = 0)
  g <- gtable::gtable_add_grob(g, grid::textGrob(ttl, gp = grid::gpar(fontsize = 8.5, fontface = "bold")),
                               t = 1, l = 1, r = ncol(g))
  heat_pages[[ct]] <- list(g = g, n = nrow(hits))
}
if (length(heat_pages)) {
  pdf(file.path(DEG_DIR, "heatmap_allDEG.pdf"), width = 6.5,
      height = 2.2 + 0.13 * min(max(sapply(heat_pages, `[[`, "n")), MAIN_LABEL_MAX), onefile = TRUE)
  for (pgi in heat_pages) { grid::grid.newpage(); grid::grid.draw(pgi$g) }
  invisible(dev.off())
}
msg("2. heatmap written (", length(heat_pages), " pages)")

## ---- 2b. log2FC heatmap: pooled | female | male --------------------------------------
# Rows: every pooled hit, sorted by pooled avg_log2FC (most up at the top). Columns: per-cell
# Wilcoxon avg_log2FC (HFD / SD) with sexes pooled, within females, within males (from
# DEG_percell_sexsplit_<group>.csv). Red = higher in HFD, blue = lower. "*" = DEG in that
# column: the pooled column holds hits only, so every row has one there; F / M use the
# within-sex cut (p_bonferroni < BONF_CUT over the group's hits, |log2FC| > LFC_CUT).
# Grey = not tested in that sex. One colour scale per page, symmetric, clipped at LFC_CLIP.
LFC_CLIP <- 2
sx_hm <- list.files(DEG_DIR, pattern = "^DEG_percell_sexsplit_.*\\.csv$", full.names = TRUE)
sx_hm <- sx_hm[!basename(sx_hm) %in% paste0("DEG_percell_sexsplit_", safe(DROP_TYPES), ".csv")]
lfc_pages <- list()
for (f in sx_hm) {
  ct <- group_of(f, "DEG_percell_sexsplit_")
  d  <- fread(f)
  if (!all(c("F_avg_log2FC", "M_avg_log2FC", "F_p_bonferroni", "M_p_bonferroni") %in% names(d))) next
  d  <- d[is_hit(d)]
  if (nrow(d) < 2) next
  lv <- d[order(avg_log2FC), gene]
  long <- rbind(
    d[, .(gene, col = "Pooled", lfc = avg_log2FC, sig = TRUE)],
    d[, .(gene, col = "Female", lfc = F_avg_log2FC,
          sig = !is.na(F_p_bonferroni) & F_p_bonferroni < BONF_CUT & abs(F_avg_log2FC) > LFC_CUT)],
    d[, .(gene, col = "Male",   lfc = M_avg_log2FC,
          sig = !is.na(M_p_bonferroni) & M_p_bonferroni < BONF_CUT & abs(M_avg_log2FC) > LFC_CUT)])
  long[, `:=`(gene = factor(gene, levels = lv), col = factor(col, levels = c("Pooled", "Female", "Male")))]
  lab_ok <- nrow(d) <= MAIN_LABEL_MAX
  n_sig  <- long[sig == TRUE, .N, by = col]
  gp <- ggplot(long, aes(col, gene, fill = lfc)) +
    geom_tile(colour = if (lab_ok) "white" else NA, linewidth = 0.2) +
    (if (lab_ok) geom_text(data = long[sig == TRUE], aes(label = "*"), size = 2.6, vjust = 0.75)) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", na.value = "grey80",
                         limits = c(-LFC_CLIP, LFC_CLIP), oob = scales::squish,
                         name = paste0("log2FC\nHFD / SD\n(clipped +/-", LFC_CLIP, ")")) +
    scale_x_discrete(position = "top",
                     labels = function(x) sprintf("%s\n(%d DEG)", x, n_sig$N[match(x, n_sig$col)])) +
    labs(x = NULL, y = NULL, title = pretty(ct)) +   # group name only, to tell the pages apart
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text.x.top = element_text(face = "bold"),
          axis.text.y = if (lab_ok) element_text(size = 6) else element_blank())
  lfc_pages[[ct]] <- list(gp = gp, n = nrow(d))
}
if (length(lfc_pages)) {
  pdf(file.path(DEG_DIR, "heatmap_log2FC_sex.pdf"), width = 5.5,
      height = 2.2 + 0.13 * min(max(sapply(lfc_pages, `[[`, "n")), MAIN_LABEL_MAX), onefile = TRUE)
  for (pg in lfc_pages) print(pg$gp)
  invisible(dev.off())
  msg("2b. log2FC heatmap written (", length(lfc_pages), " pages)")
} else msg("2b. no sexsplit tables with F/M columns; log2FC heatmap skipped")

## ---- 3. sex Venn diagrams --------------------------------------------------------------
sx_files <- list.files(DEG_DIR, pattern = "^DEG_percell_sexsplit_.*\\.csv$", full.names = TRUE)
sx_files <- sx_files[!basename(sx_files) %in% paste0("DEG_percell_sexsplit_", safe(DROP_TYPES), ".csv")]
if (length(sx_files)) {
vn <- rbindlist(lapply(sx_files, function(f) {
  d <- fread(f)
  if (!all(c("F_p_bonferroni", "M_p_bonferroni") %in% names(d))) return(NULL)
  d <- d[p_bonferroni < BONF_CUT & abs(avg_log2FC) > LFC_CUT]
  d[, .(annot = group_of(f, "DEG_percell_sexsplit_"),
        gene, F_avg_log2FC, M_avg_log2FC, F_p_bonferroni, M_p_bonferroni)]
}))
vn[, `:=`(F_deg = !is.na(F_p_bonferroni) & F_p_bonferroni < BONF_CUT & abs(F_avg_log2FC) > LFC_CUT,
          M_deg = !is.na(M_p_bonferroni) & M_p_bonferroni < BONF_CUT & abs(M_avg_log2FC) > LFC_CUT)]
vn[, region := fifelse(F_deg & M_deg, "both", fifelse(F_deg, "F only", fifelse(M_deg, "M only", "neither")))]
fwrite(vn[order(annot, region, gene)], file.path(DEG_DIR, "sex_venn_membership.csv"))

### 3a. two circles
cnt <- vn[, data.table(x = c(-1.05, 0, 1.05),
                       n = c(sum(region == "F only"), sum(region == "both"), sum(region == "M only")),
                       n_hit = .N, n_neither = sum(region == "neither")), by = annot]
cnt[, facet := sprintf("%s\n%d pooled hits; %d DEG in neither sex", pretty(annot), n_hit, n_neither)]

th <- seq(0, 2 * pi, length.out = 200)
circ <- rbind(data.table(sex = "F", x = -0.6 + cos(th), y = sin(th)),
              data.table(sex = "M", x =  0.6 + cos(th), y = sin(th)))
gp <- ggplot() +
  geom_polygon(data = circ, aes(x, y, fill = sex, group = sex), alpha = 0.25, colour = "grey30", linewidth = 0.3) +
  geom_text(data = data.table(x = c(-1.2, 1.2), y = 1.2, lab = c("Female", "Male")),
            aes(x, y, label = lab), size = 3, fontface = "bold") +
  geom_text(data = cnt, aes(x, 0, label = n), size = 4.5) +
  facet_wrap(~ facet, ncol = 4) +
  scale_fill_manual(values = c(F = "#d95f02", M = "#1b9e77"), guide = "none") +
  coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.1, 1.4)) +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(size = 7, margin = margin(2, 0, 2, 0)))
ggsave(file.path(DEG_DIR, paste0("sex_venn_F_vs_M.", FIG_EXT)), gp, width = 12, dpi = 300, bg = "white",
       height = 1.2 + 2.6 * ceiling(uniqueN(vn$annot) / 4), limitsize = FALSE)
print(dcast(vn, annot ~ region, fun.aggregate = length))

### 3b. four ellipses, by direction
# Four sets per cell type: F up, M up, F down, M down (up = log2FC > 0, higher in HFD).
# F up / F down (and M up / M down) cannot share a gene, so regions containing both are 0.
SETS <- c("F up", "M up", "F down", "M down")
vn[, `:=`(`F up`   = F_deg & F_avg_log2FC > 0, `F down` = F_deg & F_avg_log2FC < 0,
          `M up`   = M_deg & M_avg_log2FC > 0, `M down` = M_deg & M_avg_log2FC < 0)]

# Classic 4-ellipse Venn layout; every one of the 15 regions exists.
ell <- data.table(set = SETS, x0 = c(0.38, 0.47, 0.53, 0.62), y0 = c(0.45, 0.55, 0.55, 0.45),
                  a = 0.35, b = 0.2, rot = c(-45, -45, 45, 45) * pi / 180)
inside <- function(x, y, e) {   # scaled radius: < 1 inside the ellipse
  dx <- x - e$x0; dy <- y - e$y0
  u <- dx * cos(e$rot) + dy * sin(e$rot); v <- -dx * sin(e$rot) + dy * cos(e$rot)
  sqrt((u / e$a)^2 + (v / e$b)^2)
}
# Label position of each region: the grid point of that region farthest from any outline.
grid <- CJ(x = seq(0, 1, by = 0.004), y = seq(0, 1.05, by = 0.004))
r <- sapply(seq_len(nrow(ell)), function(i) inside(grid$x, grid$y, ell[i]))
grid[, key := apply(r < 1, 1, function(z) paste(as.integer(z), collapse = ""))]
grid[, depth := apply(abs(r - 1), 1, min)]
pos <- grid[key != "0000"][order(-depth), .SD[1], by = key][, .(key, x, y)]

vn[, key := paste0(as.integer(`F up`), as.integer(`M up`), as.integer(`F down`), as.integer(`M down`))]
cnt <- CJ(annot = unique(vn$annot), key = pos$key)
cnt <- merge(cnt, vn[key != "0000", .N, by = .(annot, key)], by = c("annot", "key"), all.x = TRUE)
cnt[is.na(N), N := 0L]
cnt <- merge(cnt, pos, by = "key")
strip <- vn[, .(lab = sprintf("%s\n%d pooled hits; %d DEG in neither sex", pretty(annot[1]),
                              .N, sum(key == "0000"))), by = annot]
cnt[, facet := strip$lab[match(annot, strip$annot)]]

th <- seq(0, 2 * pi, length.out = 200)
outl <- ell[, {
  u <- a * cos(th); v <- b * sin(th)
  .(x = x0 + u * cos(rot) - v * sin(rot), y = y0 + u * sin(rot) + v * cos(rot))
}, by = set]
outl[, set := factor(set, levels = SETS)]
# set names just beyond each ellipse's upper end
lab_sets <- ell[, .(set, x = x0 - 1.1 * a * cos(rot) * sign(rot) * -1, y = y0 + 1.1 * a * abs(sin(rot)) + 0.04)]
gp <- ggplot() +
  geom_polygon(data = outl, aes(x, y, fill = set, group = set), alpha = 0.2, colour = "grey30", linewidth = 0.3) +
  geom_text(data = lab_sets, aes(x, y, label = set, colour = set), size = 2.8, fontface = "bold") +
  geom_text(data = cnt, aes(x, y, label = N), size = 2.6) +
  facet_wrap(~ facet, ncol = 4) +
  scale_fill_manual(values = c("F up" = "#d95f02", "M up" = "#1b9e77", "F down" = "#7570b3", "M down" = "#e7298a"),
                    guide = "none", aesthetics = c("fill", "colour")) +
  coord_equal(xlim = c(-0.02, 1.02), ylim = c(0.02, 1.02)) +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(size = 7, margin = margin(2, 0, 2, 0)))
ggsave(file.path(DEG_DIR, paste0("sex_venn_F_vs_M_updown.", FIG_EXT)), gp, width = 12, dpi = 300, bg = "white",
       height = 1.2 + 3 * ceiling(uniqueN(vn$annot) / 4), limitsize = FALSE)
print(dcast(cnt, annot ~ key, value.var = "N"))   # key digits: F up, M up, F down, M down
msg("3. Venn diagrams written")
} else msg("3. no DEG_percell_sexsplit_*.csv; Venn skipped")

## ---- 4. pathway enrichment ------------------------------------------------------------
# Query: per group, all hits (up and down together). Background: the genes FindMarkers tested
# in that group (domain_scope = "custom"); the default, all annotated mouse genes, would make
# any pathway of brain-expressed genes look enriched. p: gost's g:SCS-adjusted p_value.
# Bar plots: one per source, x = -log10 adjusted p, TOP_N terms per group, facet by group.
if (DO_ENRICH) {
  suppressPackageStartupMessages(library(gprofiler2))
  OUTDIR   <- file.path(DEG_DIR, "enrichment")
  SOURCES  <- c(KEGG = "KEGG", WP = "WikiPathways")
  MIN_HITS <- 5     # a group needs this many hits to be queried
  TOP_N    <- 15    # terms shown per group
  MAX_TERM <- 200   # drop very large terms (often generic: "Metabolic pathways")
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  res <- list()
  for (f in deg_files) {
    ct   <- group_of(f, "DEG_percell_HFDvsSD_")
    d    <- fread(f)
    hits <- d[is_hit(d), gene]
    msg("4. ", ct, ": ", nrow(d), " tested, ", length(hits), " hits")
    if (length(hits) < MIN_HITS) next
    g <- tryCatch(gost(query = hits, organism = "mmusculus", sources = names(SOURCES),
                       custom_bg = d$gene, domain_scope = "custom", user_threshold = 0.05,
                       correction_method = "g_SCS", significant = TRUE, evcodes = TRUE),
                  error = function(err) { msg("   gost failed: ", conditionMessage(err)); NULL })
    if (is.null(g)) { msg("   no term at p < 0.05"); next }
    res[[ct]] <- as.data.table(g$result)[, .(annot = pretty(ct), source, term_id, term_name,
                                             p_value, term_size, query_size, intersection_size,
                                             precision, recall, intersection)]
  }
  if (length(res)) {
    res <- rbindlist(res)
    setorder(res, source, annot, p_value)
    fwrite(res, file.path(OUTDIR, "enrichment_all.csv"))
    for (src in names(SOURCES)) {
      d <- res[source == src & term_size <= MAX_TERM]
      if (!nrow(d)) { msg("4. ", src, ": nothing to plot"); next }
      d <- d[order(p_value), head(.SD, TOP_N), by = annot]
      d[, mlp := -log10(p_value)]
      d[, key := paste(term_name, annot, sep = "___")]   # y key unique per facet, ordered by p
      d[, key := factor(key, levels = d[order(annot, mlp), key])]
      n_rows <- d[, .N, by = annot][, max(N)]
      n_ct   <- uniqueN(d$annot)
      gp <- ggplot(d, aes(mlp, key)) +
        geom_col(fill = "grey40", width = 0.8) +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
        facet_wrap(~ annot, scales = "free_y", ncol = min(3, n_ct)) +
        scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
        labs(x = "-log10 adjusted p (g:SCS)", y = NULL) +
        theme_bw(base_size = 8) +
        theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
      ggsave(file.path(OUTDIR, paste0("enrichment_", src, ".", FIG_EXT)), gp, dpi = 300,
             width = 5 * min(3, n_ct), height = 1.5 + ceiling(n_ct / 3) * (0.6 + 0.16 * n_rows),
             bg = "white", limitsize = FALSE)
    }
    msg("4. enrichment written")
  } else msg("4. no enriched term in any group")
}
msg("done")
