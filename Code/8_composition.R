#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 8_composition.R  --  cell-type composition, run by 8_composition.sh.
#
# Reads the manually annotated object and asks whether the cell-type mix moves
# with diet. Five panels, one figure, no titles on it (user, 2026-09-01):
#
#   A  proportion of each cell type in each of the 12 animals
#   B  the same pooled by diet (SD vs HFD)
#   C  the same pooled by diet x sex
#   D  propeller test, HFD vs SD, one point per cell type
#   E  the per-animal proportions those tests are computed on
#
# The cell-type call is the `annot` column written by the manual pass at the end
# of 7_annotation_featureplots.R, NOT cluster_label/cell_class -- those come
# from the cluster_annotation.csv route the user abandoned.
#
# THE TEST IS AT THE ANIMAL LEVEL (user, 2026-09-01). scProportionTest was
# dropped for this reason: it permutes CELL labels, so with ~250k nuclei from 12
# mice it asks whether a reshuffle of nuclei could produce the split, which is a
# question about sampling depth, not about diet -- it treats one mouse's 20k
# nuclei as 20k independent observations and will call effects that are really
# one animal being different.
#
# speckle::propeller instead. n = 12 animals, one proportion per cell type per
# animal, asin-sqrt transformed to stabilise the variance, then a limma
# moderated t-test with empirical Bayes shrinkage across cell types -- the
# shrinkage is what makes 6 vs 6 usable. The design is ~ 0 + diet + sex, so sex
# is blocked out rather than left in the residual: the 2x2 is balanced 3 animals
# per cell. PropRatio in the output is HFD / SD.
#
# asin-sqrt rather than logit because a rare cell type can be absent from an
# animal, and logit needs a pseudo-count for those zeros.
#
# A per-cell-type Wilcoxon on the same 12 proportions (6 vs 6) is written to the
# same csv as a non-parametric cross-check. It borrows no information across
# cell types, so it is much less powerful -- its smallest attainable two-sided p
# is 1/462 = 0.0022, and after FDR nothing may survive. Where the two disagree,
# the disagreement is the finding: it means the effect rests on the shrinkage.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(speckle)
  library(limma)
})

source("/scratch/user/uqschauq/Insular/Code/lib_celltype_palette.R")

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
IN_UMAP <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_umap.rds")
OUTDIR  <- file.path(PROJECT, "Results", "8_composition")

ANNOT     <- "annot"
TRANSFORM <- "asin"
FDR_CUT   <- 0.05
COL_SD    <- "#0d366b"
COL_HFD   <- "#d1341f"
SEED      <- 42

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

if (!file.exists(IN))
  stop("no annotated object at ", IN,
       " -- run the saveRDS line at the end of 7_annotation_featureplots.R first")
if (file.exists(IN_UMAP) && file.mtime(IN) < file.mtime(IN_UMAP))
  msg("WARNING: ", basename(IN), " is OLDER than ", basename(IN_UMAP),
      " -- check it carries the current clustering and the manual annot")

msg("reading ", basename(IN))
obj <- readRDS(IN)
if (!ANNOT %in% names(obj@meta.data))
  stop("no '", ANNOT, "' column in the object -- the manual annotation was not saved")
stopifnot(all(c("sample", "mouse", "diet", "sex") %in% names(obj@meta.data)))
msg("loaded ", format(ncol(obj), big.mark = ","), " nuclei")

md <- as.data.table(obj@meta.data, keep.rownames = "barcode")
md[, annot := as.character(get(ANNOT))]
rm(obj); invisible(gc())

# Cell-type order: the object's own factor levels, reversed so the order reads
# the same way as the dotplot (excitatory first). Levels the manual pass never
# renamed are still bare cluster numbers; they stay in as their own category
# rather than being silently dropped.
lv <- levels(md[[ANNOT]])
ct_levels <- if (!is.null(lv)) rev(lv) else md[, .N, by = annot][order(-N), annot]
unassigned <- grep("^[0-9]+$", ct_levels, value = TRUE)
if (length(unassigned))
  msg("NOTE: ", length(unassigned), " level(s) still a bare cluster number, kept as-is: ",
      paste(unassigned, collapse = ", "))
md[, annot := factor(annot, levels = ct_levels)]

md[, diet := factor(as.character(diet), levels = c("SD", "HFD"))]
md[, sex  := factor(as.character(sex),  levels = c("F", "M"))]
md[, diet_sex := factor(paste(diet, sex), levels = c("SD F", "SD M", "HFD F", "HFD M"))]
setorder(md, diet, sex, sample)
md[, sample := factor(as.character(sample), levels = unique(as.character(sample)))]
msg(nlevels(md$annot), " cell types | ", nlevels(md$sample), " animals")
print(unique(md[, .(sample, mouse, diet, sex)]))

## ---- proportion tables -----------------------------------------------------

prop_by <- function(dt, by) {
  x <- dt[, .(n = .N), by = c(by, "annot")]
  full <- CJ(a = levels(droplevels(dt[[by]])), b = levels(dt$annot))
  setnames(full, c(by, "annot"))
  full[, (by)  := factor(get(by), levels = levels(dt[[by]]))]
  full[, annot := factor(annot,   levels = levels(dt$annot))]
  x <- x[full, on = c(by, "annot")]
  x[is.na(n), n := 0L]
  x[, prop := n / sum(n), by = by]
  x[]
}

p_sample <- prop_by(md, "sample")
p_sample <- merge(p_sample, unique(md[, .(sample, mouse, diet, sex, diet_sex)]),
                  by = "sample", sort = FALSE)
p_diet   <- prop_by(md, "diet")
p_ds     <- prop_by(md, "diet_sex")

fwrite(p_sample[order(annot, sample)], file.path(OUTDIR, "composition_per_sample.csv"))
fwrite(p_diet[order(annot, diet)],     file.path(OUTDIR, "composition_by_diet.csv"))
fwrite(p_ds[order(annot, diet_sex)],   file.path(OUTDIR, "composition_by_diet_sex.csv"))

# Animal-level summary: mean +/- sd of the per-animal proportions. These are the
# numbers to quote in text -- the pooled ones above weight animals by yield.
per_animal <- p_sample[, .(n_animals = .N, mean_prop = mean(prop), sd_prop = sd(prop),
                           min_prop = min(prop), max_prop = max(prop)),
                       by = .(annot, diet)]
per_animal <- dcast(per_animal, annot ~ diet,
                    value.var = c("mean_prop", "sd_prop", "min_prop", "max_prop"))
msg("per-animal mean proportions by diet (%):")
print(per_animal[, .(annot, SD = round(100 * mean_prop_SD, 2), sd_SD = round(100 * sd_prop_SD, 2),
                     HFD = round(100 * mean_prop_HFD, 2), sd_HFD = round(100 * sd_prop_HFD, 2))])

## ---- propeller, animal level -----------------------------------------------

props <- getTransformedProps(clusters = md$annot, sample = md$sample,
                             transform = TRANSFORM)
animals <- unique(md[, .(sample, mouse, diet, sex)])
animals <- animals[match(colnames(props$TransformedProps), as.character(sample))]
stopifnot(identical(as.character(animals$sample), colnames(props$TransformedProps)))
msg("propeller design: ", nrow(animals), " animals, ",
    paste(sprintf("%s=%d", levels(animals$diet), table(animals$diet)), collapse = " "),
    " | sex blocked (", paste(sprintf("%s=%d", levels(animals$sex), table(animals$sex)),
                              collapse = " "), ")")

design <- model.matrix(~ 0 + diet + sex, data = animals)
ctr <- setNames(rep(0, ncol(design)), colnames(design))
stopifnot(all(c("dietSD", "dietHFD") %in% names(ctr)))
ctr["dietHFD"] <- 1
ctr["dietSD"]  <- -1

res <- propeller.ttest(props, design, contrasts = ctr,
                       robust = TRUE, trend = FALSE, sort = FALSE)
res <- as.data.table(res, keep.rownames = "annot")
setnames(res, sub("^PropMean\\.diet", "PropMean_", names(res)))
res[, annot := factor(annot, levels = ct_levels)]
res[, log2FC := log2(PropRatio)]

# Non-parametric cross-check on the same 12 numbers, borrowing nothing across
# cell types.
wil <- p_sample[, .(wilcox_p = tryCatch(
  suppressWarnings(wilcox.test(prop ~ diet)$p.value), error = function(e) NA_real_)),
  by = annot]
wil[, wilcox_FDR := p.adjust(wilcox_p, "BH")]

res <- merge(res, wil, by = "annot", all.x = TRUE)
res <- merge(res, per_animal, by = "annot", all.x = TRUE)
setorder(res, P.Value)
fwrite(res, file.path(OUTDIR, "propeller_HFD_vs_SD.csv"))
msg(res[FDR < FDR_CUT, .N], " cell type(s) at propeller FDR < ", FDR_CUT,
    "; ", res[wilcox_FDR < FDR_CUT, .N], " at Wilcoxon FDR < ", FDR_CUT)
print(res[, .(annot, SD = round(100 * mean_prop_SD, 2), HFD = round(100 * mean_prop_HFD, 2),
              log2FC = round(log2FC, 3), P.Value = signif(P.Value, 3),
              FDR = signif(FDR, 3), wilcox_p = signif(wilcox_p, 3))])

## ---- figure ----------------------------------------------------------------

pal <- celltype_cols(levels(md$annot))

stacked <- function(dt, xcol, rotate = FALSE) {
  p <- ggplot(dt, aes(.data[[xcol]], prop, fill = annot)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "proportion of nuclei") +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
  if (rotate) p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
}

pA <- stacked(p_sample, "sample", rotate = TRUE)
pB <- stacked(p_diet, "diet")
pC <- stacked(p_ds, "diet_sex")

yscale <- scale_y_discrete(limits = rev(ct_levels))

res[, sig := factor(ifelse(FDR < FDR_CUT, paste("FDR <", FDR_CUT), "n.s."),
                    levels = c(paste("FDR <", FDR_CUT), "n.s."))]
pD <- ggplot(res, aes(log2FC, annot, colour = sig)) +
  geom_vline(xintercept = 0) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = setNames(c(COL_HFD, "grey60"), levels(res$sig)),
                      name = NULL, drop = FALSE) +
  yscale +
  labs(x = "log2 (HFD / SD) of mean animal proportion", y = NULL) +
  theme_bw(base_size = 11)

# The 12 numbers behind panel D. Zeros cannot go on a log axis, so an animal
# with none of a cell type is drawn at half the smallest non-zero proportion.
floor_p <- min(p_sample[prop > 0, prop]) / 2
nz <- p_sample[prop == 0, .N]
if (nz) msg(nz, " animal x cell-type combination(s) with zero nuclei, drawn at ",
            signif(floor_p, 3))
p_sample[, prop_plot := pmax(prop, floor_p)]
dodge <- position_dodge(width = 0.7)
pE <- ggplot(p_sample, aes(prop_plot, annot, colour = diet)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.55, linewidth = 0.25,
               position = dodge, show.legend = FALSE) +
  geom_point(size = 1.5, alpha = 0.85, position = dodge) +
  scale_colour_manual(values = c(SD = COL_SD, HFD = COL_HFD), name = NULL) +
  scale_x_log10(labels = scales::percent_format(accuracy = 0.1)) +
  yscale +
  labs(x = "proportion of nuclei per animal", y = NULL) +
  theme_bw(base_size = 11)

fig <- pA / (pB | pC) / (pD | pE) +
  plot_layout(heights = c(1, 1, 1.3), guides = "collect")

f <- file.path(OUTDIR, "composition_overview")
ggsave(paste0(f, ".png"), fig, width = 14, height = 15, dpi = 300)
ggsave(paste0(f, ".pdf"), fig, width = 14, height = 15)
msg("written ", basename(f), ".png/.pdf")

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
