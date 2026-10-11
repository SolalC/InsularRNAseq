#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10_5_subclass_composition.R  --  does diet shift the composition of the
# neuronal subsets? Run by 10_5_subclass_composition.sh.
#
# Step 8 asked this of the 13 `annot` cell types across the whole object. This
# asks it one level down, of the Azimuth subclasses WITHIN each neuronal
# lineage: among a given animal's excitatory nuclei, what fraction is L5 IT,
# L2/3 IT, L6 CT ...; and among its inhibitory nuclei, what fraction is Pvalb,
# Sst, Vip, Lamp5 ... The denominator is the lineage, not the whole object, so
# a change in how many neurons an animal yielded cannot masquerade as a change
# in the mix within them.
#
# THE TEST is step 8's, unchanged, so the two are directly comparable:
# speckle::propeller on animal-level proportions, n = 12, asin-sqrt transformed
# (not logit -- a rare subclass can be absent from an animal and logit needs a
# pseudo-count for those zeros), design ~ 0 + diet + sex so sex is blocked
# rather than ignored, contrast HFD - SD, limma robust empirical-Bayes
# shrinkage across subclasses. That shrinkage is what makes 6 vs 6 usable.
#
# A per-subclass Wilcoxon on the same 12 proportions is written alongside as a
# non-parametric cross-check. It borrows nothing across subclasses, so it is far
# less powerful -- its smallest attainable two-sided p is 1/462 = 0.0022, and
# after FDR nothing may survive. Where the two disagree, THAT is the finding:
# it means the effect rests on the shrinkage.
#
# WHICH NUCLEI COUNT. Only those whose Azimuth subclass is modally the
# lineage's own class, i.e. the same set 10_3 tests and the complement of what
# 10_4 reclassifies away. A nucleus sitting in the excitatory subset that maps
# to Astro or to Pvalb is not an excitatory subtype and would distort the mix.
# The excluded counts are reported per animal rather than dropped silently.
#
# No Seurat object is read: the 10_2 prediction tables already carry mouse,
# sample, diet and sex per barcode. Outputs go to Results/reporting.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(speckle)
  library(limma)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
OUTDIR  <- file.path(PROJECT, "Results", "reporting")

RES       <- 0.2
TRANSFORM <- "asin"
FDR_CUT   <- 0.05
LINEAGE_CLASS <- c(excitatory = "Glutamatergic", inhibitory = "GABAergic")
LABEL_ORDER <- list(
  excitatory = c("L2/3 IT", "L4/5 IT", "L5 IT", "L5 ET", "L5/6 NP",
                 "L6 IT", "L6 IT Car3", "L6 CT", "L6b"),
  inhibitory = c("Pvalb", "Sst", "Sst Chodl", "Vip", "Lamp5", "Sncg", "Meis2"))
DIET_COLS <- c(SD = "#0f7ba2", HFD = "#dd5129")
MOUSE_LEVELS <- paste0("Ech", 1:12)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

run_lineage <- function(lin) {
  msg("=== ", lin, " ===")
  f <- file.path(PROJECT, "Results", sprintf("10_2_%s_azimuth", lin),
                 sprintf("azimuth_predictions_res%s.rds", RES))
  if (!file.exists(f))
    stop("no Azimuth predictions for ", lin, " -- run 10_2_azimuth.sh ", lin)
  d <- as.data.table(readRDS(f))
  d <- d[, .(barcode, mouse, sample = as.character(sample), diet, sex,
             subclass = as.character(predicted.subclass),
             pred_class = as.character(predicted.class))]

  modal <- d[, .N, by = .(subclass, pred_class)][order(subclass, -N)][, .SD[1], by = subclass]
  keep  <- modal[pred_class == LINEAGE_CLASS[lin], subclass]
  excl  <- d[!subclass %in% keep, .N, by = .(mouse, subclass)]
  msg("  excluding ", format(sum(excl$N), big.mark = ","), " of ",
      format(nrow(d), big.mark = ","), " nuclei not modally ",
      LINEAGE_CLASS[lin], " (see 10_4 for where they go)")
  fwrite(excl[order(-N)],
         file.path(OUTDIR, sprintf("subclass_composition_%s_excluded.csv", lin)))
  d <- d[subclass %in% keep]

  lv <- c(intersect(LABEL_ORDER[[lin]], unique(d$subclass)),
          setdiff(unique(d$subclass), LABEL_ORDER[[lin]]))
  d[, subclass := factor(subclass, levels = lv)]
  d[, mouse := factor(mouse, levels = MOUSE_LEVELS)]
  d[, diet := factor(diet, levels = c("SD", "HFD"))]
  d[, sex  := factor(sex,  levels = c("F", "M"))]

  # per-animal proportions WITHIN the lineage
  ps <- d[, .N, by = .(sample, mouse, diet, sex, subclass)]
  ps[, prop := N / sum(N), by = sample]
  ps <- ps[CJ(sample = unique(sample), subclass = factor(lv, levels = lv), unique = TRUE),
           on = .(sample, subclass)]
  ps[is.na(N), N := 0L][is.na(prop), prop := 0]
  ps[, c("mouse", "diet", "sex") := NULL]
  ps <- merge(ps, unique(d[, .(sample, mouse, diet, sex)]), by = "sample")
  fwrite(ps[order(subclass, diet, mouse)],
         file.path(OUTDIR, sprintf("subclass_composition_%s_per_animal.csv", lin)))

  per_animal <- ps[, .(mean_prop = mean(prop), sd_prop = sd(prop)), by = .(subclass, diet)]
  per_animal <- dcast(per_animal, subclass ~ diet, value.var = c("mean_prop", "sd_prop"))

  props   <- getTransformedProps(clusters = d$subclass, sample = d$sample,
                                 transform = TRANSFORM)
  animals <- unique(d[, .(sample, mouse, diet, sex)])
  animals <- animals[match(colnames(props$TransformedProps), sample)]
  stopifnot(identical(animals$sample, colnames(props$TransformedProps)))
  msg("  design: ", nrow(animals), " animals, ",
      paste(sprintf("%s=%d", levels(animals$diet), table(animals$diet)), collapse = " "),
      " | sex blocked")

  design <- model.matrix(~ 0 + diet + sex, data = animals)
  ctr <- setNames(rep(0, ncol(design)), colnames(design))
  stopifnot(all(c("dietSD", "dietHFD") %in% names(ctr)))
  ctr["dietHFD"] <- 1; ctr["dietSD"] <- -1

  res <- as.data.table(propeller.ttest(props, design, contrasts = ctr,
                                       robust = TRUE, trend = FALSE, sort = FALSE),
                       keep.rownames = "subclass")
  setnames(res, sub("^PropMean\\.diet", "PropMean_", names(res)))
  res[, log2FC := log2(PropRatio)]

  wil <- ps[, .(wilcox_p = tryCatch(
    suppressWarnings(wilcox.test(prop ~ diet)$p.value), error = function(e) NA_real_)),
    by = subclass]
  wil[, wilcox_FDR := p.adjust(wilcox_p, "BH")]

  res <- merge(res, wil, by = "subclass", all.x = TRUE)
  res <- merge(res, per_animal, by = "subclass", all.x = TRUE)
  res[, subclass := factor(subclass, levels = lv)]
  setorder(res, P.Value)
  fwrite(res, file.path(OUTDIR, sprintf("subclass_composition_%s_propeller.csv", lin)))

  msg("  ", res[FDR < FDR_CUT, .N], " subclass(es) at propeller FDR < ", FDR_CUT,
      "; ", res[wilcox_FDR < FDR_CUT, .N], " at Wilcoxon FDR < ", FDR_CUT)
  print(res[, .(subclass, SD = round(100 * mean_prop_SD, 2),
                HFD = round(100 * mean_prop_HFD, 2), log2FC = round(log2FC, 3),
                P.Value = signif(P.Value, 3), FDR = signif(FDR, 3),
                wilcox_p = signif(wilcox_p, 3))])

  scols <- setNames(scales::hue_pal()(nlevels(ps$subclass)), levels(ps$subclass))
  p_bar <- ggplot(ps, aes(mouse, prop, fill = subclass)) +
    geom_col(width = 0.78) +
    scale_fill_manual(values = scols, name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = paste0("proportion of that animal's ", lin, " nuclei")) +
    facet_grid(~ diet, scales = "free_x", space = "free_x") +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA))

  # The propeller p and FDR for that subclass, drawn in its own facet. x = 1.5
  # is the midpoint of the two discrete diet positions; y = Inf anchors to the
  # top of each free-scaled panel.
  lab <- res[, .(subclass, prop = 0,
                 lab = sprintf("p = %s\nFDR = %s", signif(P.Value, 2), signif(FDR, 2)),
                 sig = FDR < FDR_CUT)]
  p_dot <- ggplot(ps, aes(diet, prop, colour = diet)) +
    geom_point(position = position_jitter(width = 0.14, height = 0), size = 1.9,
               show.legend = FALSE) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.3,
                 colour = "#3d3d3a", show.legend = FALSE) +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = lab, colour = sig),
              inherit.aes = FALSE, vjust = 1.25, size = 2.9, lineheight = 0.95,
              show.legend = FALSE) +
    facet_wrap(~ subclass, scales = "free_y", nrow = 2) +
    scale_colour_manual(values = c(DIET_COLS, `TRUE` = "#bb292c", `FALSE` = "#3d3d3a")) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       expand = expansion(mult = c(0.05, 0.28))) +
    labs(x = NULL, y = "proportion of the lineage") +
    guides(colour = "none") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

  fig <- p_bar / p_dot + plot_layout(heights = c(1, 1.15))
  fig <- fig & labs(title = NULL, subtitle = NULL) &
    theme(plot.title = element_blank(), plot.subtitle = element_blank())
  nm <- sprintf("fig_%s_subclass_composition_%s",
                if (lin == "excitatory") "25" else "26", lin)
  ggsave(file.path(OUTDIR, paste0(nm, ".png")), fig, width = 12, height = 9, dpi = 300)
  ggsave(file.path(OUTDIR, paste0(nm, ".pdf")), fig, width = 12, height = 9)
  msg("  ", nm)
  res[, lineage := lin][]
}

all_res <- rbindlist(lapply(names(LINEAGE_CLASS), run_lineage), fill = TRUE)
fwrite(all_res, file.path(OUTDIR, "subclass_composition_propeller_all.csv"))
writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, "sessionInfo_subclass_composition.txt"))
msg("done: ", OUTDIR)
