#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# projection_background.R  --  is a UMI threshold enough to remove ambient
# transgene RNA from the projection calls, or does it need background removal?
# Run by projection_background.sh. Asked for 2026-09-23.
#
# THE LOGIC: only excitatory neurons can carry a retrograde label. So every
# label on a NON-excitatory nucleus is background (ambient cassette RNA, see
# step 2.8), and those nuclei are an internal negative control. Their labelling
# rate, per animal x reporter x sequencing depth, is the empirical ambient
# false-positive rate. Applied to the excitatory nuclei it predicts how many
# excitatory labels are ambient too -> a false-call fraction per threshold.
#
#   control   = every annotated cell type except Excitatory Neurons and
#               Unidentified (unknown identity, so in neither group)
#   depth     = nCount_mouse (transgene UMIs excluded, step 5), cut into global
#               deciles. Ambient pick-up scales with how well a droplet was
#               captured, and excitatory nuclei are the deepest, so an
#               unstratified rate would understate their background.
#   expected  = sum over deciles of  control rate (animal, ORF, decile)
#                                   x excitatory nuclei (animal, decile)
#               a decile with no control nuclei in an animal falls back to that
#               animal x ORF's overall control rate.
#   upper     = 97.5% exact Poisson bound on the control positives, carried
#               through to `expected` by ratio. At 0 control positives the
#               bound uses the unstratified rate, so a high threshold that
#               simply ran out of data does not read as "background = 0".
#
# THE CALL is projection_figures.R's: counts summed over the three tracks
# (base + BLAST-rescued + antisense), a nucleus is labelled when EXACTLY ONE
# ORF reaches the threshold. Same counts, same rule, only k varies (1..10).
#
# Caveat on the control: glia can carry neuronal material and doublets exist,
# so the control rate may slightly OVER-state ambient -- conservative.
#
# Outputs, in Results/projection_background:
#   threshold_summary.csv          one row per threshold, pooled
#   threshold_by_animal_orf.csv    threshold x animal x ORF
#   fig_background_by_threshold.png
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(data.table)
  library(Matrix)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
OBJ_IN  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
REVISED <- file.path(PROJECT, "Results", "reporting", "reclass_annot_revised.csv.gz")
RESCUE  <- file.path(PROJECT, "Results", "2_8_rescue_assign")
ANTI    <- file.path(PROJECT, "Results", "2_10_antisense_harvest")
OUTDIR  <- file.path(PROJECT, "Results", "projection_background")

ORFS    <- c("EGFP", "tdTomato", "EBFP2")
REGIONS <- c("BLA", "NAc", "LHr")
EXC     <- "Excitatory Neurons"
NOT_CTRL <- c(EXC, "Unidentified")
KS      <- 1:10
N_BINS  <- 10L
MOUSE_LEVELS <- paste0("Ech", 1:12)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
pois_up <- function(x) qchisq(0.975, 2 * (x + 1)) / 2

## ---- object, annot revision, TG counts (as projection_figures.R) ----------

obj <- readRDS(OBJ_IN)
msg("loaded ", format(ncol(obj), big.mark = ","), " nuclei")
for (col in c("annot", "mouse", "sample", "gem_group", "nCount_mouse", REGIONS))
  if (!col %in% names(obj@meta.data)) stop("object has no '", col, "' column")

rv <- fread(REVISED)[changed == TRUE]
md <- data.table(barcode = colnames(obj), annot = as.character(obj$annot))
md[rv, on = "barcode", annot := i.annot_azimuth]
if (sum(md$barcode %in% rv$barcode) != nrow(rv)) stop("revised barcodes absent from object")

if (length(Layers(obj[["TG"]], search = "counts")) > 1) obj <- JoinLayers(obj, assay = "TG")
tgm <- GetAssayData(obj, assay = "TG", layer = "counts")
tot <- sapply(ORFS, function(f) as.numeric(tgm[f, ]))
rownames(tot) <- colnames(obj)

gem <- unique(as.data.table(obj@meta.data)[, .(sample = as.character(sample),
                                               gem_group = as.integer(gem_group))])
for (track in c("blast", "antisense")) {
  for (i in seq_len(nrow(gem))) {
    s <- gem$sample[i]
    if (track == "blast") {
      d <- fread(file.path(RESCUE, paste0(s, "_rescued_reads.tsv")), colClasses = "character")
      d <- d[in_called_cell == "1" & nzchar(CB) & nzchar(UB), .(feature = construct, CB, UB)]
    } else {
      d <- fread(file.path(ANTI, paste0(s, "_antisense_reads.tsv")), colClasses = "character")
      d <- d[nzchar(CB) & nzchar(UB), .(feature = contig, CB, UB)]
    }
    d[, feature := sub("_cassette$", "", feature)]
    d <- unique(d[feature %in% ORFS])            # a molecule, not a read
    d[, barcode := sub("-1$", paste0("-", gem$gem_group[i]), CB)]
    d <- d[barcode %in% rownames(tot), .N, by = .(barcode, feature)]
    tot[cbind(match(d$barcode, rownames(tot)), match(d$feature, ORFS))] <-
      tot[cbind(match(d$barcode, rownames(tot)), match(d$feature, ORFS))] + d$N
  }
}
msg("ORF UMIs after all three tracks: ",
    paste(ORFS, colSums(tot), sep = " ", collapse = ", "))

## ---- nucleus table --------------------------------------------------------

nuc <- data.table(barcode = colnames(obj),
                  mouse = factor(as.character(obj$mouse), levels = MOUSE_LEVELS),
                  annot = md$annot,
                  depth = as.numeric(obj$nCount_mouse))
nuc[, group := fcase(annot == EXC, "exc", annot %in% NOT_CTRL, "none", default = "ctrl")]
nuc[, bin := cut(depth, unique(quantile(depth, 0:N_BINS / N_BINS)),
                 include.lowest = TRUE, labels = FALSE)]
msg("excitatory ", format(nuc[group == "exc", .N], big.mark = ","),
    " | control ", format(nuc[group == "ctrl", .N], big.mark = ","),
    " | excluded (Unidentified) ", nuc[group == "none", .N])
msg("median depth: excitatory ", nuc[group == "exc", median(depth)],
    ", control ", nuc[group == "ctrl", median(depth)])
nuc <- nuc[group != "none"]
tot <- tot[nuc$barcode, , drop = FALSE]

## ---- per threshold --------------------------------------------------------

by_k <- rbindlist(lapply(KS, function(k) {
  det <- tot >= k
  single <- rowSums(det) == 1L
  orf <- ifelse(single, ORFS[max.col(det, ties.method = "first")], NA_character_)
  d <- copy(nuc)[, orf := orf]
  # denominators per animal x decile, positives per animal x ORF x decile
  den <- d[, .N, by = .(mouse, group, bin)]
  pos <- d[!is.na(orf), .N, by = .(mouse, group, bin, orf)]
  grid <- CJ(mouse = factor(MOUSE_LEVELS, levels = MOUSE_LEVELS), bin = sort(unique(d$bin)), orf = ORFS)
  grid <- merge(grid, dcast(den, mouse + bin ~ group, value.var = "N", fill = 0L),
                by = c("mouse", "bin"), all.x = TRUE)
  setnames(grid, c("ctrl", "exc"), c("n_ctrl", "n_exc"))
  p <- dcast(pos, mouse + bin + orf ~ group, value.var = "N", fill = 0L)
  for (g in c("ctrl", "exc")) if (!g %in% names(p)) p[, (g) := 0L]
  setnames(p, c("ctrl", "exc"), c("pos_ctrl", "pos_exc"))
  grid <- merge(grid, p, by = c("mouse", "bin", "orf"), all.x = TRUE)
  for (j in c("n_ctrl", "n_exc", "pos_ctrl", "pos_exc")) grid[is.na(get(j)), (j) := 0L]
  # animal x ORF overall control rate: fallback for an empty decile
  grid[, rate_all := sum(pos_ctrl) / sum(n_ctrl), by = .(mouse, orf)]
  grid[, rate := fifelse(n_ctrl > 0, pos_ctrl / n_ctrl, rate_all)]
  a <- grid[, .(n_exc = sum(n_exc), n_ctrl = sum(n_ctrl),
                obs_exc = sum(pos_exc), obs_ctrl = sum(pos_ctrl),
                exp_adj = sum(rate * n_exc)), by = .(mouse, orf)]
  a[, exp_unadj := obs_ctrl / n_ctrl * n_exc]
  a[, threshold := k]
  a
}))

finish <- function(x) {
  x[, exp_upper := fifelse(obs_ctrl > 0, exp_adj * pois_up(obs_ctrl) / obs_ctrl,
                           pois_up(0) / n_ctrl * n_exc)]
  x[, `:=`(
    pct_exc_of_labelled = round(100 * obs_exc / (obs_exc + obs_ctrl), 1),
    false_pct_adj   = round(100 * pmin(exp_adj,   obs_exc) / obs_exc, 1),
    false_pct_upper = round(100 * pmin(exp_upper, obs_exc) / obs_exc, 1),
    false_pct_unadj = round(100 * pmin(exp_unadj, obs_exc) / obs_exc, 1),
    exp_adj = round(exp_adj, 1), exp_upper = round(exp_upper, 1),
    exp_unadj = round(exp_unadj, 1))]
  x[, true_exc_est := round(pmax(obs_exc - exp_adj, 0), 1)]
  x[]
}

per <- finish(copy(by_k))
setcolorder(per, c("threshold", "mouse", "orf"))
setorder(per, threshold, mouse, orf)

# pooled: sum counts first, then the ratio -- never an average of ratios
pooled <- by_k[, .(n_exc = sum(n_exc) / length(ORFS), n_ctrl = sum(n_ctrl) / length(ORFS),
                   obs_exc = sum(obs_exc), obs_ctrl = sum(obs_ctrl),
                   exp_adj = sum(exp_adj), exp_unadj = sum(exp_unadj)), by = threshold]
pooled <- finish(pooled)
anim <- per[, .(true_exc_est = sum(true_exc_est)), by = .(threshold, mouse)]
pooled <- merge(pooled, anim[, .(animals_ge1_est = sum(true_exc_est >= 1),
                                 animals_ge10_est = sum(true_exc_est >= 10)),
                             by = threshold], by = "threshold")
setcolorder(pooled, c("threshold", "obs_exc", "obs_ctrl", "pct_exc_of_labelled",
                      "exp_adj", "exp_upper", "false_pct_adj", "false_pct_upper",
                      "exp_unadj", "false_pct_unadj", "true_exc_est",
                      "animals_ge1_est", "animals_ge10_est", "n_exc", "n_ctrl"))

fwrite(pooled, file.path(OUTDIR, "threshold_summary.csv"))
fwrite(per,    file.path(OUTDIR, "threshold_by_animal_orf.csv"))
msg("pooled, per threshold:")
print(pooled[, .(threshold, obs_exc, obs_ctrl, pct_exc_of_labelled, exp_adj,
                 false_pct_adj, false_pct_upper, false_pct_unadj, true_exc_est,
                 animals_ge1_est, animals_ge10_est)])
msg("per reporter, pooled over animals:")
print(per[, .(obs_exc = sum(obs_exc), obs_ctrl = sum(obs_ctrl),
              exp_adj = round(sum(exp_adj), 1),
              false_pct_adj = round(100 * min(sum(exp_adj), sum(obs_exc)) / sum(obs_exc), 1)),
          by = .(threshold, orf)][threshold %in% c(1, 2, 3, 5)])

## ---- figure: two panels, one axis each (no dual axis) ---------------------

long <- rbind(
  pooled[, .(threshold, panel = "Estimated false-call fraction among excitatory labels (%)",
             y = false_pct_adj, lo = NA_real_, hi = false_pct_upper)],
  pooled[, .(threshold, panel = "Labelled excitatory nuclei (observed)",
             y = as.numeric(obs_exc), lo = NA_real_, hi = NA_real_)])
p <- ggplot(long, aes(threshold, y)) +
  geom_linerange(aes(ymin = y, ymax = hi), colour = "grey60", linewidth = 0.6, na.rm = TRUE) +
  geom_line(colour = "#0f7ba2", linewidth = 0.7) +
  geom_point(colour = "#0f7ba2", size = 2) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_x_continuous(breaks = KS, name = "Threshold (>= UMI on one ORF, three tracks summed)") +
  labs(y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(hjust = 0))
ggsave(file.path(OUTDIR, "fig_background_by_threshold.png"), p, width = 7, height = 7, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
