#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 5_seurat_qc.R
#
# Per-nucleus QC on the ambient-corrected, doublet-flagged object from step 4.
#
# WHY THIS RUNS AFTER STEP 4 (order changed 2026-08-31, user):
#   The cutoff rule is per-animal mean +/- SD_N sd on nCount and nFeature. A
#   doublet is by construction high on both, so leaving doublets in the
#   distribution pulls the mean up and inflates the sd, pushing the upper bound
#   out past the very nuclei it exists to catch. Removing them first makes the
#   bound describe the singlet population it is meant to describe. The size of
#   that bias is measured, not asserted: qc_bounds_doublet_bias.csv gives both
#   sets of bounds side by side.
#
#   CellBender likewise belongs upstream: it needs the EMPTY droplets to learn
#   the ambient profile, so it can only run on the per-sample raw matrices.
#
# Cutoff rule, three parts:
#
#  - An ABSOLUTE floor on nFeature, NFEATURE_MIN, applied to every animal.
#    Added 2026-08-31 (user) after the sd rule failed. A dispersion rule answers
#    "is this nucleus unusual for this animal"; it cannot answer "is this a
#    nucleus at all", and the two questions need different instruments. Ech8 and
#    Ech9 carried a population of droplets with ONE detected gene (nFeature p01
#    = 1, nCount p01 = 9 and 11, log10_genes_per_umi p01 = 0). On the log10
#    scale those sit at exactly 0, which inflated the sd enough that
#    mean - 3 sd back-transformed to a floor of 1.0 gene for Ech8 and 2.5 for
#    Ech9 -- the contamination set its own floor, and ~2,000 of those nuclei
#    survived into step 6 as cluster 18 (97.6% Ech8, zero cells from ten of the
#    twelve animals).
#
#  - nCount and nFeature bounded per ANIMAL by DISPERSION on that animal's own
#    SINGLET distribution, on the log10 scale, so a shallower library is not
#    penalised against a deeper one. DISPERSION = "mad" (user, 2026-08-31):
#    median +/- SD_N MAD. The mean/sd version is still computed and written to
#    qc_cutoffs_per_sample.csv for comparison, but it is not what filters.
#    MAD is robust to the low tail that broke the sd rule: it gives Ech8 a floor
#    of 139 rather than 1.0, and Ech9 73.9 rather than 2.5.
#
#  - percent_mt is a fixed ceiling, MT_MAX -- NOT a dispersion rule, and
#    deliberately not switched to MAD with the others. After CellBender the
#    median percent_mt is 0 in all twelve animals, so the MAD is also 0 and
#    median + k*MAD collapses to a ceiling of 0, which would delete every
#    nucleus carrying a single mitochondrial read. percent_mt_hi_mad in
#    qc_cutoffs_sd_vs_mad.csv is 0 for all twelve; that column is a warning, not
#    a candidate. See the warning at MT_MAX: the value was chosen on
#    UNCORRECTED data and ambient removal moves this distribution.
#
# Transgene note: the object's feature names have already been through
# CreateSeuratObject, which rewrites "_" to "-", so AAV_shared is "AAV-shared"
# here. Every metric is computed twice -- the plain Seurat value over all 33,700
# features and a "_mouse" value over the 33,696 mouse genes. Threshold on the
# _mouse columns: AAV-shared can be large in a transgene-positive nucleus, and
# leaving it in the denominator inflates nCount and deflates percent_mt for
# exactly the nuclei the projection arm depends on.
# ---------------------------------------------------------------------------

library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)

## ---- CONFIG ---------------------------------------------------------------

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_ambient_doublets.rds")
OUTDIR  <- file.path(PROJECT, "Results", "5_seurat_qc")
OBJDIR  <- file.path(PROJECT, "Data", "Seurat")
CLEANRDS <- file.path(OBJDIR, "HFDseq_v2_filtered.rds")

MT_PATTERN <- "^mt-"          # GENCODE vM33 is lowercase; "^MT-" matches nothing
N_MT_EXPECTED <- 13
TRANSGENES <- c("EGFP", "tdTomato", "EBFP2", "AAV_shared")

SD_N <- 3
LOG_COUNTS <- TRUE            # log10 scale (user, 2026-08-29): raw-scale floors
                              # were negative in all 12 animals and removed nothing

# Dispersion rule for the nFeature/nCount bounds: "mad" (median +/- SD_N MAD) or
# "sd" (mean +/- SD_N sd). "mad" since 2026-08-31 (user). Both are computed and
# written either way; this only selects which one filters.
DISPERSION <- "mad"

# Absolute nFeature floor, applied per animal as pmax(dispersion_floor,
# NFEATURE_MIN) so it can only ever tighten a bound, never loosen one. See the
# header: this is the part that does not depend on the shape of the
# distribution it is filtering.
NFEATURE_MIN <- 200

# percent_mt ceiling, 1% (user, 2026-08-31). RE-DECIDED on corrected data --
# this supersedes the 5% chosen 2026-08-29, which had been read off a histogram
# of UNCORRECTED counts where it did essentially all the filtering (3,544 of
# 3,568 removals).
#
# In snRNA-seq a nucleus contains no mitochondria, so percent_mt IS ambient,
# and CellBender has removed it: the median is now 0 in all twelve animals and
# the distribution is squeezed against zero. Measured against the 249,671
# nuclei that passed the MAD + NFEATURE_MIN filter at MT_MAX = 5, tightening
# the ceiling costs:
#
#     mt > 3     325 nuclei  (0.13%)
#     mt > 2     799          (0.32%)
#     mt > 1   2,413          (0.97%)   <- chosen
#     mt > 0.5 5,965          (2.39%)
#
# So 5% was very nearly inert post-CellBender and 1% is where the ceiling
# starts doing real work. Flat across the design at 2%: 0.27% HFD vs 0.37% SD,
# 0.27% M vs 0.37% F, so this is not confounded with diet or sex.
MT_MAX <- 1
MT_UPPER_ONLY <- TRUE

# Remove the step 4 doublet calls before anything is measured. This is the
# reason the step order changed; setting it FALSE reproduces the old behaviour.
DROP_DOUBLETS <- TRUE

APPLY_FILTER <- TRUE
SCATTER_N <- 60000
set.seed(1)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")

## ---- READ -----------------------------------------------------------------

if (!file.exists(IN)) stop("missing step 4 object: ", IN)
obj <- readRDS(IN)
DefaultAssay(obj) <- "RNA"
msg("loaded ", ncol(obj), " nuclei x ", nrow(obj), " features from ", basename(IN))
stopifnot(all(c("sample", "mouse", "diet", "sex", "dbl_class", "dbl_score") %in%
              colnames(obj@meta.data)))

TRANSGENES_OBJ <- gsub("_", "-", TRANSGENES, fixed = TRUE)
missing_tg <- setdiff(TRANSGENES_OBJ, rownames(obj))
if (length(missing_tg))
  stop("transgene features absent: ", paste(missing_tg, collapse = ", "))

mt_genes <- grep(MT_PATTERN, rownames(obj), value = TRUE)
msg("mitochondrial genes matched: ", length(mt_genes), " (", paste(mt_genes, collapse = ", "), ")")
if (length(mt_genes) != N_MT_EXPECTED)
  stop("expected ", N_MT_EXPECTED, " genes matching ", MT_PATTERN, ", found ", length(mt_genes))

## ---- METRICS (computed on every nucleus, before any removal) ---------------

counts <- LayerData(obj, assay = "RNA", layer = "counts")
tg_umi      <- Matrix::colSums(counts[TRANSGENES_OBJ, , drop = FALSE])
tg_detected <- Matrix::colSums(counts[TRANSGENES_OBJ, , drop = FALSE] > 0)
mt_umi      <- Matrix::colSums(counts[mt_genes, , drop = FALSE])
rm(counts); invisible(gc())

obj$nCount_transgene   <- as.numeric(tg_umi)
obj$nFeature_transgene <- as.integer(tg_detected)
obj$nCount_mouse       <- as.numeric(obj$nCount_RNA)   - as.numeric(tg_umi)
obj$nFeature_mouse     <- as.integer(obj$nFeature_RNA) - as.integer(tg_detected)
obj$percent_mt         <- 100 * as.numeric(mt_umi) / pmax(obj$nCount_mouse, 1)
obj$percent_mt_all     <- 100 * as.numeric(mt_umi) / pmax(as.numeric(obj$nCount_RNA), 1)
obj$umis_per_gene      <- obj$nCount_mouse / pmax(obj$nFeature_mouse, 1)
obj$log10_genes_per_umi <- log10(pmax(obj$nFeature_mouse, 1)) / log10(pmax(obj$nCount_mouse, 10))

qc_all <- as.data.table(obj@meta.data, keep.rownames = "barcode")

## ---- BOUNDS FUNCTIONS -----------------------------------------------------

sd_bounds <- function(v, k, log_scale, upper_only = FALSE) {
  x <- if (log_scale) log10(pmax(v, 1)) else as.numeric(v)
  m <- mean(x); s <- stats::sd(x)
  lo <- m - k * s; hi <- m + k * s
  if (log_scale) { lo <- 10^lo; hi <- 10^hi }
  c(lower = if (upper_only) NA_real_ else lo, upper = hi)
}

mad_bounds <- function(v, k, log_scale, upper_only = FALSE) {
  x <- if (log_scale) log10(pmax(v, 1)) else as.numeric(v)
  m <- stats::median(x); s <- stats::mad(x)
  lo <- m - k * s; hi <- m + k * s
  if (log_scale) { lo <- 10^lo; hi <- 10^hi }
  c(lower = if (upper_only) NA_real_ else lo, upper = hi)
}

disp_bounds <- function(v, k, log_scale, upper_only = FALSE) {
  f <- if (DISPERSION == "mad") mad_bounds else sd_bounds
  f(v, k, log_scale, upper_only = upper_only)
}

bounds_for <- function(d, k = SD_N) {
  fl <- sd_bounds(d$nFeature_mouse, k, TRUE);  fr <- sd_bounds(d$nFeature_mouse, k, FALSE)
  cl <- sd_bounds(d$nCount_mouse,   k, TRUE);  cr <- sd_bounds(d$nCount_mouse,   k, FALSE)
  mt <- sd_bounds(d$percent_mt,     k, FALSE, upper_only = MT_UPPER_ONLY)
  fm <- mad_bounds(d$nFeature_mouse, k, LOG_COUNTS)
  cm <- mad_bounds(d$nCount_mouse,   k, LOG_COUNTS)
  mm <- mad_bounds(d$percent_mt,     k, FALSE, upper_only = MT_UPPER_ONLY)
  list(n_nuclei = nrow(d),
       nFeature_mean = round(mean(d$nFeature_mouse), 1),
       nFeature_sd   = round(stats::sd(d$nFeature_mouse), 1),
       nFeature_lo_log = round(fl[["lower"]], 1), nFeature_hi_log = round(fl[["upper"]], 1),
       nFeature_lo_raw = round(fr[["lower"]], 1), nFeature_hi_raw = round(fr[["upper"]], 1),
       nCount_mean = round(mean(d$nCount_mouse), 1),
       nCount_sd   = round(stats::sd(d$nCount_mouse), 1),
       nCount_lo_log = round(cl[["lower"]], 1), nCount_hi_log = round(cl[["upper"]], 1),
       nCount_lo_raw = round(cr[["lower"]], 1), nCount_hi_raw = round(cr[["upper"]], 1),
       percent_mt_mean = round(mean(d$percent_mt), 3),
       percent_mt_sd   = round(stats::sd(d$percent_mt), 3),
       percent_mt_hi_3sd = round(mt[["upper"]], 3),
       nFeature_lo_mad = round(fm[["lower"]], 1), nFeature_hi_mad = round(fm[["upper"]], 1),
       nCount_lo_mad   = round(cm[["lower"]], 1), nCount_hi_mad   = round(cm[["upper"]], 1),
       # 0 in all twelve animals post-CellBender -- reported as a warning, never used
       percent_mt_hi_mad = round(mm[["upper"]], 3),
       percent_mt_hi = MT_MAX)
}

## ---- DOUBLET REMOVAL, AND WHAT IT DOES TO THE BOUNDS -----------------------

dbl_tab <- qc_all[, .(n = .N, n_doublet = sum(dbl_class == "doublet")), by = .(sample, mouse)]
dbl_tab[, pct_doublet := round(100 * n_doublet / n, 2)]
msg("step 4 doublet calls: ", sum(dbl_tab$n_doublet), " of ", sum(dbl_tab$n),
    " (", round(100 * sum(dbl_tab$n_doublet) / sum(dbl_tab$n), 2), "%)")

# The point of the reorder, quantified: the same rule applied to the same
# animals, with and without the doublets in the distribution.
bias <- merge(
  qc_all[, bounds_for(.SD), by = sample][, .(sample,
      nFeature_hi_withdbl = nFeature_hi_log, nCount_hi_withdbl = nCount_hi_log,
      nFeature_mean_withdbl = nFeature_mean, nFeature_sd_withdbl = nFeature_sd)],
  qc_all[dbl_class == "singlet", bounds_for(.SD), by = sample][, .(sample,
      nFeature_hi_singlet = nFeature_hi_log, nCount_hi_singlet = nCount_hi_log,
      nFeature_mean_singlet = nFeature_mean, nFeature_sd_singlet = nFeature_sd)],
  by = "sample")
bias[, `:=`(nFeature_hi_shift_pct = round(100 * (nFeature_hi_withdbl - nFeature_hi_singlet) / nFeature_hi_singlet, 2),
            nCount_hi_shift_pct   = round(100 * (nCount_hi_withdbl - nCount_hi_singlet) / nCount_hi_singlet, 2))]
fwrite(bias, file.path(OUTDIR, "qc_bounds_doublet_bias.csv"))
msg("upper-bound inflation from leaving doublets in (this is why step 4 runs first):")
print(bias[, .(sample, nFeature_hi_withdbl, nFeature_hi_singlet, nFeature_hi_shift_pct,
               nCount_hi_shift_pct)])

if (DROP_DOUBLETS) {
  before <- ncol(obj)
  obj <- subset(obj, subset = dbl_class == "singlet")
  msg("doublets removed: ", before, " -> ", ncol(obj), " nuclei (",
      round(100 * ncol(obj) / before, 2), "% kept)")
} else {
  msg("DROP_DOUBLETS is FALSE - doublets retained, bounds will be inflated")
}
qc <- as.data.table(obj@meta.data, keep.rownames = "barcode")

## ---- SUMMARY TABLES -------------------------------------------------------

qmetrics <- c("nFeature_mouse", "nCount_mouse", "percent_mt",
              "umis_per_gene", "log10_genes_per_umi", "nCount_transgene")

summarise <- function(d) {
  out <- list(n_nuclei = nrow(d))
  for (m in qmetrics) {
    v <- d[[m]]; q <- quantile(v, c(0.01, 0.25, 0.5, 0.75, 0.99), names = FALSE)
    out[[paste0(m, "_mean")]]   <- round(mean(v), 3)
    out[[paste0(m, "_p01")]]    <- round(q[1], 3)
    out[[paste0(m, "_p25")]]    <- round(q[2], 3)
    out[[paste0(m, "_median")]] <- round(q[3], 3)
    out[[paste0(m, "_p75")]]    <- round(q[4], 3)
    out[[paste0(m, "_p99")]]    <- round(q[5], 3)
  }
  out
}
per_sample <- qc[, summarise(.SD), by = .(sample, mouse, diet, sex)]
overall <- cbind(data.table(sample = "ALL", mouse = "ALL", diet = NA, sex = NA),
                 as.data.table(summarise(qc)))
fwrite(rbind(per_sample, overall, fill = TRUE), file.path(OUTDIR, "qc_summary_per_sample.csv"))

cutoffs <- qc[, bounds_for(.SD), by = .(sample, mouse, diet, sex)]
stopifnot(DISPERSION %in% c("mad", "sd"))
sfx <- if (DISPERSION == "mad") "_mad" else if (LOG_COUNTS) "_log" else "_raw"
cutoffs[, `:=`(nFeature_lo_disp = get(paste0("nFeature_lo", sfx)),
               nFeature_hi = get(paste0("nFeature_hi", sfx)),
               nCount_lo   = get(paste0("nCount_lo", sfx)),
               nCount_hi   = get(paste0("nCount_hi", sfx)),
               scale_used  = if (LOG_COUNTS) "log10" else "raw",
               sd_n = SD_N, dispersion = DISPERSION, nfeature_min = NFEATURE_MIN,
               mt_rule = "fixed", doublets_dropped = DROP_DOUBLETS)]

# The absolute floor can only tighten a per-animal bound, never loosen one.
cutoffs[, nFeature_lo := pmax(nFeature_lo_disp, NFEATURE_MIN)]
cutoffs[, floor_from := fifelse(NFEATURE_MIN > nFeature_lo_disp, "absolute", DISPERSION)]
fwrite(cutoffs, file.path(OUTDIR, "qc_cutoffs_per_sample.csv"))
msg("per-sample ", SD_N, "-", DISPERSION, " cutoffs (",
    if (LOG_COUNTS) "log10" else "raw", " scale), nFeature floor >= ", NFEATURE_MIN, ":")
print(cutoffs[, .(sample, n_nuclei, nFeature_lo_disp, nFeature_lo, floor_from,
                  nFeature_hi, nCount_lo, nCount_hi, percent_mt_hi)])
msg("nFeature floor set by the absolute minimum in ",
    cutoffs[floor_from == "absolute", .N], " of ", nrow(cutoffs), " animals")

# Guard. The previous version tested nFeature_lo <= 0, which under
# LOG_COUNTS = TRUE is unsatisfiable: the bound is back-transformed as 10^x and
# is therefore strictly positive for every finite x. It was written on
# 2026-08-29 for the raw scale and became dead code when log10 was switched on,
# which is why the Ech8 floor of 1.0 gene passed without a word in the log.
# Test against the absolute floor instead, and stop rather than warn.
bad <- cutoffs[!is.finite(nFeature_lo) | !is.finite(nFeature_hi) |
               !is.finite(nCount_lo)   | !is.finite(nCount_hi) |
               nFeature_lo >= nFeature_hi | nCount_lo >= nCount_hi]
if (nrow(bad)) {
  print(bad[, .(sample, nFeature_lo, nFeature_hi, nCount_lo, nCount_hi)])
  stop("non-finite or inverted bounds in ", nrow(bad), " sample(s) - refusing to filter")
}
soft <- cutoffs[nFeature_lo_disp < NFEATURE_MIN / 4]
if (nrow(soft)) {
  msg("NOTE: the ", DISPERSION, " floor fell below ", round(NFEATURE_MIN / 4),
      " genes in ", nrow(soft), " sample(s) - the absolute floor is doing the work there:")
  print(soft[, .(sample, nFeature_lo_disp, nFeature_lo)])
}

## ---- FLAG -----------------------------------------------------------------

qc <- merge(qc, cutoffs[, .(sample, nFeature_lo, nFeature_hi, nCount_lo, nCount_hi, percent_mt_hi)],
            by = "sample", all.x = TRUE, sort = FALSE)
qc[, qc_pass := nFeature_mouse >= nFeature_lo & nFeature_mouse <= nFeature_hi &
                nCount_mouse   >= nCount_lo   & nCount_mouse   <= nCount_hi]
if (!is.na(MT_MAX)) {
  qc[, qc_pass := qc_pass & percent_mt <= MT_MAX]
  msg("percent_mt ceiling applied at ", MT_MAX, "%")
} else {
  msg("percent_mt is measured only - no ceiling")
}

setkey(qc, barcode)
qc <- qc[colnames(obj)]
stopifnot(identical(qc$barcode, colnames(obj)))
obj$qc_pass <- qc$qc_pass

# What each criterion costs on its own, so the filter is not a single number.
crit <- qc[, .(mito = sum(percent_mt > MT_MAX, na.rm = TRUE),
               below_nfeature_min = sum(nFeature_mouse < NFEATURE_MIN),
               below_floor = sum(nFeature_mouse < nFeature_lo | nCount_mouse < nCount_lo),
               above_ceiling = sum(nFeature_mouse > nFeature_hi | nCount_mouse > nCount_hi),
               total = .N)]
msg("removed by criterion (overlapping): mito ", crit$mito,
    ", below nFeature ", NFEATURE_MIN, " ", crit$below_nfeature_min,
    ", below floor ", crit$below_floor,
    ", above ceiling ", crit$above_ceiling, " of ", crit$total)
fwrite(crit, file.path(OUTDIR, "qc_removal_by_criterion.csv"))

retained <- qc[, .(total = .N, kept = sum(qc_pass)), by = sample]
retained[, pct_kept := round(100 * kept / total, 2)]
retained <- merge(retained, cutoffs[, .(sample, mouse, diet, sex)], by = "sample", sort = FALSE)
fwrite(retained, file.path(OUTDIR, "qc_retention_per_sample.csv"))
msg("retained at ", SD_N, " ", DISPERSION, ": ", sum(retained$kept), " / ", sum(retained$total), " nuclei (",
    round(100 * sum(retained$kept) / sum(retained$total), 2), "%)")
print(retained[, .(sample, total, kept, pct_kept)])

sens <- rbindlist(lapply(c(2, 2.5, 3, 3.5, 4), function(k) {
  b <- qc[, { f <- disp_bounds(nFeature_mouse, k, LOG_COUNTS)
              cc <- disp_bounds(nCount_mouse, k, LOG_COUNTS)
              flo <- max(f[["lower"]], NFEATURE_MIN)
              .(kept = sum(nFeature_mouse >= flo & nFeature_mouse <= f[["upper"]] &
                           nCount_mouse >= cc[["lower"]] & nCount_mouse <= cc[["upper"]]), total = .N) }, by = sample]
  b[, `:=`(sd_n = k, dispersion = DISPERSION, pct_kept = round(100 * kept / total, 2))]; b }))
fwrite(sens[order(sd_n, sample)], file.path(OUTDIR, "qc_sd_sensitivity.csv"))

mt_levels <- sort(unique(c(0.1, 0.25, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, if (!is.na(MT_MAX)) MT_MAX)))
exceed <- rbindlist(lapply(mt_levels, function(L)
  qc[, .(mt_level = L, n_above = sum(percent_mt > L), total = .N,
         pct_above = round(100 * mean(percent_mt > L), 3)), by = sample]))
fwrite(exceed[order(mt_level, sample)], file.path(OUTDIR, "percent_mt_exceedance.csv"))
msg("percent_mt after ambient removal: mean ", round(mean(qc$percent_mt), 3),
    ", median ", round(median(qc$percent_mt), 4),
    ", p99 ", round(quantile(qc$percent_mt, 0.99), 3),
    "  [MT_MAX = ", MT_MAX, "%, re-decided on corrected data 2026-08-31]")

# Both rules side by side. bounds_for() computes them together now, so this is
# a projection of `cutoffs` rather than a second pass over the nuclei.
fwrite(cutoffs[, .(sample, dispersion = DISPERSION, nfeature_min = NFEATURE_MIN,
                   nFeature_lo_applied = nFeature_lo, floor_from,
                   nFeature_lo_log, nFeature_hi_log, nCount_lo_log, nCount_hi_log,
                   nFeature_lo_mad, nFeature_hi_mad, nCount_lo_mad, nCount_hi_mad,
                   percent_mt_hi, percent_mt_hi_3sd, percent_mt_hi_mad)],
       file.path(OUTDIR, "qc_cutoffs_sd_vs_mad.csv"))

fwrite(qc[, c("barcode", "sample", "mouse", "diet", "sex", "gem_group",
              "nFeature_RNA", "nCount_RNA", qmetrics, "percent_mt_all",
              "dbl_score", "dbl_class",
              "nFeature_lo", "nFeature_hi", "nCount_lo", "nCount_hi",
              "percent_mt_hi", "qc_pass"), with = FALSE],
       file.path(OUTDIR, "qc_metrics_per_cell.csv.gz"), compress = "gzip")

## ---- PLOTS ----------------------------------------------------------------

theme_set(theme_bw(base_size = 11))
rot <- theme(axis.text.x = element_text(angle = 45, hjust = 1))

viol <- function(metric, ylab, logy, lo_col = NULL, hi_col = NULL) {
  p <- ggplot(qc, aes(x = sample, y = .data[[metric]], fill = diet)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2, fill = "white")
  if (!is.null(hi_col)) {
    bd <- cutoffs[, c("sample", lo_col, hi_col), with = FALSE]
    setnames(bd, c("sample", if (!is.null(lo_col)) "lo", "hi"))
    bd <- melt(bd, id.vars = "sample", value.name = "y")[!is.na(y)]
    if (logy) bd <- bd[y > 0]
    p <- p + geom_crossbar(data = bd, inherit.aes = FALSE,
                           aes(x = sample, y = y, ymin = y, ymax = y),
                           colour = "red", linewidth = 0.3, width = 0.8)
  }
  if (logy) ylab <- paste0(ylab, "  (log10 axis)")
  p <- p + labs(x = NULL, y = ylab,
                title = paste0(metric, "  (red = per-animal mean +/- ", SD_N,
                               " sd, log10 scale, doublets already removed)")) + rot
  if (logy) p <- p + scale_y_log10()
  ggsave(file.path(OUTDIR, paste0("qc_violin_", metric, ".png")), p, width = 10, height = 5, dpi = 150)
}
viol("nFeature_mouse", "unique genes per nucleus", TRUE, "nFeature_lo", "nFeature_hi")
viol("nCount_mouse",   "UMIs per nucleus",         TRUE, "nCount_lo",   "nCount_hi")

sub <- qc[sample(.N, min(SCATTER_N, .N))]
fits <- qc[, { m <- stats::lm(log10(pmax(nFeature_mouse, 1)) ~ log10(pmax(nCount_mouse, 1)))
               cf <- stats::coef(m)
  .(intercept = unname(cf[1]), slope = unname(cf[2]), r2 = summary(m)$r.squared, n = .N) }, by = sample]
fits[, lab := sprintf("log10(genes) = %.3f x log10(UMIs) %+.3f\nR2 = %.4f  n = %s",
                      slope, intercept, r2, format(n, big.mark = ","))]
fwrite(fits, file.path(OUTDIR, "qc_count_vs_feature_fits.csv"))
msg("genes-vs-UMIs log-log slope: ", round(min(fits$slope), 3), " - ", round(max(fits$slope), 3),
    ", R2 ", round(min(fits$r2), 4), " - ", round(max(fits$r2), 4))

# y = Inf, not -Inf: the panel is log10 and log10(-Inf) is NaN, which silently
# drops every label ("Removed 12 rows containing missing values").
ggsave(file.path(OUTDIR, "qc_scatter_count_vs_feature.png"),
       ggplot(sub, aes(nCount_mouse, nFeature_mouse)) +
         geom_point(aes(colour = percent_mt), size = 0.15, alpha = 0.3) +
         geom_abline(data = fits, aes(intercept = intercept, slope = slope),
                     colour = "red", linewidth = 0.4) +
         geom_text(data = fits, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
                   hjust = 1.03, vjust = 1.35, size = 2.4, lineheight = 0.95, colour = "grey20") +
         scale_x_log10() + scale_y_log10() +
         scale_colour_viridis_c(option = "magma", limits = c(0, 5), oob = scales::squish) +
         facet_wrap(~ sample, ncol = 4) +
         labs(x = "UMIs per nucleus, mouse genes  (log10 axis)",
              y = "unique genes per nucleus  (log10 axis)", colour = "% mt",
              title = "Unique genes vs total UMIs per nucleus, by animal",
              subtitle = "ambient-corrected, doublet-free; red = per-animal least-squares fit on all nuclei"),
       width = 13, height = 9.5, dpi = 150)

mt_ecdf <- ggplot(qc, aes(percent_mt, colour = sample)) + stat_ecdf(linewidth = 0.4) +
  coord_cartesian(xlim = c(0, 25)) + labs(x = "% mitochondrial UMIs", y = "cumulative fraction")
if (!is.na(MT_MAX)) mt_ecdf <- mt_ecdf +
  geom_vline(xintercept = MT_MAX, colour = "red", linetype = 2, linewidth = 0.4) +
  labs(subtitle = paste0("red = MT_MAX ceiling at ", MT_MAX, "% (chosen pre-CellBender)"))
ggsave(file.path(OUTDIR, "qc_ecdf_percent_mt.png"), mt_ecdf, width = 8, height = 5, dpi = 150)

# x capped at the pooled 99th percentile. 0.999 is not extreme enough at this n:
# it is still the ~300th worst nucleus and on a heavy tail caps nothing.
MT_HIST_Q <- 0.99
mt_cap <- max(ceiling(as.numeric(quantile(qc$percent_mt, MT_HIST_Q))), 5)
n_beyond <- qc[, .(n = sum(percent_mt > mt_cap), total = .N, med = median(percent_mt)), by = sample]
n_beyond[, lab := paste0("median ", sprintf("%.3f", med), "%\nbeyond ", mt_cap, "%: n = ", n, " of ", total)]
ggsave(file.path(OUTDIR, "qc_hist_percent_mt.png"),
       ggplot(qc[percent_mt <= mt_cap], aes(percent_mt, fill = diet)) +
         geom_histogram(binwidth = mt_cap / 100, colour = NA) +
         geom_vline(data = n_beyond, aes(xintercept = med), colour = "red",
                    linetype = 2, linewidth = 0.35) +
         geom_text(data = n_beyond, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
                   hjust = 1.05, vjust = 1.3, size = 2.6, lineheight = 0.95, colour = "grey25") +
         facet_wrap(~ sample, ncol = 4, scales = "free_y") + scale_y_log10() +
         labs(x = "mitochondrial UMIs as % of the nucleus' mouse-gene UMIs",
              y = "nuclei  (log10 axis)",
              title = "Mitochondrial fraction per nucleus, by animal -- AFTER ambient removal",
              subtitle = paste0("x capped at the pooled ", MT_HIST_Q * 100, "th percentile = ", mt_cap,
                                "%; red dashed = that animal's median.\n",
                                "snRNA-seq: these reads are ambient, and CellBender has already removed ",
                                "ambient, so re-read MT_MAX (currently ", MT_MAX, "%) off this plot.")),
       width = 12, height = 8, dpi = 150)

## ---- FILTER AND SAVE ------------------------------------------------------

if (APPLY_FILTER) {
  before <- ncol(obj)
  obj <- subset(obj, subset = qc_pass)
  msg("filter applied per animal at ", SD_N, " ", DISPERSION, " (log10), nFeature >= ",
      NFEATURE_MIN, ", percent_mt <= ", MT_MAX, "%: ",
      before, " -> ", ncol(obj), " nuclei (", round(100 * ncol(obj) / before, 2), "%)")

  drop <- intersect(c("nFeature_lo", "nFeature_hi", "nCount_lo", "nCount_hi",
                      "percent_mt_hi", "qc_pass"), colnames(obj@meta.data))
  for (d in drop) obj@meta.data[[d]] <- NULL
  msg("clean copy metadata: ", paste(colnames(obj@meta.data), collapse = ", "))

  saveRDS(obj, CLEANRDS, compress = FALSE)
  msg("clean filtered object written: ", CLEANRDS, " (",
      round(file.size(CLEANRDS) / 2^30, 2), " GB, ", ncol(obj), " nuclei, ", nrow(obj), " features)")

  kept <- as.data.table(obj@meta.data)[, .(kept = .N), by = sample]
  fwrite(merge(retained[, .(sample, total, kept_flag = kept)], kept, by = "sample"),
         file.path(OUTDIR, "filtered_object_composition.csv"))
  print(kept[order(sample)])
} else {
  msg("APPLY_FILTER is FALSE - no clean copy written")
}

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
