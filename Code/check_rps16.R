#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# check_rps16.R  --  diagnostic, not part of the numbered pipeline. Run by
# check_rps16.sh.
#
# 10_3_deg.R returned exactly one diet-associated gene across all
# eight Azimuth subclasses: Rps16 in L5 ET, model B, dietHFD +4.20 (FDR 4.7e-4)
# with a near-equal-and-opposite interaction dietHFD:sexM -4.41 (FDR 1.8e-3).
# A ribosomal protein gene moving 4 log2 units in one subclass and nowhere else
# is not a credible biological diet effect, so this script asks what produced
# it. Four questions, in order of how decisive they are:
#
#   1. WHAT DOES THE GENE ACTUALLY DO ACROSS THE 12 LIBRARIES. Raw count, CPM
#      and per-nucleus detection rate per animal in L5 ET. If one animal is an
#      order of magnitude off the rest, the coefficient is that animal.
#   2. LEAVE ONE ANIMAL OUT. Refit model B twelve times, dropping one animal
#      each time, and report the dietHFD logFC and FDR for Rps16. A real effect
#      is stable; a one-animal artefact collapses when that animal leaves.
#   3. IS IT THE GENE OR THE ANIMAL. Same per-animal CPM for every other
#      ribosomal protein gene (^Rp[sl]) in L5 ET. If the odd animal is odd for
#      all of them it is a library-composition effect; if only Rps16 moves it is
#      specific to this gene's counting.
#   4. IS IT L5 ET OR EVERYWHERE. Rps16 CPM per animal in all eight subclasses.
#
# It also writes the pseudobulk count matrix that 10_3 built and discarded, so
# any further question about a specific gene can be answered without reloading
# the 11.7 GB object.
#
# The pseudobulk construction is copied from 10_3_deg.R deliberately
# and must stay identical to it -- same barcode join, same non-glutamatergic
# drop, same transgene removal -- or the numbers here would not be the numbers
# being explained.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(data.table)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_sub_excitatory.rds")
AZDIR   <- file.path(PROJECT, "Results", "10_2_excitatory_azimuth")
DEGDIR  <- file.path(PROJECT, "Results", "10_3_excitatory_deg")
OUTDIR  <- file.path(PROJECT, "Results", "check_rps16")

GROUP      <- "predicted.subclass"
GENE       <- "Rps16"
FOCUS      <- "L5 ET"
TRANSGENES <- gsub("_", "-", c("EGFP", "tdTomato", "EBFP2", "AAV_shared"), fixed = TRUE)
MIN_NUCLEI <- 10
LAYER_ORDER <- c("L2/3 IT", "L4/5 IT", "L5 IT", "L5 ET", "L5/6 NP",
                 "L6 IT", "L6 IT Car3", "L6 CT", "L6b")
RES  <- 0.2
PRED <- file.path(AZDIR, sprintf("azimuth_predictions_res%s.rds", RES))

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---- rebuild exactly what 10_3 built ----------------------------------------

msg("reading ", basename(IN))
obj <- readRDS(IN)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]], search = "counts")) > 1) obj <- JoinLayers(obj)

pred <- as.data.table(readRDS(PRED))
md <- as.data.table(obj@meta.data, keep.rownames = "barcode")
p  <- pred[, .(barcode, group = get(GROUP), pred_class = predicted.class)]
md[p, on = "barcode", `:=`(group = i.group, pred_class = i.pred_class)]
stopifnot(identical(md$barcode, colnames(obj)), !anyNA(md$group))

modal <- md[, .N, by = .(group, pred_class)][order(group, -N)][, .SD[1], by = group]
glut  <- modal[pred_class == "Glutamatergic", as.character(group)]
md    <- md[group %in% glut]
sizes <- md[, .N, by = group][order(-N)]
levs  <- c(intersect(LAYER_ORDER, sizes$group), setdiff(sizes$group, LAYER_ORDER))
md[, group := factor(as.character(group), levels = levs)]
md[, diet := factor(as.character(diet), levels = c("SD", "HFD"))]
md[, sex  := factor(as.character(sex),  levels = c("F", "M"))]
md[, sample := as.character(sample)]
msg(nlevels(md$group), " subclasses, ", format(nrow(md), big.mark = ","), " nuclei")

cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")[, md$barcode, drop = FALSE]
cnt <- cnt[setdiff(rownames(cnt), TRANSGENES), , drop = FALSE]
rm(obj); invisible(gc())
if (!GENE %in% rownames(cnt)) stop(GENE, " is not a feature of the object")

## per-nucleus detection, before any summing -----------------------------------
gv <- as.numeric(cnt[GENE, ])
md[, gene_umi := gv]
cell_tab <- md[, .(n_nuclei = .N,
                   pct_detected = round(100 * mean(gene_umi > 0), 2),
                   mean_umi = round(mean(gene_umi), 3),
                   max_umi = max(gene_umi)),
               by = .(group, sample, mouse, diet, sex)]
setorder(cell_tab, group, diet, sex, sample)
fwrite(cell_tab, file.path(OUTDIR, sprintf("%s_per_nucleus_detection.csv", GENE)))

key <- factor(paste(md$group, md$sample, sep = "__"))
ind <- Matrix::sparse.model.matrix(~ 0 + key)
colnames(ind) <- levels(key)
pb <- round(as.matrix(cnt %*% ind))
rm(cnt, ind); invisible(gc())
msg("pseudobulk: ", nrow(pb), " genes x ", ncol(pb), " libraries")

# the matrix 10_3 discarded -- written once, here, so no later question needs
# the 11.7 GB object again
fwrite(as.data.table(pb, keep.rownames = "gene"),
       file.path(DEGDIR, "pseudobulk_counts.csv.gz"), compress = "gzip")

lib <- data.table(library = colnames(pb))
lib[, c("group", "sample") := tstrsplit(library, "__", fixed = TRUE)]
lib <- merge(lib, unique(md[, .(sample, mouse, diet, sex)]), by = "sample", sort = FALSE)
lib <- merge(lib, md[, .N, by = .(group, sample)][, group := as.character(group)],
             by = c("group", "sample"), all.x = TRUE)
setnames(lib, "N", "n_nuclei")
lib[, lib_size := colSums(pb)[library]]
lib[, group := factor(group, levels = levs)]
setorder(lib, group, diet, sex, sample)

cpm_all <- edgeR::cpm(pb, log = FALSE)

## ---- Q1: the gene across the 12 L5 ET libraries -----------------------------

focus <- lib[group == FOCUS]
focus[, raw_count := pb[GENE, library]]
focus[, cpm := round(cpm_all[GENE, library], 2)]
focus <- merge(focus, cell_tab[group == FOCUS, .(sample, pct_detected, mean_umi, max_umi)],
               by = "sample", sort = FALSE)
setorder(focus, diet, sex, sample)
fwrite(focus, file.path(OUTDIR, sprintf("%s_%s_per_library.csv", GENE,
                                        gsub("[ /]", "", FOCUS))))
msg("--- Q1: ", GENE, " in ", FOCUS, ", one row per animal ---")
print(focus[, .(mouse, diet, sex, n_nuclei, lib_size, raw_count, cpm,
                pct_detected, mean_umi, max_umi)])

## ---- Q2: leave one animal out -----------------------------------------------

fit_one <- function(l) {
  des <- model.matrix(~ diet * sex, data = l)
  if (qr(des)$rank < ncol(des) || nrow(l) - ncol(des) < 2) return(NULL)
  y <- DGEList(counts = pb[, l$library, drop = FALSE])
  keep <- filterByExpr(y, design = des)
  if (!keep[GENE]) return(data.table(kept = FALSE))
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y); y <- estimateDisp(y, des, robust = TRUE)
  fit <- glmQLFit(y, des, robust = TRUE)
  out <- lapply(c("dietHFD", "dietHFD:sexM"), function(cf) {
    tt <- topTags(glmQLFTest(fit, coef = cf), n = Inf, sort.by = "none")$table
    data.table(contrast = cf, logFC = round(tt[GENE, "logFC"], 2),
               PValue = signif(tt[GENE, "PValue"], 3),
               FDR = signif(tt[GENE, "FDR"], 3))
  })
  cbind(data.table(kept = TRUE), rbindlist(out))
}

l0 <- lib[group == FOCUS & n_nuclei >= MIN_NUCLEI]
full_fit <- fit_one(copy(l0))[, dropped := "none (all 12)"]
loo <- rbindlist(lapply(l0$mouse, function(m) {
  r <- fit_one(copy(l0[mouse != m]))
  if (is.null(r)) return(data.table(dropped = m, kept = NA))
  r[, dropped := m][]
}), fill = TRUE)
loo <- rbind(full_fit, loo, fill = TRUE)
loo <- merge(loo, unique(l0[, .(dropped = mouse, diet, sex)]), by = "dropped",
             all.x = TRUE, sort = FALSE)
fwrite(loo, file.path(OUTDIR, sprintf("%s_leave_one_out.csv", GENE)))
msg("--- Q2: leave-one-animal-out, model B, ", GENE, " in ", FOCUS, " ---")
print(loo[order(contrast, FDR)])

## ---- Q3: is it the gene or the animal ---------------------------------------

ribo <- grep("^Rp[sl][0-9]", rownames(pb), value = TRUE)
ribo <- setdiff(ribo, grep("-ps|Rpl21-|rt$", ribo, value = TRUE))
rc <- as.data.table(cpm_all[ribo, focus$library, drop = FALSE], keep.rownames = "gene")
rl <- melt(rc, id.vars = "gene", variable.name = "library", value.name = "cpm")
rl <- merge(rl, focus[, .(library, mouse, diet, sex)], by = "library")
# per gene, how far is each animal from that gene's median across animals
rl[, z := (cpm - median(cpm)) / (mad(cpm) + 1e-9), by = gene]
ribo_summ <- rl[, .(median_abs_z = round(median(abs(z)), 2),
                    max_abs_z = round(max(abs(z)), 2)), by = mouse][order(-median_abs_z)]
gene_z <- rl[gene == GENE][order(-abs(z)),
                           .(mouse, diet, sex, cpm = round(cpm, 1), z = round(z, 1))]
fwrite(rl, file.path(OUTDIR, "ribosomal_cpm_L5ET_long.csv"))
fwrite(ribo_summ, file.path(OUTDIR, "ribosomal_animal_outlierness.csv"))
msg("--- Q3a: across ", length(ribo), " ribosomal genes in ", FOCUS,
    ", which animal is the outlier ---")
print(ribo_summ)
msg("--- Q3b: ", GENE, " itself, per animal, robust z within the gene ---")
print(gene_z)

## ---- Q4: the same gene in the other subclasses ------------------------------

wide <- lib[, .(library, group, mouse, diet, sex)]
wide[, cpm := round(cpm_all[GENE, library], 1)]
tab <- dcast(wide, mouse + diet + sex ~ group, value.var = "cpm")
setorder(tab, diet, sex, mouse)
fwrite(tab, file.path(OUTDIR, sprintf("%s_cpm_by_subclass.csv", GENE)))
msg("--- Q4: ", GENE, " CPM per animal in every subclass ---")
print(tab)

## ---- figure ------------------------------------------------------------------

wide[, group := factor(group, levels = levs)]
pl <- ggplot(wide, aes(diet, cpm, colour = sex)) +
  geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.2) +
  facet_wrap(~ group, nrow = 2, scales = "free_y") +
  scale_colour_manual(values = c(F = "#bb292c", M = "#376795")) +
  labs(x = NULL, y = paste0(GENE, "  pseudobulk CPM"),
       title = paste0(GENE, " per animal, by subclass"),
       subtitle = "one point = one animal; the DE hit was L5 ET") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(file.path(OUTDIR, sprintf("%s_cpm_by_subclass.png", GENE)), pl,
       width = 11, height = 6, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
