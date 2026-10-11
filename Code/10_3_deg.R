#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10_3_deg.R  --  pseudobulk differential expression inside a neuronal subset,
# grouped by the Azimuth transferred subclass. Run by 10_3_deg.sh.
#
# The lineage is the first argument: `excitatory` (default) or `inhibitory`.
# It selects the subset object, the Azimuth predictions written by 10_2, and
# which reference class is kept -- Glutamatergic for the excitatory lineage,
# GABAergic for the inhibitory one.
#
# This is step 9 re-run one level down. Step 9 tested HFD vs SD within each of
# the 13 `annot` cell types on the whole object; here the unit is the Azimuth
# `predicted.subclass` call inside one neuronal lineage, so the question becomes
# which cortical subclass carries the diet effect.
#
# Everything about the statistics is unchanged from 9_pseudobulk_deg.R:
# counts summed to one library per ANIMAL per group, edgeR quasi-likelihood
# F-tests, two models
#
#   A  ~ diet              HFD vs SD, sex ignored
#   B  ~ diet * sex        diet, sex, and the diet:sex interaction
#
# QL on animal-level libraries rather than a per-nucleus test because the unit
# of replication is the animal, n = 12 -- see [[insular-animal-level-stats]].
# With ~ 50k L5 IT nuclei a per-nucleus test would return thousands of "DE"
# genes that are one mouse. The inhibitory lineage is ~ 4x smaller, so its
# subclasses have fewer nuclei per animal but exactly the same 12 replicates --
# the power that matters here is n = 12 either way, and what changes between
# lineages is the precision of each library, not the sample size.
#
# READING MODEL B. R uses treatment coding, so with ~ diet * sex the `dietHFD`
# coefficient is HFD vs SD IN FEMALES (the reference level), not an average over
# sexes, and `dietHFD:sexM` asks whether the diet effect differs between sexes.
# Model A's `dietHFD` is the sex-averaged effect but pays for it by leaving sex
# variance in the residual. Neither is a corrected version of the other.
#
# WHERE THE GROUPS COME FROM. The per-nucleus Azimuth predictions written by
# 10_2_excitatory_azimuth.R, joined onto the excitatory object by barcode. The
# label is per NUCLEUS, not per cluster: nothing about RNA_snn_res.0.2 is used
# here, so a cluster that Azimuth splits across two subclasses contributes to
# both. No prediction-score or mapping-score filter is applied (user, 2026-09-07)
# -- the median predicted.subclass score is 1.0 in almost every cluster, and a
# mapping-score gate would remove the L6 IT Car3 and small-cluster nuclei
# preferentially rather than at random.
#
# Azimuth also assigns a handful of nuclei to subclasses of the other lineages
# (in the excitatory subset: Micro-PVM, VLMC, Vip, ...; ~850 of 137k). Those are
# dropped: a subclass is kept only when its nuclei are modally the lineage's own
# predicted.class. They are counted in the log rather than silently discarded.
#
# SANITY CHECK BUILT IN: the `sexM` coefficient in model B must recover Xist
# (down in males) and the Y genes Ddx3y, Uty, Eif2s3y, Kdm5d (up in males) in
# any subclass with real power. If they are absent, the pseudobulk matrix or the
# animal metadata is wrong and nothing else in the output should be believed.
#
# The 4 transgene features are dropped before DE: they are the experimental
# label, not biology, and their counts are presence/absence calls.
#
# A subclass is skipped when fewer than MIN_ANIMALS animals clear MIN_NUCLEI
# nuclei, or when the surviving animals cannot support the design -- both are
# reported rather than silently producing a rank-deficient fit.
#
# Usage:
#   sbatch 10_3_deg.sh                     # excitatory, res 0.2 predictions
#   sbatch 10_3_deg.sh inhibitory 0.2
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

PROJECT <- "/scratch/user/uqschauq/Insular"

# Per lineage: the reference class its subclasses must be modally assigned to,
# and the order its subclasses are reported in. Cortical depth for the
# excitatory lineage; the standard MOp GABAergic order (which is developmental
# origin -- MGE-derived Pvalb/Sst first, then CGE-derived Vip/Lamp5/Sncg) for
# the inhibitory one. Anything Azimuth returns that is not listed is appended by
# size, so a new label cannot silently vanish from the output.
LINEAGE_SPEC <- list(
  excitatory = list(
    class = "Glutamatergic",
    order = c("L2/3 IT", "L4/5 IT", "L5 IT", "L5 ET", "L5/6 NP",
              "L6 IT", "L6 IT Car3", "L6 CT", "L6b")),
  inhibitory = list(
    class = "GABAergic",
    order = c("Pvalb", "Sst", "Sst Chodl", "Vip", "Lamp5", "Sncg", "Meis2")))

GROUP       <- "predicted.subclass"
TRANSGENES  <- gsub("_", "-", c("EGFP", "tdTomato", "EBFP2", "AAV_shared"), fixed = TRUE)
MIN_NUCLEI  <- 10     # nuclei an animal needs to contribute a pseudobulk library
MIN_ANIMALS <- 8      # animals a subclass needs before it is tested at all
FDR_CUT     <- 0.05
TOP_LABEL   <- 5      # genes labelled per volcano facet
SEX_GENES   <- c("Xist", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d")
SEED        <- 42

args    <- commandArgs(trailingOnly = TRUE)
LINEAGE <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "excitatory"
RES     <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(args[2]) else 0.2
if (!LINEAGE %in% names(LINEAGE_SPEC))
  stop("first argument must be one of: ", paste(names(LINEAGE_SPEC), collapse = ", "))
if (!is.finite(RES)) stop("second argument must be a resolution, e.g. 0.2")

KEEP_CLASS  <- LINEAGE_SPEC[[LINEAGE]]$class
LABEL_ORDER <- LINEAGE_SPEC[[LINEAGE]]$order
IN     <- file.path(PROJECT, "Data", "Seurat", sprintf("HFDseq_v2_sub_%s.rds", LINEAGE))
AZDIR  <- file.path(PROJECT, "Results", sprintf("10_2_%s_azimuth", LINEAGE))
OUTDIR <- file.path(PROJECT, "Results", sprintf("10_3_%s_deg", LINEAGE))
PRED   <- file.path(AZDIR, sprintf("azimuth_predictions_res%s.rds", RES))

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

if (!file.exists(IN))  stop("no ", LINEAGE, " subset at ", IN)
if (!file.exists(PRED))
  stop("no Azimuth predictions at ", PRED, " -- run 10_2_azimuth.sh ", LINEAGE, " first")

## ---- object + transferred labels --------------------------------------------

msg("reading ", basename(IN), " (", round(file.size(IN) / 2^30, 1), " GB)")
obj <- readRDS(IN)
DefaultAssay(obj) <- "RNA"
if (length(Layers(obj[["RNA"]], search = "counts")) > 1) obj <- JoinLayers(obj)
stopifnot(all(c("sample", "mouse", "diet", "sex") %in% names(obj@meta.data)))
msg("loaded ", format(ncol(obj), big.mark = ","), " ", LINEAGE, " nuclei")

pred <- as.data.table(readRDS(PRED))
if (!GROUP %in% names(pred)) stop("no '", GROUP, "' column in ", basename(PRED))
msg("predictions: ", format(nrow(pred), big.mark = ","), " nuclei, ",
    uniqueN(pred[[GROUP]]), " ", GROUP, " labels")

md <- as.data.table(obj@meta.data, keep.rownames = "barcode")
miss <- setdiff(md$barcode, pred$barcode)
if (length(miss))
  stop(length(miss), " object barcodes are absent from ", basename(PRED),
       " -- the predictions were made on a different object")
# An update join, not merge(): it writes into md in place and so cannot reorder
# the rows away from the column order of the count matrix.
p <- pred[, .(barcode, group = get(GROUP), pred_class = predicted.class,
              diet_pred = diet, sex_pred = sex)]
md[p, on = "barcode", `:=`(group = i.group, pred_class = i.pred_class,
                           diet_pred = i.diet_pred, sex_pred = i.sex_pred)]
stopifnot(identical(md$barcode, colnames(obj)), !anyNA(md$group))

# The prediction table carries its own copy of the animal metadata. It must
# agree with the object's, or the two were built from different subsets and the
# design matrix would be attached to the wrong libraries.
if (!identical(as.character(md$diet), as.character(md$diet_pred)) ||
    !identical(as.character(md$sex),  as.character(md$sex_pred)))
  stop("diet/sex disagree between the object and ", basename(PRED))
md[, c("diet_pred", "sex_pred") := NULL]

## ---- keep the subclasses of this lineage's own class ------------------------

cls <- md[, .N, by = .(group, pred_class)][order(group, -N)]
modal <- cls[, .SD[1], by = group]
keepers <- modal[pred_class == KEEP_CLASS, as.character(group)]
if (!length(keepers))
  stop("no subclass is modally ", KEEP_CLASS, " -- wrong lineage or wrong predictions file")
drop <- md[!group %in% keepers, .N, by = group][order(-N)]
if (nrow(drop))
  msg("dropping ", format(sum(drop$N), big.mark = ","), " nuclei in ", nrow(drop),
      " subclass(es) not modally ", KEEP_CLASS, ": ",
      paste(sprintf("%s (%d)", drop$group, drop$N), collapse = ", "))
fwrite(cls, file.path(OUTDIR, "subclass_class_composition.csv"))
md <- md[group %in% keepers]

sizes  <- md[, .N, by = group][order(-N)]
levs   <- c(intersect(LABEL_ORDER, sizes$group), setdiff(sizes$group, LABEL_ORDER))
md[, group := factor(as.character(group), levels = levs)]
md[, diet := factor(as.character(diet), levels = c("SD", "HFD"))]
md[, sex  := factor(as.character(sex),  levels = c("F", "M"))]
md[, sample := as.character(sample)]
msg("testing ", nlevels(md$group), " subclasses over ",
    format(nrow(md), big.mark = ","), " nuclei and ", uniqueN(md$sample), " animals")

## ---- pseudobulk: sum counts per animal per subclass -------------------------

cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
cnt <- cnt[, md$barcode, drop = FALSE]
cnt <- cnt[setdiff(rownames(cnt), TRANSGENES), , drop = FALSE]
rm(obj); invisible(gc())

key <- factor(paste(md$group, md$sample, sep = "__"))
ind <- Matrix::sparse.model.matrix(~ 0 + key)
colnames(ind) <- levels(key)
pb <- as.matrix(cnt %*% ind)
rm(cnt, ind); invisible(gc())

# CellBender writes integer counts but stores them as doubles; edgeR wants whole
# numbers, so assert rather than assume before rounding the storage type away.
off <- max(abs(pb - round(pb)))
if (off > 1e-6) stop("pseudobulk counts are not whole numbers (max deviation ", off, ")")
pb <- round(pb)
msg("pseudobulk: ", nrow(pb), " genes x ", ncol(pb), " libraries")

lib <- data.table(library = colnames(pb))
lib[, c("group", "sample") := tstrsplit(library, "__", fixed = TRUE)]
lib <- merge(lib, unique(md[, .(sample, mouse, diet, sex)]), by = "sample", sort = FALSE)
ncell <- md[, .N, by = .(group, sample)]
ncell[, group := as.character(group)]
lib <- merge(lib, ncell, by = c("group", "sample"), all.x = TRUE)
setnames(lib, "N", "n_nuclei")
lib[is.na(n_nuclei), n_nuclei := 0L]
lib[, lib_size := colSums(pb)[library]]
lib[, group := factor(group, levels = levs)]
setorder(lib, group, diet, sex, sample)
fwrite(lib, file.path(OUTDIR, "pseudobulk_library_qc.csv"))
msg("nuclei per library: median ", median(lib$n_nuclei),
    ", min ", min(lib$n_nuclei), ", ", lib[n_nuclei < MIN_NUCLEI, .N],
    " below MIN_NUCLEI = ", MIN_NUCLEI)

## ---- DE, per subclass -------------------------------------------------------

MODELS <- list(
  A = list(formula = ~ diet,       coefs = "dietHFD"),
  B = list(formula = ~ diet * sex, coefs = c("dietHFD", "sexM", "dietHFD:sexM")))

run_one <- function(ct, model_name) {
  spec <- MODELS[[model_name]]
  l <- lib[group == ct & n_nuclei >= MIN_NUCLEI]
  if (nrow(l) < MIN_ANIMALS) {
    msg("  ", ct, " / ", model_name, ": SKIPPED, only ", nrow(l),
        " animal(s) with >= ", MIN_NUCLEI, " nuclei")
    return(NULL)
  }
  l[, diet := droplevels(diet)][, sex := droplevels(sex)]
  if (nlevels(l$diet) < 2) {
    msg("  ", ct, " / ", model_name, ": SKIPPED, diet not represented in both groups")
    return(NULL)
  }
  des <- model.matrix(spec$formula, data = l)
  if (qr(des)$rank < ncol(des)) {
    msg("  ", ct, " / ", model_name, ": SKIPPED, design is rank-deficient (",
        paste(sprintf("%s=%d", levels(interaction(l$diet, l$sex)),
                      table(interaction(l$diet, l$sex))), collapse = " "), ")")
    return(NULL)
  }
  if (nrow(l) - ncol(des) < 2) {
    msg("  ", ct, " / ", model_name, ": SKIPPED, only ", nrow(l) - ncol(des),
        " residual df")
    return(NULL)
  }
  missing_cf <- setdiff(spec$coefs, colnames(des))
  if (length(missing_cf))
    stop(ct, " / ", model_name, ": coefficient(s) absent from the design: ",
         paste(missing_cf, collapse = ", "))

  y <- DGEList(counts = pb[, l$library, drop = FALSE], samples = as.data.frame(l))
  keep <- filterByExpr(y, design = des)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  y <- estimateDisp(y, des, robust = TRUE)
  fit <- glmQLFit(y, des, robust = TRUE)

  out <- rbindlist(lapply(spec$coefs, function(cf) {
    tt <- topTags(glmQLFTest(fit, coef = cf), n = Inf, sort.by = "none")$table
    data.table(subclass = ct, model = model_name, contrast = cf,
               gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM,
               PValue = tt$PValue, FDR = tt$FDR)
  }))
  attr(out, "meta") <- data.table(
    subclass = ct, model = model_name, n_libraries = nrow(l),
    n_genes_tested = sum(keep), residual_df = nrow(l) - ncol(des),
    common_bcv = round(sqrt(y$common.dispersion), 3),
    n_SD = sum(l$diet == "SD"), n_HFD = sum(l$diet == "HFD"),
    min_nuclei = min(l$n_nuclei), median_nuclei = as.numeric(median(l$n_nuclei)))
  out
}

res <- list(); meta <- list()
for (mn in names(MODELS)) {
  msg("model ", mn, ": ", deparse(MODELS[[mn]]$formula))
  for (ct in levs) {
    r <- run_one(ct, mn)
    if (is.null(r)) next
    res[[length(res) + 1]]   <- r
    meta[[length(meta) + 1]] <- attr(r, "meta")
  }
}
if (!length(res)) stop("no subclass produced a fit")
de   <- rbindlist(res)
meta <- rbindlist(meta)
de[, subclass := factor(subclass, levels = levs)]
meta[, subclass := factor(subclass, levels = levs)]

fwrite(de[order(model, contrast, subclass, PValue)],
       file.path(OUTDIR, "pseudobulk_de_all.csv.gz"), compress = "gzip")
fwrite(meta[order(model, subclass)], file.path(OUTDIR, "pseudobulk_fit_summary.csv"))

sig <- de[FDR < FDR_CUT, .(n_de = .N,
                           n_up = sum(logFC > 0), n_down = sum(logFC < 0)),
          by = .(model, contrast, subclass)]
full <- CJ(model = unique(de$model), contrast = unique(de$contrast),
           subclass = factor(levs, levels = levs), unique = TRUE)
full <- full[unique(de[, .(model, contrast)]), on = .(model, contrast)]
sig <- sig[full, on = .(model, contrast, subclass)]
sig[is.na(n_de), `:=`(n_de = 0L, n_up = 0L, n_down = 0L)]
sig[, label := factor(fifelse(model == "A", "A: ~ diet",
                       fifelse(contrast == "dietHFD",     "B: diet (in females)",
                       fifelse(contrast == "sexM",        "B: sex (in SD)",
                                                          "B: diet x sex"))),
                      levels = c("A: ~ diet", "B: diet (in females)",
                                 "B: sex (in SD)", "B: diet x sex"))]
fwrite(sig[order(label, subclass)], file.path(OUTDIR, "pseudobulk_de_counts.csv"))
msg("DE genes at FDR < ", FDR_CUT, ":")
print(dcast(sig, subclass ~ label, value.var = "n_de"))

msg("sex-gene sanity check (model B, sexM coefficient):")
sx <- de[model == "B" & contrast == "sexM" & gene %in% SEX_GENES,
         .(subclass, gene, logFC = round(logFC, 2), FDR = signif(FDR, 2))]
print(dcast(sx, subclass ~ gene, value.var = "logFC"))
hit <- sx[FDR < FDR_CUT, uniqueN(subclass)]
if (!nrow(sx) || hit == 0)
  msg("WARNING: no sex gene reaches FDR < ", FDR_CUT, " in any subclass -- ",
      "check the pseudobulk matrix and the animal metadata before reading anything else")

## ---- figures ----------------------------------------------------------------

CONTRAST_COLS <- c("A: ~ diet" = "#dd5129", "B: diet (in females)" = "#0f7ba2",
                   "B: sex (in SD)" = "#43b284", "B: diet x sex" = "#fab255")

pA <- ggplot(sig, aes(subclass, n_de, fill = label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = CONTRAST_COLS, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = paste0("genes at FDR < ", FDR_CUT)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())

volcano <- function(dd, xlab) {
  dd <- copy(dd)
  dd[, sig := FDR < FDR_CUT]
  lab <- dd[sig == TRUE][order(FDR), head(.SD, TOP_LABEL), by = subclass]
  p <- ggplot(dd, aes(logFC, -log10(PValue))) +
    geom_point(aes(colour = sig), size = 0.5, alpha = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#bb292c"),
                        name = NULL, labels = c("n.s.", paste0("FDR < ", FDR_CUT))) +
    facet_wrap(~ subclass, ncol = 4, scales = "free_y") +
    labs(x = xlab, y = "-log10 p") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 8, face = "bold"))
  if (nrow(lab) && requireNamespace("ggrepel", quietly = TRUE))
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 2.4,
                                      max.overlaps = 20, segment.size = 0.2)
  p
}

vA <- volcano(de[model == "A" & contrast == "dietHFD"],
              "log2FC, HFD vs SD   (model A: ~ diet)")
fig <- pA / vA + plot_layout(heights = c(1, 2.1))
ggsave(file.path(OUTDIR, "pseudobulk_overview.png"), fig, width = 13, height = 12, dpi = 300)
ggsave(file.path(OUTDIR, "pseudobulk_overview.pdf"), fig, width = 13, height = 12)
msg("written pseudobulk_overview.png/.pdf")

vI <- volcano(de[model == "B" & contrast == "dietHFD:sexM"],
              "log2FC of the diet effect in males vs females   (model B: ~ diet * sex)")
ggsave(file.path(OUTDIR, "pseudobulk_volcano_interaction.png"), vI,
       width = 13, height = 8, dpi = 300)
ggsave(file.path(OUTDIR, "pseudobulk_volcano_interaction.pdf"), vI, width = 13, height = 8)

vS <- volcano(de[model == "B" & contrast == "sexM"],
              "log2FC, M vs F   (model B: ~ diet * sex, sex coefficient)")
ggsave(file.path(OUTDIR, "pseudobulk_volcano_sex.png"), vS,
       width = 13, height = 8, dpi = 300)
msg("written interaction and sex volcanoes")

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
