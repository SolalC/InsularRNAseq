#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 9_pseudobulk_deg.R  --  pseudobulk differential expression, run by
# 9_pseudobulk_deg.sh.
#
# One DE analysis per cell type, on counts summed to one library per ANIMAL per
# cell type (13 x 12 = 156 pseudobulk libraries). Two models, as asked:
#
#   A  ~ diet              HFD vs SD, sex ignored
#   B  ~ diet * sex        diet, sex, and the diet:sex interaction
#
# Both are edgeR quasi-likelihood F-tests (glmQLFit + glmQLFTest, robust). QL
# rather than LRT or a per-nucleus Wilcoxon because the unit of replication here
# is the animal, n = 12, and the QL F-test is the one that keeps type I error
# under control at that n -- see [[insular-animal-level-stats]]: a per-nucleus
# test would treat one mouse's ~20k nuclei as 20k replicates and return
# thousands of "DE" genes that are one animal.
#
# READING MODEL B. R uses treatment coding, so with ~ diet * sex the `dietHFD`
# coefficient is HFD vs SD IN FEMALES (the reference level), not an average over
# sexes, and `dietHFD:sexM` asks whether the diet effect differs between sexes.
# Model A's `dietHFD` is the sex-averaged effect but pays for it by leaving sex
# variance in the residual. They answer different questions; neither is a
# corrected version of the other.
#
# SANITY CHECK BUILT IN: the `sexM` coefficient in model B must recover Xist
# (down in males) and the Y genes Ddx3y, Uty, Eif2s3y, Kdm5d (up in males) in
# any cell type with real power. The script reports them per cell type. If they
# are absent, the pseudobulk matrix or the animal metadata is wrong and nothing
# else in the output should be believed.
#
# The 4 transgene features are dropped before DE: they are the experimental
# label, not biology, and their counts are presence/absence calls.
#
# A cell type is skipped when fewer than MIN_ANIMALS animals clear MIN_NUCLEI
# nuclei, or when the surviving animals cannot support the design -- both are
# reported rather than silently producing a rank-deficient fit.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

source("/scratch/user/uqschauq/Insular/Code/lib_celltype_palette.R")

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
IN_UMAP <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_umap.rds")
OUTDIR  <- file.path(PROJECT, "Results", "9_pseudobulk_deg")

ANNOT       <- "annot"
TRANSGENES  <- gsub("_", "-", c("EGFP", "tdTomato", "EBFP2", "AAV_shared"), fixed = TRUE)
MIN_NUCLEI  <- 10     # nuclei an animal needs to contribute a pseudobulk library
MIN_ANIMALS <- 8      # animals a cell type needs before it is tested at all
FDR_CUT     <- 0.05
TOP_LABEL   <- 5      # genes labelled per volcano facet
SEX_GENES   <- c("Xist", "Ddx3y", "Uty", "Eif2s3y", "Kdm5d")
SEED        <- 42

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

if (!file.exists(IN))
  stop("no annotated object at ", IN,
       " -- run the saveRDS line at the end of 7_annotation_featureplots.R first")
if (file.exists(IN_UMAP) && file.mtime(IN) < file.mtime(IN_UMAP))
  msg("WARNING: ", basename(IN), " is OLDER than ", basename(IN_UMAP))

msg("reading ", basename(IN))
obj <- readRDS(IN)
if (!ANNOT %in% names(obj@meta.data))
  stop("no '", ANNOT, "' column in the object")
stopifnot(all(c("sample", "mouse", "diet", "sex") %in% names(obj@meta.data)))

md <- as.data.table(obj@meta.data, keep.rownames = "barcode")
lv <- levels(md[[ANNOT]])
ct_levels <- if (!is.null(lv)) rev(lv) else sort(unique(as.character(md[[ANNOT]])))
md[, annot := factor(as.character(get(ANNOT)), levels = ct_levels)]
md[, diet := factor(as.character(diet), levels = c("SD", "HFD"))]
md[, sex  := factor(as.character(sex),  levels = c("F", "M"))]
md[, sample := as.character(sample)]
msg("loaded ", format(nrow(md), big.mark = ","), " nuclei | ",
    nlevels(md$annot), " cell types | ", uniqueN(md$sample), " animals")

## ---- pseudobulk: sum counts per animal per cell type ------------------------

cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")
stopifnot(identical(colnames(cnt), md$barcode))
cnt <- cnt[setdiff(rownames(cnt), TRANSGENES), , drop = FALSE]
rm(obj); invisible(gc())

key <- factor(paste(md$annot, md$sample, sep = "__"))
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
lib[, c("annot", "sample") := tstrsplit(library, "__", fixed = TRUE)]
lib <- merge(lib, unique(md[, .(sample, mouse, diet, sex)]), by = "sample", sort = FALSE)
ncell <- md[, .N, by = .(annot, sample)]
lib <- merge(lib, ncell, by.x = c("annot", "sample"), by.y = c("annot", "sample"),
             all.x = TRUE)
setnames(lib, "N", "n_nuclei")
lib[is.na(n_nuclei), n_nuclei := 0L]
lib[, lib_size := colSums(pb)[library]]
lib[, annot := factor(annot, levels = ct_levels)]
setorder(lib, annot, diet, sex, sample)
fwrite(lib, file.path(OUTDIR, "pseudobulk_library_qc.csv"))
msg("nuclei per library: median ", median(lib$n_nuclei),
    ", min ", min(lib$n_nuclei), ", ", lib[n_nuclei < MIN_NUCLEI, .N],
    " below MIN_NUCLEI = ", MIN_NUCLEI)

## ---- DE, per cell type ------------------------------------------------------

MODELS <- list(
  A = list(formula = ~ diet,       coefs = "dietHFD"),
  B = list(formula = ~ diet * sex, coefs = c("dietHFD", "sexM", "dietHFD:sexM")))

run_one <- function(ct, model_name) {
  spec <- MODELS[[model_name]]
  l <- lib[annot == ct & n_nuclei >= MIN_NUCLEI]
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
    data.table(annot = ct, model = model_name, contrast = cf,
               gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM,
               PValue = tt$PValue, FDR = tt$FDR)
  }))
  attr(out, "meta") <- data.table(
    annot = ct, model = model_name, n_libraries = nrow(l),
    n_genes_tested = sum(keep), residual_df = nrow(l) - ncol(des),
    common_bcv = round(sqrt(y$common.dispersion), 3),
    n_SD = sum(l$diet == "SD"), n_HFD = sum(l$diet == "HFD"),
    min_nuclei = min(l$n_nuclei), median_nuclei = as.numeric(median(l$n_nuclei)))
  out
}

res <- list(); meta <- list()
for (mn in names(MODELS)) {
  msg("model ", mn, ": ", deparse(MODELS[[mn]]$formula))
  for (ct in ct_levels) {
    r <- run_one(ct, mn)
    if (is.null(r)) next
    res[[length(res) + 1]]  <- r
    meta[[length(meta) + 1]] <- attr(r, "meta")
  }
}
if (!length(res)) stop("no cell type produced a fit")
de   <- rbindlist(res)
meta <- rbindlist(meta)
de[, annot := factor(annot, levels = ct_levels)]

fwrite(de[order(model, contrast, annot, PValue)],
       file.path(OUTDIR, "pseudobulk_de_all.csv.gz"), compress = "gzip")
fwrite(meta, file.path(OUTDIR, "pseudobulk_fit_summary.csv"))

sig <- de[FDR < FDR_CUT, .(n_de = .N,
                           n_up = sum(logFC > 0), n_down = sum(logFC < 0)),
          by = .(model, contrast, annot)]
full <- CJ(model = unique(de$model), contrast = unique(de$contrast),
           annot = factor(ct_levels, levels = ct_levels), unique = TRUE)
full <- full[unique(de[, .(model, contrast)]), on = .(model, contrast)]
sig <- sig[full, on = .(model, contrast, annot)]
sig[is.na(n_de), `:=`(n_de = 0L, n_up = 0L, n_down = 0L)]
sig[, label := factor(fifelse(model == "A", "A: ~ diet",
                       fifelse(contrast == "dietHFD",     "B: diet (in females)",
                       fifelse(contrast == "sexM",        "B: sex (in SD)",
                                                          "B: diet x sex"))),
                      levels = c("A: ~ diet", "B: diet (in females)",
                                 "B: sex (in SD)", "B: diet x sex"))]
fwrite(sig[order(label, annot)], file.path(OUTDIR, "pseudobulk_de_counts.csv"))
msg("DE genes at FDR < ", FDR_CUT, ":")
print(dcast(sig, annot ~ label, value.var = "n_de"))

msg("sex-gene sanity check (model B, sexM coefficient):")
sx <- de[model == "B" & contrast == "sexM" & gene %in% SEX_GENES,
         .(annot, gene, logFC = round(logFC, 2), FDR = signif(FDR, 2))]
print(dcast(sx, annot ~ gene, value.var = "logFC"))
hit <- sx[FDR < FDR_CUT, uniqueN(annot)]
if (!nrow(sx) || hit == 0)
  msg("WARNING: no sex gene reaches FDR < ", FDR_CUT, " in any cell type -- ",
      "check the pseudobulk matrix and the animal metadata before reading anything else")

## ---- figures ----------------------------------------------------------------

CONTRAST_COLS <- c("A: ~ diet" = "#dd5129", "B: diet (in females)" = "#0f7ba2",
                   "B: sex (in SD)" = "#43b284", "B: diet x sex" = "#fab255")

pA <- ggplot(sig, aes(annot, n_de, fill = label)) +
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
  lab <- dd[sig == TRUE][order(FDR), head(.SD, TOP_LABEL), by = annot]
  p <- ggplot(dd, aes(logFC, -log10(PValue))) +
    geom_point(aes(colour = sig), size = 0.5, alpha = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#bb292c"),
                        name = NULL, labels = c("n.s.", paste0("FDR < ", FDR_CUT))) +
    facet_wrap(~ annot, ncol = 5, scales = "free_y") +
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
ggsave(file.path(OUTDIR, "pseudobulk_overview.png"), fig, width = 14, height = 12, dpi = 300)
ggsave(file.path(OUTDIR, "pseudobulk_overview.pdf"), fig, width = 14, height = 12)
msg("written pseudobulk_overview.png/.pdf")

vI <- volcano(de[model == "B" & contrast == "dietHFD:sexM"],
              "log2FC of the diet effect in males vs females   (model B: ~ diet * sex)")
ggsave(file.path(OUTDIR, "pseudobulk_volcano_interaction.png"), vI,
       width = 14, height = 8, dpi = 300)
ggsave(file.path(OUTDIR, "pseudobulk_volcano_interaction.pdf"), vI, width = 14, height = 8)

vS <- volcano(de[model == "B" & contrast == "sexM"],
              "log2FC, M vs F   (model B: ~ diet * sex, sex coefficient)")
ggsave(file.path(OUTDIR, "pseudobulk_volcano_sex.png"), vS,
       width = 14, height = 8, dpi = 300)
msg("written interaction and sex volcanoes")

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
