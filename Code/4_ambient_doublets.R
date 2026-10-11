#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 4_ambient_doublets.R  --  phase 2 of step 4, run by 4_ambient_doublets.sh.
#
# Concatenates the 12 CellBender-corrected matrices into one object and calls
# doublets on it. Three things here are deliberate and easy to get wrong:
#
#  1. ORDER. The 12 matrices are concatenated in aggregation.csv row order and
#     the barcode suffix -1..-12 is assigned from that order, so the object is
#     byte-compatible with the metadata convention every later script uses.
#     Never glob the directory instead; the glob is alphabetical and would put
#     Ech10 third.
#
#  2. TRANSGENES ARE NOT DECONTAMINATED. CellBender corrects all 33,700
#     features, but the projection call is presence/absence at >= 1 UMI, and a
#     generic decontaminator will zero an unknown mix of true and false single
#     UMIs with no way to tell them apart afterwards. The 4 transgene rows are
#     therefore overwritten with their RAW counts for the retained barcodes.
#     The ambient-transgene question is handled separately, by the AAV_shared
#     denominator, which is what it was kept in the reference for.
#
#  3. scDblFinder SEES MOUSE GENES ONLY. AAV_shared can carry a large UMI count
#     in a transgene-positive nucleus, so leaving the transgenes in would let a
#     4-feature block drive the HVG/PCA step the doublet score is built on.
#
# CellBender re-does cell calling, so the nucleus set here is CellBender's, not
# cellranger's 323,898. That is expected, not a bug.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(BiocParallel)
  library(data.table)
  library(ggplot2)
})

source("/scratch/user/uqschauq/Insular/Code/lib_read_cellbender.R")

PROJECT <- "/scratch/user/uqschauq/Insular"
AGGCSV  <- file.path(PROJECT, "Data", "Aggr_v2", "HFDseq_aggr_v2", "outs", "aggregation.csv")
OUTDIR  <- file.path(PROJECT, "Results", "4_ambient_doublets")
CBDIR   <- file.path(OUTDIR, "cellbender")
OBJDIR  <- file.path(PROJECT, "Data", "Seurat")
OBJOUT  <- file.path(OBJDIR, "HFDseq_v2_ambient_doublets.rds")

TRANSGENES <- c("EGFP", "tdTomato", "EBFP2", "AAV_shared")
WORKERS <- 6
SEED    <- 42
SEQ_HI  <- "#0d366b"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---- read the 12 corrected matrices ---------------------------------------

agg <- fread(AGGCSV)
setnames(agg, gsub("\r", "", names(agg), fixed = TRUE))
agg[, gem_group := .I]
norm_fluor <- function(x) {
  x <- trimws(gsub("\r", "", x, fixed = TRUE))
  ifelse(x == "EBFP", "EBFP2", x)
}
for (r in c("BLA", "NAc", "LHr")) if (r %in% names(agg)) set(agg, j = r, value = norm_fluor(agg[[r]]))
msg("aggregation.csv: ", nrow(agg), " samples in order ", paste(agg$sample_id, collapse = ", "))

mats <- vector("list", nrow(agg))
ambient <- vector("list", nrow(agg))

for (i in seq_len(nrow(agg))) {
  s   <- agg$sample_id[i]
  cbh5 <- file.path(CBDIR, paste0(s, "_filtered.h5"))
  rawh5 <- file.path(PROJECT, "Data", "Counts_v2", s, "outs", "raw_feature_bc_matrix.h5")
  if (!file.exists(cbh5))  stop("missing CellBender output: ", cbh5)
  if (!file.exists(rawh5)) stop("missing raw matrix: ", rawh5)

  cb <- read_cellbender_h5(cbh5)

  raw <- Read10X_h5(rawh5)
  if (is.list(raw)) raw <- raw[["Gene Expression"]]

  # CellBender keeps cellranger's feature order, but assert rather than assume:
  # a mismatch here would silently transpose the transgene rows.
  if (!identical(rownames(cb), rownames(raw)))
    stop("feature names differ between CellBender and raw for ", s)
  missing_tg <- setdiff(TRANSGENES, rownames(cb))
  if (length(missing_tg))
    stop("transgenes absent for ", s, ": ", paste(missing_tg, collapse = ", "))

  keep <- colnames(cb)
  if (!all(keep %in% colnames(raw)))
    stop("CellBender kept barcodes absent from the raw matrix for ", s)
  raw_keep <- raw[, keep, drop = FALSE]

  # Ambient accounting BEFORE the transgene rows are put back, so the number
  # describes what CellBender actually did to the mouse transcriptome.
  mouse_rows <- setdiff(rownames(cb), TRANSGENES)
  raw_umi <- sum(Matrix::colSums(raw_keep[mouse_rows, , drop = FALSE]))
  cb_umi  <- sum(Matrix::colSums(cb[mouse_rows, , drop = FALSE]))
  ambient[[i]] <- data.table(
    sample = s, mouse = agg$mouse[i], diet = agg$diet[i], sex = agg$sex[i],
    n_cells_cellbender = length(keep),
    n_cells_cellranger = NA_integer_,
    umi_raw = raw_umi, umi_corrected = cb_umi,
    pct_removed = round(100 * (raw_umi - cb_umi) / raw_umi, 2))

  # (2) transgene rows back to raw.
  cb[TRANSGENES, ] <- raw_keep[TRANSGENES, ]

  colnames(cb) <- paste0(sub("-[0-9]+$", "", keep), "-", agg$gem_group[i])
  mats[[i]] <- cb
  rm(raw, raw_keep, cb); invisible(gc())
  msg("  ", s, " (gem ", agg$gem_group[i], "): ", length(keep), " nuclei, ",
      ambient[[i]]$pct_removed, "% of mouse UMIs removed as ambient")
}

ambient <- rbindlist(ambient)
counts <- do.call(cbind, mats)
rm(mats); invisible(gc())
if (anyDuplicated(colnames(counts))) stop("duplicate barcodes after concatenation")
msg("combined: ", nrow(counts), " features x ", ncol(counts), " nuclei, ",
    format(length(counts@x), big.mark = ","), " non-zeros")

ambient[, pct_of_total := round(100 * n_cells_cellbender / sum(n_cells_cellbender), 2)]
fwrite(ambient, file.path(OUTDIR, "cellbender_ambient_summary.csv"))
msg("ambient removed per animal: ", round(min(ambient$pct_removed), 2), "% - ",
    round(max(ambient$pct_removed), 2), "% (pooled ",
    round(100 * (sum(ambient$umi_raw) - sum(ambient$umi_corrected)) / sum(ambient$umi_raw), 2), "%)")
print(ambient[, .(sample, n_cells_cellbender, pct_removed)])

## ---- metadata from the barcode suffix -------------------------------------

bc  <- colnames(counts)
gem <- as.integer(sub("^.*-", "", bc))
if (anyNA(gem)) stop("barcodes without a -N GEM suffix")
if (!setequal(unique(gem), agg$gem_group)) stop("GEM groups do not match aggregation.csv")

md <- data.frame(
  gem_group = gem,
  sample    = factor(agg$sample_id[gem], levels = agg$sample_id),
  mouse     = agg$mouse[gem],
  diet      = factor(agg$diet[gem], levels = c("SD", "HFD")),
  sex       = factor(agg$sex[gem],  levels = c("F", "M")),
  BLA = agg$BLA[gem], NAc = agg$NAc[gem], LHr = agg$LHr[gem],
  row.names = bc, stringsAsFactors = FALSE)
msg("nuclei per sample:"); print(table(md$sample))

## ---- doublets, on mouse genes only ----------------------------------------

mouse_genes <- setdiff(rownames(counts), TRANSGENES)
msg("scDblFinder input: ", length(mouse_genes), " mouse genes (", length(TRANSGENES),
    " transgenes withheld)")

sce <- SingleCellExperiment(list(counts = counts[mouse_genes, ]))
sce$mouse <- md$mouse

ncell <- table(sce$mouse)
msg("per-animal nuclei and implied default dbr:")
for (m in names(ncell))
  cat(sprintf("    %-8s %6d   dbr ~ %.3f\n", m, ncell[[m]], ncell[[m]] / 1000 * 0.01))

msg("running scDblFinder (per animal, random artificial doublets, ", WORKERS, " workers)")
sce <- scDblFinder(sce, samples = "mouse",
                   BPPARAM = MulticoreParam(workers = WORKERS, RNGseed = SEED))
stopifnot(identical(colnames(sce), bc))
msg("done: ", sum(sce$scDblFinder.class == "doublet"), " doublets of ", ncol(sce),
    " (", round(100 * mean(sce$scDblFinder.class == "doublet"), 2), "%)")

md$dbl_score <- as.numeric(sce$scDblFinder.score)
md$dbl_class <- as.character(sce$scDblFinder.class)
rm(sce); invisible(gc())

## ---- object ---------------------------------------------------------------

obj <- CreateSeuratObject(counts = counts, project = "HFDseq_v2",
                          meta.data = md, min.cells = 0, min.features = 0)
renamed <- data.table(matrix_name = rownames(counts), object_name = rownames(obj))
renamed <- renamed[matrix_name != object_name]
if (nrow(renamed)) {
  msg("feature names rewritten by Seurat: ",
      paste(renamed[, paste0(matrix_name, " -> ", object_name)], collapse = ", "))
  fwrite(renamed, file.path(OUTDIR, "feature_name_changes.csv"))
}
rm(counts); invisible(gc())

dt <- as.data.table(obj@meta.data, keep.rownames = "barcode")
fwrite(dt[, .(barcode, sample, mouse, diet, sex, dbl_score, dbl_class)],
       file.path(OUTDIR, "doublet_calls.csv.gz"), compress = "gzip")

by_mouse <- dt[, .(n = .N, n_doublet = sum(dbl_class == "doublet"),
                   median_score = round(median(dbl_score), 4)), by = .(sample, mouse, diet, sex)]
by_mouse[, pct_doublet := round(100 * n_doublet / n, 2)]
setorder(by_mouse, -pct_doublet)
fwrite(by_mouse, file.path(OUTDIR, "doublet_by_mouse.csv"))
msg("doublets by animal:"); print(by_mouse)

ggsave(file.path(OUTDIR, "doublet_score_hist.png"),
       ggplot(dt, aes(dbl_score)) +
         geom_histogram(bins = 80, fill = SEQ_HI) +
         facet_wrap(~ sample, ncol = 4, scales = "free_y") +
         labs(x = "scDblFinder score", y = "nuclei",
              title = "Doublet score distribution per animal",
              subtitle = "after CellBender ambient removal; mouse genes only") +
         theme_bw(base_size = 9) + theme(strip.background = element_blank()),
       width = 12, height = 8, dpi = 150)

ggsave(file.path(OUTDIR, "cellbender_pct_removed.png"),
       ggplot(ambient, aes(reorder(sample, pct_removed), pct_removed, fill = diet)) +
         geom_col(width = 0.75) +
         geom_text(aes(label = sprintf("%.1f%%", pct_removed)), hjust = -0.15, size = 3) +
         coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
         labs(x = NULL, y = "% of mouse-gene UMIs removed as ambient",
              title = "CellBender ambient removal per animal",
              subtitle = paste0("fpr = 0; transgene features excluded from this total ",
                                "and carried through uncorrected")) +
         theme_classic(base_size = 11),
       width = 9, height = 5, dpi = 150)

msg("saving ", OBJOUT)
saveRDS(obj, OBJOUT, compress = FALSE)
msg("written: ", round(file.size(OBJOUT) / 2^30, 2), " GB, ", ncol(obj), " nuclei, ",
    nrow(obj), " features -- NOTHING FILTERED, doublets are flagged only")

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
