# Azimuth reference mapping of a neuronal subset onto the Allen mouse
# motor-cortex reference, summarised against a chosen subclustering resolution.
# The lineage is the first argument: `excitatory` (the default, res 0.2 = the
# 16-cluster solution chosen in 10_1) or `inhibitory` (res 0.2 = 16 clusters,
# the only swept resolution with no cluster below the size floor).
#
# WHAT THIS IS AND IS NOT. RunAzimuth SCTransforms the query, finds transfer
# anchors against the reference's stored SPCA reduction, transfers the reference
# labels and projects the query into the reference UMAP. It is a per-NUCLEUS
# prediction: nothing about our clustering constrains it. That is the point --
# the agreement between the transferred labels and our 16 clusters is the
# evidence, so this script never re-labels a cluster automatically. It writes
# the cross-tabulation and a majority call per cluster, and leaves the naming to
# the user, exactly as step 7 was done by hand.
#
# TWO MISMATCHES TO KEEP IN MIND WHEN READING THE OUTPUT.
#   1. The reference is primary motor cortex (MOp, 159,738 cells); this is
#      anterior insula. The neocortical excitatory hierarchy (L2/3 IT, L4/5 IT,
#      L5 IT, L5 ET, L5/6 NP, L6 CT, L6b, Car3) is shared across the neocortex,
#      so `subclass` should transfer sensibly, but insula-specific identity --
#      and anything agranular/claustral -- has no matching reference label and
#      will be forced onto the nearest MOp type. Low mapping.score is the flag.
#   2. The reference is 10x v3 WHOLE CELLS; this is single NUCLEI. Nuclear
#      transcriptomes are intron-rich and cytoplasm-poor, so prediction and
#      mapping scores run systematically lower than the reference's own. Compare
#      scores BETWEEN our clusters, not against an absolute cutoff.
#
# WHY THIS DOES NOT CALL RunAzimuth. Azimuth 0.5.1's RunAzimuth() hands the
# query straight to FindTransferAnchors with normalization.method = "SCT" and
# query.assay = DefaultAssay(query) -- without ever SCTransforming the query.
# Seurat 5.3.0 refuses that combination outright:
#   "Reference assay is SCT, but query assay is RNA. Mixing SCT and non-SCT in
#    FindTransferAnchors is not supported."
# (job 27868892, died in 43 s). Azimuth's own Shiny app never takes that path:
# server.R SCTransforms the query against the REFERENCE'S STORED SCT MODEL into
# an assay called refAssay first, and then calls FindTransferAnchors with
# recompute.residuals = FALSE -- a flag RunAzimuth does not expose.
#
# So the app's sequence is reproduced here directly (server.R:1085 and 1151,
# then the tail of RunAzimuth verbatim: TransferData -> IntegrateEmbeddings ->
# query_ref.nn -> NNTransform -> project refUMAP -> MappingScore). Same
# algorithm, same parameters, same reference -- just without the wrapper that
# assumes an older Seurat.
#
# Usage:
#   sbatch 10_2_azimuth.sh                            # excitatory, res 0.2
#   sbatch 10_2_azimuth.sh inhibitory 0.2
#   sbatch 10_2_azimuth.sh inhibitory 0.2 smoke       # 8k nuclei into .smoke/
#
# The reference carries the GABAergic subclasses (Pvalb, Sst, Sst Chodl, Vip,
# Lamp5, Sncg, Meis2) alongside the glutamatergic ones, so the same reference
# and the same call serve both lineages; only the query changes.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratData)
  library(Azimuth)
  library(ggplot2)
  library(data.table)
})

PROJECT  <- "/scratch/user/uqschauq/Insular"
LINEAGES <- c("excitatory", "inhibitory")

SEQ_LO <- "#eef4fc"
SEQ_HI <- "#0d366b"
SEED   <- 42

args    <- commandArgs(trailingOnly = TRUE)
LINEAGE <- if (length(args) >= 1 && nzchar(args[1])) tolower(args[1]) else "excitatory"
RES     <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(args[2]) else 0.2
SMOKE   <- length(args) >= 3 && identical(tolower(args[3]), "smoke")
SMOKE_N <- 8000L
if (!LINEAGE %in% LINEAGES)
  stop("first argument must be one of: ", paste(LINEAGES, collapse = ", "))
if (!is.finite(RES)) stop("second argument must be a resolution, e.g. 0.2")

IN     <- file.path(PROJECT, "Data", "Seurat",
                    sprintf("HFDseq_v2_sub_%s.rds", LINEAGE))
OUTDIR <- file.path(PROJECT, "Results", sprintf("10_2_%s_azimuth", LINEAGE))
if (SMOKE) OUTDIR <- file.path(OUTDIR, ".smoke")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
options(future.globals.maxSize = 200 * 1024^3)

refdir <- system.file("azimuth", package = "mousecortexref.SeuratData")
if (!nzchar(refdir) || !all(c("ref.Rds", "idx.annoy") %in% list.files(refdir)))
  stop("mouse cortex reference not installed -- run Code/install_azimuth.sh first")
msg("reference: ", refdir)

stopifnot(file.exists(IN))
msg("reading ", basename(IN), " (", round(file.size(IN) / 2^30, 1), " GB)")
sub <- readRDS(IN)
msg("loaded ", format(ncol(sub), big.mark = ","), " ", LINEAGE, " nuclei")

res_col <- paste0("RNA_snn_res.", RES)
if (!res_col %in% names(sub@meta.data))
  stop("no '", res_col, "' column in the object")
cl <- sub@meta.data[[res_col]]
cl <- factor(as.character(cl), levels = as.character(sort(as.integer(unique(as.character(cl))))))
sub$dot_cluster <- cl
K <- nlevels(cl)
msg(res_col, ": ", K, " clusters")

DefaultAssay(sub) <- "RNA"
if (length(Layers(sub[["RNA"]], search = "counts")) > 1) sub <- JoinLayers(sub)

if (SMOKE) {
  set.seed(SEED)
  keep <- sample(colnames(sub), min(SMOKE_N, ncol(sub)))
  sub <- subset(sub, cells = keep)
  sub$dot_cluster <- droplevels(sub$dot_cluster)
  cl <- sub$dot_cluster; K <- nlevels(cl)
  msg("SMOKE: ", ncol(sub), " nuclei, ", K, " clusters -> ", OUTDIR)
}

## ---- the mapping ------------------------------------------------------------

refs <- LoadReference(path = refdir)
ref  <- refs$map
dims <- as.double(slot(ref, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)
msg("reference: ", ncol(ref), " cells, ", dims, " dims")

# Seurat 5.3.0's FindTransferAnchors hardcodes reference[["SCT"]] when
# normalization.method = "SCT" -- it reads the model list from an assay named
# literally "SCT", ignoring the reference.assay argument it was given, and its
# own DietSeurat call has already dropped every other assay by then. Azimuth
# references name their assay "refAssay", so the lookup dies with
# "'SCT' not found in this Seurat object" (job 27868972). Renaming the assay on
# the reference side satisfies the hardcoded name; the reductions are re-pointed
# at it the same way Azimuth's reference builder points them at refAssay.
ref[["SCT"]] <- ref[["refAssay"]]
DefaultAssay(ref) <- "SCT"
for (rd in Reductions(ref)) DefaultAssay(ref[[rd]]) <- "SCT"
ref[["refAssay"]] <- NULL
msg("reference assay renamed refAssay -> SCT for FindTransferAnchors")

# Gene names. Both sides are mouse Gene.name, so this is expected to be a no-op;
# it needs the network for the homolog table, so a failure is logged and the
# original names kept rather than killing the run.
sub <- tryCatch(
  ConvertGeneNames(sub, reference.names = rownames(ref),
                   homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds"),
  error = function(e) { msg("ConvertGeneNames failed (", conditionMessage(e),
                            ") -- keeping query gene names"); sub })
ov <- length(intersect(rownames(sub), rownames(ref)))
msg("query shares ", ov, " / ", nrow(ref), " reference features")
if (ov < 0.5 * nrow(ref))
  stop("query and reference share too few features (", ov, ") -- gene names do not match")

# SCTransform against the reference's own SCT model, into refAssay.
model <- slot(ref[["SCT"]], "SCTModel.list")[["refmodel"]]
sct <- function(method) SCTransform(sub, residual.features = rownames(ref),
                                    reference.SCT.model = model, method = method,
                                    do.correct.umi = FALSE, do.scale = FALSE,
                                    do.center = TRUE, new.assay.name = "refAssay",
                                    verbose = TRUE)
msg("SCTransform against the reference model (glmGamPoi)")
sub <- tryCatch(suppressWarnings(sct("glmGamPoi")),
                error = function(e) { msg("glmGamPoi failed (", conditionMessage(e),
                                          ") -- falling back to poisson")
                                      suppressWarnings(sct("poisson")) })
DefaultAssay(sub) <- "refAssay"
msg("query normalised: ", nrow(sub[["refAssay"]]), " residual features")

msg("FindTransferAnchors")
anchors <- FindTransferAnchors(
  reference = ref, query = sub, k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "SCT", query.assay = "refAssay",
  reference.reduction = "refDR", normalization.method = "SCT",
  recompute.residuals = FALSE,
  features = rownames(Loadings(ref[["refDR"]])),
  dims = 1:dims, n.trees = 20, mapping.score.k = 100, verbose = TRUE)
msg("anchors: ", nrow(slot(anchors, "anchors")))

annotation.levels <- grep("^nCount|^nFeature|^ori", names(slot(ref, "meta.data")),
                          invert = TRUE, value = TRUE)
msg("transferring: ", paste(annotation.levels, collapse = ", "))
refdata <- lapply(annotation.levels, function(x) ref[[x, drop = TRUE]])
names(refdata) <- annotation.levels

sub <- TransferData(reference = ref, query = sub, query.assay = "refAssay",
                    dims = 1:dims, anchorset = anchors, refdata = refdata,
                    n.trees = 20, store.weights = TRUE, k.weight = 50, verbose = TRUE)
sub <- IntegrateEmbeddings(anchorset = anchors, reference = ref, query = sub,
                           query.assay = "refAssay", reductions = "pcaproject",
                           reuse.weights.matrix = TRUE, verbose = TRUE)
sub[["query_ref.nn"]] <- FindNeighbors(object = Embeddings(ref[["refDR"]]),
                                       query = Embeddings(sub[["integrated_dr"]]),
                                       return.neighbor = TRUE, l2.norm = TRUE)
# The shipped reference is downsampled relative to the one the UMAP model was
# fit on; NNTransform corrects the neighbour indices for that.
sub <- Azimuth:::NNTransform(object = sub, meta.data = ref[[]])
sub[["ref.umap"]] <- RunUMAP(object = sub[["query_ref.nn"]],
                             reduction.model = ref[["refUMAP"]],
                             reduction.key = "UMAP_", verbose = FALSE)
sub <- AddMetaData(sub, metadata = MappingScore(anchors = anchors, ndim = dims),
                   col.name = "mapping.score")
msg("mapping done")

pred_cols <- grep("^predicted\\.", names(sub@meta.data), value = TRUE)
levels_avail <- sub("^predicted\\.", "", grep("\\.score$", pred_cols, invert = TRUE,
                                              value = TRUE))
msg("annotation levels transferred: ", paste(levels_avail, collapse = ", "))
if (!length(levels_avail)) stop("RunAzimuth returned no predicted.* columns")

md <- as.data.table(sub@meta.data, keep.rownames = "barcode")
keep <- c("barcode", "dot_cluster", "mouse", "sample", "diet", "sex",
          "parent_cluster", pred_cols, "mapping.score")
pred <- md[, intersect(keep, names(md)), with = FALSE]
fwrite(pred, file.path(OUTDIR, sprintf("azimuth_predictions_res%s.csv.gz", RES)))
saveRDS(pred, file.path(OUTDIR, sprintf("azimuth_predictions_res%s.rds", RES)))
msg("wrote per-nucleus predictions (", nrow(pred), " rows)")

if ("ref.umap" %in% Reductions(sub))
  fwrite(as.data.table(Embeddings(sub, "ref.umap"), keep.rownames = "barcode"),
         file.path(OUTDIR, "ref_umap_embeddings.csv.gz"))

## ---- cluster x transferred label, per annotation level ----------------------

for (lvl in levels_avail) {
  pcol <- paste0("predicted.", lvl)
  scol <- paste0("predicted.", lvl, ".score")
  dt <- pred[, .(cluster = dot_cluster, label = get(pcol),
                 score = if (scol %in% names(pred)) get(scol) else NA_real_,
                 mscore = mapping.score)]

  tab <- dcast(dt[, .N, by = .(cluster, label)], cluster ~ label,
               value.var = "N", fill = 0)
  fwrite(tab, file.path(OUTDIR, sprintf("crosstab_%s_counts_res%s.csv", lvl, RES)))

  pct <- dt[, .N, by = .(cluster, label)][, pct := 100 * N / sum(N), by = cluster]
  fwrite(dcast(pct, cluster ~ label, value.var = "pct", fill = 0),
         file.path(OUTDIR, sprintf("crosstab_%s_percent_res%s.csv", lvl, RES)))

  # Majority call per cluster -- reported, never written back onto the object.
  maj <- pct[order(cluster, -pct)][, .(
    top_label   = label[1], top_pct = round(pct[1], 1),
    second_label = if (.N > 1) label[2] else NA_character_,
    second_pct   = if (.N > 1) round(pct[2], 1) else NA_real_,
    n_labels_over_10pct = sum(pct >= 10)), by = cluster]
  sc <- dt[, .(n = .N, median_pred_score = round(median(score, na.rm = TRUE), 3),
               median_mapping_score = round(median(mscore, na.rm = TRUE), 3)),
           by = cluster]
  maj <- merge(maj, sc, by = "cluster", sort = FALSE)[order(cluster)]
  fwrite(maj, file.path(OUTDIR, sprintf("majority_call_%s_res%s.csv", lvl, RES)))
  msg("--- ", lvl, " ---"); print(maj)

  # Heatmap: % of each cluster assigned to each reference label.
  ord <- pct[, .(tot = sum(pct)), by = label][order(-tot), as.character(label)]
  pct[, label := factor(as.character(label), levels = ord)]
  h <- ggplot(pct, aes(label, cluster, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = fifelse(pct >= 5, sprintf("%.0f", pct), ""),
                  colour = pct > 55), size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#3d3d3a"), guide = "none") +
    scale_fill_gradient(low = SEQ_LO, high = SEQ_HI, limits = c(0, 100),
                        name = "% of\ncluster") +
    scale_y_discrete(limits = rev(levels(pct$cluster))) +
    labs(x = paste("Azimuth predicted", lvl), y = "subcluster",
         title = sprintf("%s subclusters (res %s) vs mouse MOp reference",
                         tools::toTitleCase(LINEAGE), RES),
         subtitle = "row-normalised; whole-cell MOp reference mapped onto insula nuclei") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(OUTDIR, sprintf("heatmap_%s_res%s.png", lvl, RES)), h,
         width = 0.28 * uniqueN(pct$label) + 5, height = 0.30 * K + 3,
         dpi = 200, limitsize = FALSE)
}

## ---- UMAPs ------------------------------------------------------------------

lvl_main <- if ("subclass" %in% levels_avail) "subclass" else levels_avail[1]
for (rd in intersect(c("umap", "ref.umap"), Reductions(sub))) {
  ggsave(file.path(OUTDIR, sprintf("%s_predicted_%s.png", gsub("\\.", "_", rd), lvl_main)),
         DimPlot(sub, reduction = rd, group.by = paste0("predicted.", lvl_main),
                 label = TRUE, label.size = 3, repel = TRUE, raster = TRUE,
                 raster.dpi = c(1400, 1400), pt.size = 0.4) + coord_equal(),
         width = 11, height = 8.5, dpi = 200)
  ggsave(file.path(OUTDIR, sprintf("%s_subcluster.png", gsub("\\.", "_", rd))),
         DimPlot(sub, reduction = rd, group.by = "dot_cluster", label = TRUE,
                 label.size = 4, raster = TRUE, raster.dpi = c(1400, 1400),
                 pt.size = 0.4) + coord_equal() + NoLegend(),
         width = 9, height = 8.5, dpi = 200)
}

score_feats <- intersect(c("mapping.score", paste0("predicted.", lvl_main, ".score")),
                         names(sub@meta.data))
ggsave(file.path(OUTDIR, "umap_scores.png"),
       FeaturePlot(sub, reduction = "umap", features = score_feats, raster = TRUE,
                   raster.dpi = c(1200, 1200), pt.size = 0.4,
                   ncol = length(score_feats)) &
         scale_colour_gradient(low = SEQ_LO, high = SEQ_HI),
       width = 6 * length(score_feats), height = 5.5, dpi = 150)

vp <- ggplot(pred, aes(dot_cluster, mapping.score)) +
  geom_violin(fill = SEQ_LO, colour = "#3d3d3a", linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.3) +
  labs(x = "subcluster", y = "Azimuth mapping score",
       subtitle = "low = the reference has no good match for these nuclei") +
  theme_bw(base_size = 10)
ggsave(file.path(OUTDIR, sprintf("mapping_score_by_cluster_res%s.png", RES)), vp,
       width = 0.45 * K + 3, height = 5, dpi = 200)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)
