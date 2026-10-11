# Two-layer cluster markers, lifted from 6_seurat_umap_markers.R so that
# 10_neuronal_subclustering.R runs the SAME test on subclusters as the pipeline
# ran on the parent clusters. Step 6 still carries its own inline copy and was
# deliberately not rewired -- it is a finished step and its outputs are on disk.
# If the statistics here change, change them there too.
#
#   layer 1  Seurat one-vs-rest Wilcoxon (presto-accelerated). avg_log2FC and
#            the pct.1/pct.2 detection rates annotation keys on.
#   layer 2  Pseudobulk edgeR, cluster vs rest, PAIRED WITHIN MOUSE. Wilcoxon
#            treats every nucleus as a replicate, so at n ~ 20k its p-values are
#            uninterpretable and a "marker" can be one or two animals. Layer 2
#            re-tests each candidate with the mouse as the unit and reports how
#            many animals support it.
#
# `concordant` is the column to rank on: significant in the pseudobulk test AND
# up in at least 75% of the animals tested.

cluster_markers <- function(obj, clusters, mouse,
                            prefilt_pct = 0.10, min_pct = 0.25, logfc_thr = 0.25,
                            min_cells = 10, min_mice = 4, log = message) {
  stopifnot(requireNamespace("edgeR", quietly = TRUE))
  clusters <- droplevels(as.factor(clusters))
  mouse    <- droplevels(as.factor(mouse))
  counts   <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  stopifnot(ncol(counts) == length(clusters), length(clusters) == length(mouse))

  # Detection rate per gene per cluster, via a sparse indicator multiply.
  # Binarise by overwriting the nonzero slot rather than `counts > 0`, which
  # would materialise a second matrix of the same nnz.
  bin <- counts
  bin@x <- rep.int(1, length(bin@x))
  ind_c <- Matrix::sparseMatrix(i = seq_len(ncol(counts)), j = as.integer(clusters),
                                x = 1, dims = c(ncol(counts), nlevels(clusters)),
                                dimnames = list(NULL, levels(clusters)))
  det <- (bin %*% ind_c) %*% Matrix::Diagonal(x = 1 / as.numeric(table(clusters)))
  keep_genes <- rownames(counts)[apply(as.matrix(det), 1, max) >= prefilt_pct]
  rm(bin, det); invisible(gc())
  log("gene pre-filter: ", length(keep_genes), " / ", nrow(counts),
      " genes detected in >=", prefilt_pct * 100, "% of at least one cluster")

  Seurat::Idents(obj) <- clusters
  log("layer 1: FindAllMarkers (wilcox, presto = ",
      requireNamespace("presto", quietly = TRUE), ")")
  mk <- Seurat::FindAllMarkers(obj, features = keep_genes, only.pos = TRUE,
                               min.pct = min_pct, logfc.threshold = logfc_thr,
                               test.use = "wilcox", verbose = FALSE)
  data.table::setDT(mk)
  data.table::setnames(mk, "cluster", "cluster_id")
  mk[, cluster_id := as.character(cluster_id)]
  log("layer 1: ", nrow(mk), " marker rows across ", data.table::uniqueN(mk$cluster_id),
      " clusters")

  log("layer 2: aggregating pseudobulk by (cluster x mouse)")
  grp   <- interaction(clusters, mouse, drop = FALSE, sep = "|")
  ind_g <- Matrix::sparseMatrix(i = seq_len(ncol(counts)), j = as.integer(grp),
                                x = 1, dims = c(ncol(counts), nlevels(grp)),
                                dimnames = list(NULL, levels(grp)))
  pb    <- as.matrix(counts[keep_genes, ] %*% ind_g)
  ncell <- as.integer(table(grp)); names(ncell) <- levels(grp)

  key   <- data.table::tstrsplit(colnames(pb), "|", fixed = TRUE)
  pb_cl <- factor(key[[1]], levels = levels(clusters))
  pb_ms <- factor(key[[2]], levels = levels(mouse))
  tot_ms <- sapply(levels(mouse), function(m) rowSums(pb[, pb_ms == m, drop = FALSE]))
  log("layer 2: ", ncol(pb), " units, median ", median(ncell), " cells, ",
      sum(ncell < min_cells), " below ", min_cells)

  pb_one <- function(k) {
    sel <- pb_cl == k & ncell >= min_cells
    mice_k <- as.character(pb_ms[sel])
    if (length(mice_k) < min_mice) return(NULL)
    inmat <- pb[, sel, drop = FALSE]
    restm <- tot_ms[, mice_k, drop = FALSE] - inmat
    colnames(inmat) <- paste0("in_", mice_k); colnames(restm) <- paste0("rest_", mice_k)
    mat  <- cbind(inmat, restm)
    grp2 <- factor(rep(c("in", "rest"), each = length(mice_k)), levels = c("rest", "in"))
    ms2  <- factor(rep(mice_k, 2))

    y <- edgeR::DGEList(counts = mat, group = grp2)
    design <- stats::model.matrix(~ ms2 + grp2)
    y <- y[edgeR::filterByExpr(y, design = design), , keep.lib.sizes = FALSE]
    if (nrow(y) < 10) return(NULL)
    y <- edgeR::normLibSizes(y)
    y <- edgeR::estimateDisp(y, design, robust = TRUE)
    fit <- edgeR::glmQLFit(y, design, robust = TRUE)
    res <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = ncol(design)),
                          n = Inf, sort.by = "none")$table
    cpm_all <- edgeR::cpm(y, log = FALSE)
    sup <- rowSums(cpm_all[, grp2 == "in", drop = FALSE] >
                   cpm_all[, grp2 == "rest", drop = FALSE])
    data.table::data.table(cluster_id = as.character(k), gene = rownames(res),
                           pb_log2FC = res$logFC, pb_PValue = res$PValue,
                           pb_FDR = res$FDR, n_mice_supporting = as.integer(sup),
                           n_mice_tested = length(mice_k))
  }

  pbres <- data.table::rbindlist(lapply(levels(clusters), function(k) {
    r <- tryCatch(pb_one(k), error = function(e) {
      log("  cluster ", k, ": pseudobulk failed -- ", conditionMessage(e)); NULL })
    if (!is.null(r)) log("  cluster ", k, ": ", r$n_mice_tested[1], " mice, ",
                         nrow(r), " genes tested")
    r
  }), fill = TRUE)

  if (nrow(pbres) == 0L) {
    pbres <- data.table::data.table(
      cluster_id = character(), gene = character(), pb_log2FC = numeric(),
      pb_PValue = numeric(), pb_FDR = numeric(),
      n_mice_supporting = integer(), n_mice_tested = integer())
    log("layer 2: WARNING -- no cluster produced pseudobulk stats")
  }
  skipped <- setdiff(levels(clusters), unique(pbres$cluster_id))
  if (length(skipped))
    log("layer 2: no pseudobulk stats for cluster(s) ", paste(skipped, collapse = ", "),
        " (< ", min_mice, " usable animals) -- rely on layer 1 there")

  out <- merge(mk, pbres, by = c("cluster_id", "gene"), all.x = TRUE, sort = FALSE)
  out[, concordant := !is.na(pb_FDR) & pb_FDR < 0.05 & pb_log2FC > 0 &
        n_mice_supporting >= pmax(min_mice, ceiling(0.75 * n_mice_tested))]
  out[, cluster_id := factor(cluster_id, levels = levels(clusters))]
  data.table::setorder(out, cluster_id, -avg_log2FC)
  data.table::setcolorder(out, c("cluster_id", "gene", "avg_log2FC", "pct.1", "pct.2",
                                 "p_val", "p_val_adj", "pb_log2FC", "pb_PValue",
                                 "pb_FDR", "n_mice_supporting", "n_mice_tested",
                                 "concordant"))
  out
}
