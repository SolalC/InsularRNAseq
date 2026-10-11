#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# projection_figures.R  --  the projection (retro-labelling) arm on the final
# object: one UMAP coloured by projection target, and one table of projection-
# positive nuclei per animal per cell type. Run by projection_figures.sh.
#
# HOUSEKEEPING, NOT A PIPELINE STEP, so it carries no number: it consumes the
# finished annotated object and the step-2 read tables, and changes nothing
# under Data/.
#
# ALL THREE MOLECULE TRACKS ARE FOLDED IN (asked for, 2026-09-08):
#   base   the TG assay of the object -- Cell Ranger's transgene counts
#   +resc  BLAST-rescued molecules from the unmapped pool  (steps 2.7 / 2.8)
#   +anti  antisense molecules already aligned to the cassettes  (step 2.10)
# A nucleus is projection-positive on the SUM of the three. The origin is not
# drawn on the UMAP -- a nucleus is shown as its target, however its molecules
# were found -- but the table breaks every count down by the track that first
# made that nucleus ORF-positive, so the added tracks can always be subtracted
# back out.
#
# Both added tracks bypass a filter the mouse genes passed: BLAST-rescued reads
# failed STAR's length/score cutoff, and antisense reads are not counted by
# design. Step 2.11's header is explicit that they are valid for PRESENCE /
# ABSENCE projection labelling only, never for an expression-level comparison.
# Nothing here uses them as a level.
#
# WHY THE CALLS ARE RECOMPUTED RATHER THAN READ FROM 2_1 / 2_11
# Both ran on Cell Ranger's filtered matrices -- before CellBender, before
# scDblFinder, before QC. The UMAP holds 247k nuclei, a subset of those
# barcodes, and their per-animal counts therefore do not describe it. The same
# classification is applied here to the nuclei that are actually plotted.
#
# THE CALL, unchanged from 2_1 and 2_11 and from the handover:
#   presence/absence at >= 1 UMI, NEVER an expression level. Exactly one ORF
#   over threshold -> that ORF's region for that animal; two or three ->
#   multi_ORF (ambiguous, EGFP and EBFP2 share a 98 bp identical stretch);
#   none but AAV-shared positive -> shared_only, i.e. transgene-positive but
#   projection-unassigned; nothing -> unlabelled.
#
# The fluorophore -> region map is randomised per animal and is carried on the
# object itself as the BLA / NAc / LHr metadata columns (put there by step 4
# from Data/sample_metadata.csv), so no CSV is re-read and the map cannot drift
# out of step with the object.
#
# NO PLOT CARRIES A TITLE OR SUBTITLE, matching figures_qc_to_annotation.R.
#
# TRAPS BAKED IN:
#   - transgene features are `AAV-shared` inside a Seurat object (dash, not
#     underscore) -- CreateSeuratObject rewrites `_` to `-`. The read tables use
#     `AAV_shared` and `<ORF>_cassette`, so both are normalised on the way in.
#   - TG counts are the RAW ones, never decontaminated and never normalised
#     (handover section 5), so layer = "counts" is the only correct layer.
#   - a molecule is a unique UB within a (barcode, feature), not a read. The
#     read tables are read-level and are collapsed here, exactly as 2_11 does.
#   - the rescue table is filtered on `in_called_cell`, 2_11's own filter, so
#     the numbers stay comparable with the track that was validated. Barcodes
#     CellBender called but Cell Ranger did not are therefore not rescued; the
#     count that costs is logged rather than left invisible.
#   - `annot` is revised by 10_4 from the CSV in Results/reporting, exactly as
#     figures_qc_to_annotation.R does it, or the cell types in this table would
#     disagree with fig_21 / fig_23.
#   - ~4k of 247k nuclei carry a call. Drawn in one pass they vanish under the
#     grey; the labelled nuclei are drawn as a second, larger layer on top.
#
# Outputs, in Results/reporting -- ONE figure and TWO csvs:
#   fig_27_umap_projection.png             UMAP coloured by projection target
#   projection_by_celltype_per_animal.csv  12 animals x 3 targets x cell types,
#                                          complete rectangle, zeros included,
#                                          each count split by molecule origin
#   projection_calls_per_nucleus.csv       one row per projection-positive
#                                          nucleus: target, fluorophore, ORF UMIs
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(data.table)
  library(Matrix)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
OBJ_IN  <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_annotated.rds")
OUTDIR  <- file.path(PROJECT, "Results", "reporting")
REVISED <- file.path(OUTDIR, "reclass_annot_revised.csv.gz")
RESCUE  <- file.path(PROJECT, "Results", "2_8_rescue_assign")
ANTI    <- file.path(PROJECT, "Results", "2_10_antisense_harvest")

ORFS    <- c("EGFP", "tdTomato", "EBFP2")
SHARED  <- "AAV-shared"
FEATS   <- c(ORFS, SHARED)
REGIONS <- c("BLA", "NAc", "LHr")
TRACKS  <- c("base", "blast", "antisense")
# threshold is set by the wrapper (MIN_UMI env var, default 1). Anything but 1
# tags every output `_ge<N>` so the >= 1 UMI figure and table are never overwritten.
MIN_UMI <- as.integer(Sys.getenv("MIN_UMI", "1"))
if (is.na(MIN_UMI) || MIN_UMI < 1L) stop("MIN_UMI must be an integer >= 1")
TAG <- if (MIN_UMI == 1L) "" else paste0("_ge", MIN_UMI)
SEED    <- 42

REGION_COLS <- c(BLA = "#dd5129", NAc = "#0f7ba2", LHr = "#43b284")
MOUSE_LEVELS <- paste0("Ech", 1:12)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

source(file.path(PROJECT, "Code", "lib_celltype_palette.R"))

nt <- function(p) p + labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank())

save_fig <- function(p, name, w, h, dpi = 300) {
  ggsave(file.path(OUTDIR, paste0(name, ".png")), nt(p), width = w, height = h,
         dpi = dpi, limitsize = FALSE)
  msg("  ", name, ".png")
}

# `EGFP_cassette` / `AAV_shared` in the read tables -> the object's feature names
norm_feat <- function(x) {
  x <- sub("_cassette$", "", x)
  ifelse(x == "AAV_shared", SHARED, x)
}

## ---- object ---------------------------------------------------------------

msg("reading ", basename(OBJ_IN), " (", round(file.size(OBJ_IN) / 2^30, 1), " GB)")
obj <- readRDS(OBJ_IN)
msg("loaded ", format(ncol(obj), big.mark = ","), " nuclei; threshold >= ", MIN_UMI, " UMI")

for (col in c("annot", "mouse", "sample", "gem_group", "diet", "sex", REGIONS))
  if (!col %in% names(obj@meta.data)) stop("object has no '", col, "' column")
if (!"TG" %in% Assays(obj)) stop("object has no TG assay")

# 10_4 revises `annot` where Azimuth disagrees with the manual call. Applied the
# same way figures_qc_to_annotation.R applies it, so the cell types in this
# table are the ones fig_21 and fig_23 show.
if (file.exists(REVISED)) {
  rv <- fread(REVISED)[changed == TRUE]
  md <- data.table(barcode = colnames(obj), annot = as.character(obj$annot))
  md[rv, on = "barcode", annot := i.annot_azimuth]
  n_hit <- sum(md$barcode %in% rv$barcode)
  if (n_hit != nrow(rv))
    stop(nrow(rv) - n_hit, " revised barcodes are absent from the object")
  obj$annot <- md$annot
  msg("applied Azimuth reclassification to ", format(n_hit, big.mark = ","), " nuclei")
} else {
  msg("no reclassification file -- cell types are the manual annot unchanged")
}

## ---- track 1: base, the TG assay ------------------------------------------

if (length(Layers(obj[["TG"]], search = "counts")) > 1) obj <- JoinLayers(obj, assay = "TG")
tgm <- GetAssayData(obj, assay = "TG", layer = "counts")
missing <- setdiff(FEATS, rownames(tgm))
if (length(missing)) stop("TG assay is missing: ", paste(missing, collapse = ", "))

cnt <- data.table(barcode = colnames(obj))
for (f in FEATS) set(cnt, j = paste0(f, "__base"), value = as.numeric(tgm[f, ]))

## ---- tracks 2 and 3: rescued and antisense molecules ----------------------

# Object barcodes carry the GEM-group suffix the aggregation gave them; the read
# tables carry Cell Ranger's per-sample `-1`. Translated per sample rather than
# assumed, so a reordered aggregation cannot silently mis-join.
gem <- unique(as.data.table(obj@meta.data)[, .(sample = as.character(sample),
                                               gem_group = as.integer(gem_group))])
if (nrow(gem) != uniqueN(gem$sample)) stop("a sample maps to more than one GEM group")
in_obj <- colnames(obj)

read_track <- function(track) {
  out <- vector("list", nrow(gem))
  n_read <- n_drop_filter <- n_off_object <- 0L
  for (i in seq_len(nrow(gem))) {
    s <- gem$sample[i]
    if (track == "blast") {
      f <- file.path(RESCUE, paste0(s, "_rescued_reads.tsv"))
      d <- fread(f, colClasses = "character")
      n_read <- n_read + nrow(d)
      keep <- d[in_called_cell == "1" & nzchar(CB) & nzchar(UB)]
      n_drop_filter <- n_drop_filter + (nrow(d) - nrow(keep))
      d <- keep[, .(feature = norm_feat(construct), CB, UB)]
    } else {
      f <- file.path(ANTI, paste0(s, "_antisense_reads.tsv"))
      d <- fread(f, colClasses = "character")
      n_read <- n_read + nrow(d)
      d <- d[nzchar(CB) & nzchar(UB), .(feature = norm_feat(contig), CB, UB)]
    }
    if (!nrow(d)) next
    if (length(setdiff(d$feature, FEATS)))
      stop(s, " ", track, ": unknown feature ", paste(setdiff(d$feature, FEATS), collapse = ", "))
    # a molecule is a unique UB within a (barcode, feature), not a read
    d <- unique(d)
    d[, barcode := sub("-1$", paste0("-", gem$gem_group[i]), CB)]
    n_off_object <- n_off_object + d[!barcode %in% in_obj, .N]
    out[[i]] <- d[barcode %in% in_obj, .(n = .N), by = .(barcode, feature)]
  }
  out <- rbindlist(out)
  msg("  ", track, ": ", format(n_read, big.mark = ","), " reads -> ",
      format(sum(out$n), big.mark = ","), " molecules on ",
      format(uniqueN(out$barcode), big.mark = ","), " object nuclei (",
      format(n_drop_filter, big.mark = ","), " reads dropped by in_called_cell, ",
      format(n_off_object, big.mark = ","), " molecules on barcodes not in the object)")
  out
}

msg("folding in the rescued and antisense tracks")
for (track in c("blast", "antisense")) {
  d <- read_track(track)
  w <- dcast(d, barcode ~ feature, value.var = "n", fill = 0L)
  for (f in FEATS) {
    v <- if (f %in% names(w)) w[[f]] else rep(0L, nrow(w))
    add <- setNames(as.numeric(v), w$barcode)[cnt$barcode]
    add[is.na(add)] <- 0
    set(cnt, j = paste0(f, "__", track), value = add)
  }
}

# a barcode that Cell Ranger did not call cannot be rescued, so the filter above
# costs the object something; say how much rather than leave it implicit.
msg("cumulative UMIs per feature (base | +blast | +antisense):")
print(data.table(feature = FEATS, t(vapply(FEATS, function(f)
  cumsum(vapply(TRACKS, function(k) sum(cnt[[paste0(f, "__", k)]]), 0)),
  numeric(3)))))

## ---- classify, cumulatively over the three tracks -------------------------

fluor <- as.data.table(obj@meta.data)[, c("mouse", REGIONS), with = FALSE]
fluor <- unique(fluor)
bad <- fluor[, .(ok = setequal(unlist(.SD), ORFS)), by = mouse, .SDcols = REGIONS][ok == FALSE]
if (nrow(bad)) stop("fluorophore map is not a permutation of the three ORFs for: ",
                    paste(bad$mouse, collapse = ", "))
if (nrow(fluor) != uniqueN(fluor$mouse)) stop("more than one fluorophore map per animal")

# fluorophore -> region, per nucleus, off the object's own metadata
fluor_of <- vapply(REGIONS, function(r) as.character(obj@meta.data[[r]]),
                   character(ncol(obj)))

classify <- function(upto) {
  tot <- vapply(FEATS, function(f)
    rowSums(as.matrix(cnt[, paste0(f, "__", TRACKS[seq_len(upto)]), with = FALSE])),
    numeric(nrow(cnt)))
  det   <- tot[, ORFS, drop = FALSE] >= MIN_UMI
  n_orf <- as.integer(rowSums(det))
  orf   <- ifelse(n_orf == 1L, ORFS[max.col(det, ties.method = "first")], NA_character_)
  region <- rep(NA_character_, nrow(cnt))
  for (j in seq_along(REGIONS)) region[!is.na(orf) & orf == fluor_of[, j]] <- REGIONS[j]
  if (any(n_orf == 1L & is.na(region)))
    stop(sum(n_orf == 1L & is.na(region)),
         " nuclei have a single ORF that matches no region for their animal")
  data.table(n_orf = n_orf, region = region,
             class = fcase(n_orf == 1L, region,
                           n_orf >  1L, "multi_ORF",
                           tot[, SHARED] >= MIN_UMI, "shared_only",
                           default = "unlabelled"))
}

st <- lapply(seq_along(TRACKS), classify)
final <- st[[length(TRACKS)]]

# origin = the earliest track at which this nucleus became ORF-positive. It is
# about the nucleus, not the molecule: a nucleus already positive in base is
# `base` even if the later tracks add molecules to it.
origin <- rep(NA_character_, nrow(cnt))
for (i in rev(seq_along(TRACKS))) origin[st[[i]]$n_orf > 0L] <- TRACKS[i]
origin[is.na(final$region)] <- NA_character_

msg("class at each cumulative track:")
print(dcast(rbindlist(lapply(seq_along(TRACKS), function(i)
  st[[i]][, .(track = TRACKS[i], .N), by = class]), ),
  class ~ factor(track, levels = TRACKS), value.var = "N", fill = 0L)[order(-base)])
msg("assigned nuclei by the track that first made them ORF-positive:")
print(table(origin, useNA = "no"))
changed <- sum(!is.na(st[[1]]$region) & st[[1]]$region != final$class)
if (changed) msg("note: ", changed, " nuclei assigned in `base` change class once the ",
                 "added tracks push them over a second ORF (they become multi_ORF)")

calls <- data.table(
  barcode = colnames(obj),
  sample  = as.character(obj$sample),
  mouse   = factor(as.character(obj$mouse), levels = MOUSE_LEVELS),
  diet    = as.character(obj$diet),
  sex     = as.character(obj$sex),
  annot   = as.character(obj$annot),
  region  = final$region, class = final$class, origin = origin)

# one row per projection-positive nucleus: its target, and the ORF UMIs (three
# tracks summed) the call rests on. Unlabelled / multi_ORF / shared_only nuclei
# carry no target and are left out.
per_nuc <- calls[!is.na(region)]
per_nuc[, fluorophore := fcase(region == "BLA", as.character(obj$BLA[match(barcode, colnames(obj))]),
                               region == "NAc", as.character(obj$NAc[match(barcode, colnames(obj))]),
                               region == "LHr", as.character(obj$LHr[match(barcode, colnames(obj))]))]
for (f in ORFS)
  per_nuc[, paste0("umi_", f) := rowSums(as.matrix(
    cnt[match(per_nuc$barcode, cnt$barcode), paste0(f, "__", TRACKS), with = FALSE]))]
per_nuc[, class := NULL]
setnames(per_nuc, c("region", "origin"), c("projection_target", "first_track"))
setorder(per_nuc, mouse, projection_target, barcode)
fwrite(per_nuc, file.path(OUTDIR, paste0("projection_calls_per_nucleus", TAG, ".csv")))
msg("  projection_calls_per_nucleus", TAG, ".csv (", nrow(per_nuc), " nuclei)")

## ---- fig 27: UMAP coloured by projection target ---------------------------

# The origin is deliberately absent here (asked for, 2026-09-08): a nucleus is
# drawn as its target, however its molecules were found.
um <- as.data.table(Embeddings(obj, "umap"), keep.rownames = "barcode")
setnames(um, 2:3, c("UMAP_1", "UMAP_2"))
stopifnot(identical(um$barcode, calls$barcode))
um[, region := calls$region]
hit <- um[!is.na(region)]
hit[, region := factor(region, levels = REGIONS)]
setorder(hit, -region)   # rarest target drawn last, so it is not buried

p27 <- ggplot(um, aes(UMAP_1, UMAP_2)) +
  geom_point(colour = "grey88", size = 0.30, shape = 16) +
  geom_point(data = hit, aes(colour = region), size = 0.85, shape = 16) +
  scale_colour_manual(
    values = REGION_COLS, name = NULL, drop = FALSE,
    labels = sprintf("%s  (n = %s)", REGIONS,
                     format(as.integer(table(hit$region)[REGIONS]), big.mark = ",")),
    guide = guide_legend(override.aes = list(size = 3))) +
  coord_equal() + theme_bw(base_size = 11)
save_fig(p27, paste0("fig_27_umap_projection", TAG), 8.6, 8)

## ---- the table ------------------------------------------------------------

# A complete rectangle -- 12 animals x 3 targets x every cell type, zeros
# included -- so the table can be pivoted without an absent combination being
# read as a missing measurement rather than as a zero. n_nuclei splits into the
# three n_from_* columns by construction.
ct_levels <- c(names(CELLTYPE_COLS),
               setdiff(sort(unique(calls$annot)), names(CELLTYPE_COLS)))
ct_levels <- ct_levels[ct_levels %in% unique(calls$annot)]

animals <- unique(calls[, .(mouse, sample, diet, sex)])
animals <- merge(animals, fluor, by = "mouse")

per_animal <- calls[, .(n_nuclei_animal = .N,
                        n_multi_ORF     = sum(class == "multi_ORF"),
                        n_shared_only   = sum(class == "shared_only"),
                        n_projection_animal = sum(!is.na(region))), by = mouse]
per_region <- calls[!is.na(region), .N, by = .(mouse, region)]
per_cell   <- dcast(calls[!is.na(region)], mouse + region + annot ~ origin,
                    value.var = "barcode", fun.aggregate = length, fill = 0L)
for (k in TRACKS) if (!k %in% names(per_cell)) per_cell[, (k) := 0L]
setnames(per_cell, TRACKS, paste0("n_from_", TRACKS))

tab <- CJ(mouse = factor(MOUSE_LEVELS, levels = MOUSE_LEVELS),
          region = REGIONS, annot = ct_levels, unique = TRUE)
tab <- merge(tab, animals,    by = "mouse", all.x = TRUE)
tab <- merge(tab, per_cell,   by = c("mouse", "region", "annot"), all.x = TRUE)
tab <- merge(tab, per_region, by = c("mouse", "region"), all.x = TRUE)
setnames(tab, "N", "n_region_animal")
tab <- merge(tab, per_animal, by = "mouse", all.x = TRUE)
for (k in c(paste0("n_from_", TRACKS), "n_region_animal"))
  tab[is.na(get(k)), (k) := 0L]
tab[, n_nuclei := n_from_base + n_from_blast + n_from_antisense]

tab[, fluorophore := fifelse(region == "BLA", BLA, fifelse(region == "NAc", NAc, LHr))]
tab[, pct_of_region_animal := fifelse(n_region_animal > 0,
                                      round(100 * n_nuclei / n_region_animal, 2), NA_real_)]
tab[, annot  := factor(annot, levels = ct_levels)]
tab[, region := factor(region, levels = REGIONS)]
setorder(tab, mouse, region, annot)

setcolorder(tab, c("mouse", "sample", "diet", "sex", "region", "fluorophore",
                   "annot", "n_nuclei", "n_from_base", "n_from_blast",
                   "n_from_antisense", "pct_of_region_animal", "n_region_animal",
                   "n_projection_animal", "n_multi_ORF", "n_shared_only",
                   "n_nuclei_animal"))
tab[, c(REGIONS) := NULL]
stopifnot(sum(tab$n_nuclei) == calls[!is.na(region), .N])
fwrite(tab, file.path(OUTDIR, paste0("projection_by_celltype_per_animal", TAG, ".csv")))
msg("  projection_by_celltype_per_animal", TAG, ".csv (", nrow(tab), " rows)")

msg("projection-positive nuclei per animal:")
print(merge(animals[, .(mouse, diet, sex)], per_animal, by = "mouse")[order(mouse)])
msg("projection-positive nuclei per cell type and target:")
print(dcast(calls[!is.na(region)], annot ~ region, value.var = "barcode",
            fun.aggregate = length, fill = 0L)[order(-BLA)])

writeLines(capture.output(sessionInfo()),
           file.path(OUTDIR, paste0("sessionInfo_projection", TAG, ".txt")))
msg("done: ", OUTDIR)
