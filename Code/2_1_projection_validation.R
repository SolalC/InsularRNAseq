#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 2_1_projection_validation.R
#
# Reads the cellranger count matrices for all 12 HFDseq samples, counts how many
# called cells carry each transgene above a UMI threshold, and translates the
# fluorophore calls into projection regions (BLA / NAc / LHr).
#
# The fluorophore -> region assignment is RANDOMISED PER SAMPLE, so the mapping
# is read per sample from Data/sample_metadata.csv (columns BLA, NAc, LHr each
# hold the fluorophore used for that region in that mouse). Never assume a fixed
# fluorophore -> region correspondence.
#
# Only the four transgene rows are pulled out of each matrix.mtx.gz (via awk),
# so this never loads the full 33,700 x ~28,000 matrix into memory.
#
# Reads either cellranger run, selected by RUN below:
#   v2  Data/Counts_v2, reference GRCm39_cassette_v2. Reporter features are
#       whole-cassette contigs and the shared element is "AAV_shared".
#   v1  Data/Counts, reference GRCm39_transgenes, shared element
#       "AAV_WPRE_SV40pA_shared". Kept so the two runs can be compared.
# ---------------------------------------------------------------------------

# data.table only - nothing here uses tidyverse (no pipes, no dplyr verbs), and
# dropping it removes a dependency that is not present in the r/4.4.x modules.
library(data.table)

## ---- CONFIG ---------------------------------------------------------------

# Which cellranger run to read: "v2" (cassette reference) or "v1" (original).
RUN <- "v2"

# A cell counts as positive when its UMI count is >= this. THIS IS INCLUSIVE:
# MIN_UMI <- 1 means "one or more UMIs", which is the presence/absence call the
# projection design rests on. The earlier version of this script compared with
# `>`, so MIN_UMI <- 1 there silently meant ">= 2 UMIs" - output files from that
# version are tagged "minumi<N>", these are tagged "ge<N>", so the two cannot be
# confused.
MIN_UMI <- 1

PROJECT <- "/scratch/user/uqschauq/Insular"
META    <- file.path(PROJECT, "Data", "sample_metadata.csv")

if (RUN == "v2") {
  COUNTS <- file.path(PROJECT, "Data", "Counts_v2")
  OUTDIR <- file.path(PROJECT, "Results", "2_1_projection_validation_v2")
  SHARED <- "AAV_shared"
} else if (RUN == "v1") {
  COUNTS <- file.path(PROJECT, "Data", "Counts")
  OUTDIR <- file.path(PROJECT, "Results", "2_1_projection_validation")
  SHARED <- "AAV_WPRE_SV40pA_shared"
} else {
  stop("RUN must be \"v1\" or \"v2\", got: ", RUN)
}

# >= 0 would call every cell positive; the comparison is inclusive now.
if (MIN_UMI < 1) stop("MIN_UMI must be >= 1 (the comparison is inclusive)")

# "filtered_feature_bc_matrix" = cellranger-called cells (what we want here).
# "raw_feature_bc_matrix" = all barcodes, useful for ambient-transgene estimates.
MATRIX <- "filtered_feature_bc_matrix"

# The three projection reporters. The shared AAV element is deliberately NOT one
# of them: it carries no projection information, it is only a transgene-presence
# readout used to sanity-check the ORF threshold.
# In v2 these are gene_ids on the whole-cassette contigs (EGFP_cassette etc.),
# so the names are unchanged; SHARED is set with the run above.
ORFS   <- c("EGFP", "tdTomato", "EBFP2")

REGIONS <- c("BLA", "NAc", "LHr")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---- SAMPLE METADATA / RANDOMISED FLUOROPHORE MAP -------------------------

# sample_metadata.csv has CRLF line endings: strip \r from names and values or
# string comparisons silently fail on the last column.
strip_cr <- function(x) if (is.character(x)) gsub("\r", "", x, fixed = TRUE) else x

meta <- fread(META)
setnames(meta, gsub("\r", "", names(meta), fixed = TRUE))
meta[, names(meta) := lapply(.SD, strip_cr)]

# Metadata writes the blue reporter "EBFP"; the reference contig is "EBFP2".
norm_fluor <- function(x) {
  x <- trimws(gsub("\r", "", x, fixed = TRUE))
  ifelse(x == "EBFP", "EBFP2", x)
}

fluor_map <- meta[, .(
  sample = fastq_prefix,
  mouse  = sample_id,
  sex    = sex,
  diet   = diet,
  BLA    = norm_fluor(BLA),
  NAc    = norm_fluor(NAc),
  LHr    = norm_fluor(LHr)
)]

fluor_long <- melt(
  fluor_map,
  id.vars       = c("sample", "mouse", "sex", "diet"),
  measure.vars  = REGIONS,
  variable.name = "region",
  value.name    = "fluor"
)
fluor_long[, region := as.character(region)]

# Guard: every sample must use each of the three ORFs exactly once, otherwise
# the region translation below would be silently wrong.
bad <- fluor_long[, .(ok = setequal(fluor, ORFS)), by = sample][ok == FALSE]
if (nrow(bad) > 0) {
  stop("Fluorophore map is not a permutation of the three ORFs for: ",
       paste(bad$sample, collapse = ", "))
}

samples <- fluor_map$sample

cat("Fluorophore -> region assignment per sample (randomised):\n")
print(fluor_map[, .(sample, sex, diet, BLA, NAc, LHr)])
cat("\nRun:", RUN, "  counts:", COUNTS, "\n")
cat("Detection threshold: count >=", MIN_UMI, "UMIs\n")
cat("Matrix:", MATRIX, "\n\n")

## ---- PER-SAMPLE TRANSGENE EXTRACTION --------------------------------------

# Pull only the transgene rows out of the MatrixMarket file. Row indices are
# looked up per sample from features.tsv.gz rather than hardcoded, so this stays
# correct if the reference is ever rebuilt with a different contig order.
read_transgene_counts <- function(sample) {

  dir <- file.path(COUNTS, sample, "outs", MATRIX)
  f_features <- file.path(dir, "features.tsv.gz")
  f_matrix   <- file.path(dir, "matrix.mtx.gz")

  stopifnot(file.exists(f_features), file.exists(f_matrix))

  features <- fread(cmd = paste("zcat", shQuote(f_features)),
                    header = FALSE, sep = "\t", select = 1:2,
                    col.names = c("feature_id", "feature_name"))
  features[, row := .I]

  want <- features[feature_id %in% c(ORFS, SHARED)]
  if (nrow(want) != length(ORFS) + 1L) {
    stop(sample, ": expected ", length(ORFS) + 1L,
         " transgene features, found ", nrow(want))
  }

  # Total called cells = number of columns, from the MatrixMarket dims line.
  dims <- fread(cmd = paste("zcat", shQuote(f_matrix), "| sed -n '3p'"),
                header = FALSE)
  n_cells <- as.integer(dims$V2)

  # Stream the matrix, keeping only the transgene rows.
  cond <- paste(sprintf("$1==%d", want$row), collapse = " || ")
  hits <- fread(
    cmd = paste0("zcat ", shQuote(f_matrix), " | awk 'NR>3 && (", cond, ")'"),
    header = FALSE, col.names = c("row", "col", "count")
  )

  # Cells with no transgene UMI at all never appear in `hits`; they are added
  # back as zeros via n_cells when the summary is assembled.
  if (nrow(hits) == 0) {
    wide <- data.table(col = integer(0))
    for (g in c(ORFS, SHARED)) wide[, (g) := integer(0)]
  } else {
    hits <- merge(hits, want[, .(row, feature_id)], by = "row")
    wide <- dcast(hits, col ~ feature_id, value.var = "count",
                  fill = 0L, fun.aggregate = sum)
    for (g in c(ORFS, SHARED)) if (!g %in% names(wide)) wide[, (g) := 0L]
  }

  list(sample = sample, n_cells = n_cells, wide = wide)
}

## ---- CLASSIFY CELLS -------------------------------------------------------

# Each cell falls into exactly one class:
#   BLA / NAc / LHr  - exactly one ORF over threshold, translated via the map
#   multi_ORF        - two or three ORFs over threshold (ambiguous, not counted
#                      into any region; watch this, EGFP and EBFP2 share a 98 bp
#                      identical stretch and can cross-map)
#   shared_only      - shared AAV element over threshold but no ORF:
#                      "transgene-positive, projection-unassigned"
#   unlabelled       - nothing over threshold
classify_sample <- function(res) {

  this_sample <- res$sample
  wide        <- res$wide
  n_cells     <- res$n_cells

  smap <- fluor_long[sample == this_sample]
  fluor_to_region <- setNames(smap$region, smap$fluor)

  per_cell <- data.table(col = wide$col)
  for (g in c(ORFS, SHARED)) per_cell[, (g) := wide[[g]]]

  det <- as.matrix(per_cell[, ..ORFS]) >= MIN_UMI
  per_cell[, n_orf := as.integer(rowSums(det))]
  per_cell[, shared_pos := get(SHARED) >= MIN_UMI]

  # Name of the single ORF that passed, when exactly one did.
  if (nrow(per_cell)) {
    idx <- max.col(det, ties.method = "first")
    per_cell[, orf := ifelse(n_orf == 1L, ORFS[idx], NA_character_)]
  } else {
    per_cell[, orf := character(0)]
  }
  per_cell[, region := fluor_to_region[orf]]

  per_cell[, class := fcase(
    n_orf == 1L, region,
    n_orf >  1L, "multi_ORF",
    shared_pos,  "shared_only",
    default = "unlabelled"
  )]

  n_zero_cell <- n_cells - nrow(per_cell)  # cells with no transgene UMI at all

  tally <- per_cell[, .N, by = class]
  get_n <- function(k) { v <- tally[class == k, N]; if (length(v)) v else 0L }

  summary_row <- data.table(
    sample      = this_sample,
    total_cells = n_cells,
    BLA         = get_n("BLA"),
    NAc         = get_n("NAc"),
    LHr         = get_n("LHr"),
    multi_ORF   = get_n("multi_ORF"),
    shared_only = get_n("shared_only"),
    unlabelled  = get_n("unlabelled") + n_zero_cell
  )
  summary_row[, assigned := BLA + NAc + LHr]
  summary_row[, pct_assigned := round(100 * assigned / total_cells, 2)]

  # Per-fluorophore detection, independent of the region translation - useful
  # for spotting a reporter that simply failed in a given mouse.
  fluor_row <- data.table(sample = this_sample)
  for (g in c(ORFS, SHARED)) fluor_row[, (g) := sum(per_cell[[g]] >= MIN_UMI)]

  # Which ORF combinations drive the multi_ORF calls.
  combos <- data.table()
  if (any(per_cell$n_orf > 1L)) {
    mc  <- det[per_cell$n_orf > 1L, , drop = FALSE]
    lab <- apply(mc, 1, function(r) paste(ORFS[r], collapse = "+"))
    combos <- data.table(sample = this_sample, combo = lab)[, .N, by = .(sample, combo)]
  }

  per_cell[, sample := this_sample]
  setnames(per_cell, SHARED, "shared")

  list(summary = summary_row, fluor = fluor_row, combos = combos,
       per_cell = per_cell[, .(sample, col, EGFP, tdTomato, EBFP2, shared,
                               n_orf, orf, region, class)])
}

## ---- PREFLIGHT ------------------------------------------------------------

# Fail before reading anything if the cellranger run is not finished, rather
# than part way through the loop with half the samples already processed.
missing <- samples[!file.exists(file.path(COUNTS, samples, "outs",
                                          "metrics_summary.csv"))]
if (length(missing)) {
  stop("run ", RUN, ": no completed cellranger output for ",
       length(missing), " sample(s): ", paste(missing, collapse = ", "),
       "\n  expected ", file.path(COUNTS, "<sample>", "outs",
                                  "metrics_summary.csv"))
}
cat("Preflight OK:", length(samples), "completed count directories under",
    COUNTS, "\n\n")

## ---- RUN ------------------------------------------------------------------

all_summary <- list(); all_fluor <- list(); all_combos <- list(); all_cells <- list()

for (s in samples) {
  cat("Reading", s, "... ")
  res <- read_transgene_counts(s)
  cls <- classify_sample(res)
  all_summary[[s]] <- cls$summary
  all_fluor[[s]]   <- cls$fluor
  all_combos[[s]]  <- cls$combos
  all_cells[[s]]   <- cls$per_cell
  cat(res$n_cells, "cells,", cls$summary$assigned, "assigned\n")
}

summary_dt <- rbindlist(all_summary)
fluor_dt   <- rbindlist(all_fluor)
combos_dt  <- rbindlist(all_combos, fill = TRUE)
cells_dt   <- rbindlist(all_cells)

# Attach sex/diet and the per-sample fluorophore map to the summary.
summary_dt <- merge(
  fluor_map[, .(sample, mouse, sex, diet, BLA_fluor = BLA,
                NAc_fluor = NAc, LHr_fluor = LHr)],
  summary_dt, by = "sample", sort = FALSE
)
setcolorder(summary_dt, c("sample", "mouse", "sex", "diet", "total_cells",
                          "BLA", "NAc", "LHr", "assigned", "pct_assigned",
                          "multi_ORF", "shared_only", "unlabelled",
                          "BLA_fluor", "NAc_fluor", "LHr_fluor"))

## ---- REPORT ---------------------------------------------------------------

cat("\n=====================================================================\n")
cat("Cells per projection region, threshold >= ", MIN_UMI, " UMIs (run ", RUN, ")\n", sep = "")
cat("=====================================================================\n\n")
print(summary_dt[, .(sample, sex, diet, total_cells, BLA, NAc, LHr,
                     assigned, pct_assigned, multi_ORF, shared_only)])

cat("\n--- Raw per-fluorophore positives (before region translation) ---\n")
print(merge(fluor_dt, fluor_map[, .(sample, BLA, NAc, LHr)],
            by = "sample", sort = FALSE))

if (nrow(combos_dt)) {
  cat("\n--- multi_ORF breakdown (possible cross-mapping / doublets) ---\n")
  print(dcast(combos_dt, sample ~ combo, value.var = "N", fill = 0))
}

cat("\n--- Totals by region across all samples ---\n")
print(summary_dt[, .(
  BLA = sum(BLA), NAc = sum(NAc), LHr = sum(LHr),
  assigned = sum(assigned), multi_ORF = sum(multi_ORF),
  shared_only = sum(shared_only), total_cells = sum(total_cells)
)])

cat("\n--- Region totals split by diet and sex ---\n")
print(summary_dt[, .(BLA = sum(BLA), NAc = sum(NAc), LHr = sum(LHr),
                     total_cells = sum(total_cells)), by = .(diet, sex)])

## ---- WRITE ----------------------------------------------------------------

tag <- paste0("ge", MIN_UMI)
fwrite(summary_dt, file.path(OUTDIR, paste0("projection_summary_", tag, ".tsv")), sep = "\t")
fwrite(fluor_dt,   file.path(OUTDIR, paste0("fluorophore_positives_", tag, ".tsv")), sep = "\t")
fwrite(cells_dt,   file.path(OUTDIR, paste0("per_cell_calls_", tag, ".tsv.gz")), sep = "\t")
if (nrow(combos_dt)) {
  fwrite(combos_dt, file.path(OUTDIR, paste0("multi_orf_combos_", tag, ".tsv")), sep = "\t")
}

cat("\nWritten to", OUTDIR, "\n")
