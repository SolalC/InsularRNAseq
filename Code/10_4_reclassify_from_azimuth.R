#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10_4_reclassify_from_azimuth.R  --  revise `annot` where the Azimuth mapping
# disagrees with the manual call. Run by 10_4_reclassify_from_azimuth.sh.
#
# THE PROBLEM. `annot` was called by hand in 7_annotation_featureplots.R on the
# 42-cluster embedding, and the two neuronal lineages were then subclustered and
# mapped onto the Allen mouse MOp reference in 10_2. That mapping is per
# NUCLEUS, and it disagrees with the manual call for a minority of them: nuclei
# sitting inside the excitatory or inhibitory blob that Azimuth assigns to
# Astro, Oligo, Micro-PVM, Endo, VLMC, OPC, Peri or SMC -- i.e. not neurons at
# all -- plus nuclei assigned to the OTHER neuronal lineage. 10_3 dropped both
# kinds from the DEG; this script gives them the label Azimuth gives them so
# they are reclassified rather than merely excluded.
#
# WHY THE PRIOR LABEL NEEDS NO LOOKUP. The excitatory prediction file contains
# exactly the nuclei `annot` called "Excitatory Neurons" and the inhibitory file
# exactly those it called "Inhibitory Neurons" -- that is what 10_neuronal_
# subclustering.R subset on. So the prior label is implied by which file a
# barcode is in, and the 17 GB annotated object never has to be read here.
#
# HOW A NUCLEUS IS REASSIGNED. By its subclass's MODAL predicted.class, not by
# its own predicted.class -- the two levels are transferred independently, so a
# per-nucleus class can disagree with its own subclass and would split a
# coherent group. This is the same rule 10_3 filters on, so the two agree by
# construction.
#
#   subclass modally Glutamatergic -> "Excitatory Neurons"
#   subclass modally GABAergic     -> "Inhibitory Neurons"
#   subclass modally Non-Neuronal  -> NON_NEURONAL_MAP below
#
# Micro-PVM is the one genuinely ambiguous case: the reference subclass pools
# microglia with perivascular macrophages, and `annot` carries both Microglia
# and Macrophage. It maps to Microglia, the dominant component. If that matters
# for a downstream claim, split it on Mrc1/Pf4/F13a1 rather than trusting this.
#
# WHAT IT WRITES. Into Results/reporting, prefixed `reclass_`: a per-barcode
# CSV, the move tables, and one figure. The canonical object
# HFDseq_v2_annotated.rds is NOT modified -- it holds a hand-made annotation
# that took a plot-by-plot pass to produce, and silently rewriting it would
# destroy work that cannot be regenerated from a script. Consumers join the CSV
# on barcode; figures_qc_to_annotation.R does this automatically when the file
# exists. To bake it in, assign the column and re-save the object deliberately.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
OUTDIR  <- file.path(PROJECT, "Results", "reporting")   # asked for, 2026-09-07
RES     <- 0.2

LINEAGE_PRIOR <- c(excitatory = "Excitatory Neurons",
                   inhibitory = "Inhibitory Neurons")
CLASS_TO_ANNOT <- c(Glutamatergic = "Excitatory Neurons",
                    GABAergic     = "Inhibitory Neurons")
NON_NEURONAL_MAP <- c(
  Astro       = "Astrocytes",
  Oligo       = "Oligodendrocytes",
  OPC         = "OPC",
  `Micro-PVM` = "Microglia",
  Endo        = "Endothelial Cell",
  VLMC        = "VLCM",
  Peri        = "Pericytes",
  SMC         = "Smooth Muscle")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
source(file.path(PROJECT, "Code", "lib_celltype_palette.R"))

## ---- read both mappings -----------------------------------------------------

read_lineage <- function(lin) {
  f <- file.path(PROJECT, "Results", sprintf("10_2_%s_azimuth", lin),
                 sprintf("azimuth_predictions_res%s.rds", RES))
  if (!file.exists(f))
    stop("no Azimuth predictions for ", lin, " at ", f,
         " -- run 10_2_azimuth.sh ", lin, " first")
  d <- as.data.table(readRDS(f))
  d[, .(barcode, lineage = lin, mouse, sample, diet, sex,
        subclass = as.character(predicted.subclass),
        pred_class = as.character(predicted.class),
        subclass_score = predicted.subclass.score,
        mapping_score = mapping.score)]
}
pred <- rbindlist(lapply(names(LINEAGE_PRIOR), read_lineage))
pred[, annot_prior := LINEAGE_PRIOR[lineage]]
msg("read ", format(nrow(pred), big.mark = ","), " mapped nuclei (",
    paste(sprintf("%s %s", names(table(pred$lineage)),
                  format(as.integer(table(pred$lineage)), big.mark = ",")),
          collapse = ", "), ")")

## ---- the modal class of each subclass ---------------------------------------

modal <- pred[, .N, by = .(subclass, pred_class)][order(subclass, -N)][
  , .SD[1], by = subclass][, .(subclass, modal_class = pred_class, n_modal = N)]
sz <- pred[, .N, by = subclass]
modal <- merge(modal, sz, by = "subclass")[order(-N)]
fwrite(modal, file.path(OUTDIR, "reclass_subclass_modal_class.csv"))

resolve <- function(subclass, modal_class) {
  out <- rep(NA_character_, length(subclass))
  neuro <- modal_class %in% names(CLASS_TO_ANNOT)
  out[neuro] <- CLASS_TO_ANNOT[modal_class[neuro]]
  nn <- !neuro & subclass %in% names(NON_NEURONAL_MAP)
  out[nn] <- NON_NEURONAL_MAP[subclass[nn]]
  out
}
modal[, annot_azimuth := resolve(subclass, modal_class)]

unmapped <- modal[is.na(annot_azimuth)]
if (nrow(unmapped))
  stop("no annot label for subclass(es): ",
       paste(sprintf("%s (modal %s, n=%d)", unmapped$subclass,
                     unmapped$modal_class, unmapped$N), collapse = ", "),
       " -- add them to NON_NEURONAL_MAP rather than letting them fall through")

pred <- merge(pred, modal[, .(subclass, modal_class, annot_azimuth)],
              by = "subclass", sort = FALSE)
pred[, changed := annot_azimuth != annot_prior]

## ---- report and write -------------------------------------------------------

msg(format(sum(pred$changed), big.mark = ","), " of ",
    format(nrow(pred), big.mark = ","), " mapped nuclei change label (",
    sprintf("%.2f%%", 100 * mean(pred$changed)), ")")

moves <- pred[changed == TRUE, .N, by = .(annot_prior, annot_azimuth, subclass)][order(-N)]
msg("--- what moves, and where to ---")
print(moves)
fwrite(moves, file.path(OUTDIR, "reclass_moves.csv"))

summ <- pred[changed == TRUE, .(n = .N,
                                median_subclass_score = round(median(subclass_score), 3),
                                median_mapping_score = round(median(mapping_score), 3)),
             by = .(annot_prior, annot_azimuth)][order(-n)]
msg("--- summary, with the confidence of the calls that move ---")
print(summ)
fwrite(summ, file.path(OUTDIR, "reclass_summary.csv"))

fwrite(pred[, .(barcode, lineage, mouse, sample, diet, sex, subclass, modal_class,
                subclass_score, mapping_score, annot_prior, annot_azimuth, changed)],
       file.path(OUTDIR, "reclass_annot_revised.csv.gz"), compress = "gzip")
msg("wrote reclass_annot_revised.csv.gz -- the canonical object was NOT modified")

## ---- figure: where the reassigned nuclei come from --------------------------

pm <- pred[changed == TRUE, .N, by = .(mouse, annot_azimuth)]
pm[, mouse := factor(mouse, levels = paste0("Ech", 1:12))]
acols <- celltype_cols(sort(unique(pm$annot_azimuth)))
p <- ggplot(pm, aes(mouse, N, fill = annot_azimuth)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = acols, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "nuclei reassigned from a neuronal label") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major.x = element_blank())
ggsave(file.path(OUTDIR, "reclass_per_animal.png"), p,
       width = 10, height = 5, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo_reclass.txt"))
msg("done: ", OUTDIR)
