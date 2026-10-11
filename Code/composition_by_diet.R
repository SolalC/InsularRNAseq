#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# composition_by_diet.R  --  fig_23 (cell-type composition) summarised by
# condition instead of by animal: diet, then diet x sex. Run by composition_by_diet.sh. Asked 2026-09-23.
#
# Housekeeping, no number. It reads annotation_composition_per_animal.csv, the
# table figures_qc_to_annotation.R writes alongside fig_23 from the SAME
# Azimuth-revised `annot`, so the numbers are fig_23's; each condition bar is
# the mean of its animals' bars. The 17 GB object is not read.
#
# The csv's column order is the factor order of `annot` (dcast keeps it), which
# is also fig_23's stacking order, so it is kept as-is.
#
# Outputs, in Results/reporting -- one bar per condition, animals averaged:
#   fig_23_annotation_composition_diet.png      SD | HFD
#   fig_23_annotation_composition_diet_sex.png  SD F | SD M | HFD F | HFD M
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
OUTDIR  <- file.path(PROJECT, "Results", "reporting")
source(file.path(PROJECT, "Code", "lib_celltype_palette.R"))

w <- fread(file.path(OUTDIR, "annotation_composition_per_animal.csv"))
ct <- setdiff(names(w), c("mouse", "diet", "sex"))
comp <- melt(w, id.vars = c("mouse", "diet", "sex"), variable.name = "annot",
             value.name = "pct", variable.factor = FALSE)
comp[, annot := factor(annot, levels = ct)]
comp[, mouse := factor(mouse, levels = paste0("Ech", 1:12))]
comp[, diet  := factor(diet,  levels = c("SD", "HFD"))]
comp[, sex   := factor(sex,   levels = c("F", "M"))]
chk <- comp[, sum(pct), by = mouse]
if (any(abs(chk$V1 - 100) > 1e-6)) stop("an animal's percentages do not sum to 100")
print(unique(comp[, .(mouse, diet, sex)])[order(diet, sex, mouse)])

# One bar per condition. Each animal's percentages are averaged, so every
# animal weighs the same whatever its nucleus count (the animal is the unit,
# n = 3 per diet x sex, 6 per diet) -- nuclei are NOT pooled.
acols <- celltype_cols(ct)
cond_bar <- function(by, file) {
  m <- comp[, .(pct = mean(pct), n_animals = .N), by = c(by, "annot")]
  m[, condition := factor(do.call(paste, .SD), levels = unique(do.call(paste, .SD)[order(diet, if ("sex" %in% by) sex else diet)])),
    .SDcols = by]
  print(dcast(m, annot ~ condition, value.var = "pct")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
  n <- unique(m[, .(condition, n_animals)])
  p <- ggplot(m, aes(condition, pct, fill = annot)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = acols, name = NULL) +
    scale_x_discrete(labels = setNames(sprintf("%s\n(n = %d animals)", n$condition, n$n_animals),
                                       n$condition)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "% of nuclei (mean of animals)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.x = element_blank())
  ggsave(file.path(OUTDIR, file), p, width = 3 + 1.1 * nrow(n), height = 6, dpi = 300)
}
cond_bar("diet",          "fig_23_annotation_composition_diet.png")
cond_bar(c("diet", "sex"), "fig_23_annotation_composition_diet_sex.png")
cat("done\n")
