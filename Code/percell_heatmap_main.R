# Main heatmap of the per-cell DEG (HFD vs SD), all cell types: per-mouse expression of the
# hits for SD mice (left) and HFD mice (right), the per-cell avg_log2FC in the middle.
# Reads only the CSVs written by percell_deg_diet_sexfollowup.R; no Seurat object needed.
# Same block as the "main heatmap" section of percell_deg_plots.R.
# Usage: Rscript percell_heatmap_main.R   (or paste into R)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
DEG_DIR    <- "/QRISdata/Q9732/scratch/Results/percell_wilcox_diet/all_celltypes"
BONF_CUT   <- 0.05
LFC_CUT    <- 0.5
HEAT_MAX   <- 50
DROP_TYPES <- "Unidentified"
safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
msg  <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## ---- main heatmap: per-mouse expression (SD | log2FC | HFD) ------------------------
# One figure per cell type: left, the SD mice; right, the HFD mice; middle, the per-cell
# DEG avg_log2FC (HFD / SD) of each gene. Rows: the HEAT_MAX hits with the lowest p_val,
# sorted by avg_log2FC (most up at the top).
# Expression: mouse_means_hits_<ct>.csv (DEG script), the mean of expm1(RNA data) over the
#   mouse's nuclei of that cell type, i.e. the same per-mouse average Seurat's avg_log2FC is
#   built from. Shown as log2(mean + 1), then per gene as set by MAIN_SCALE:
#     "zSD"  : (x - mean of SD mice) / sd of SD mice   (SD mice centre on 0)
#     "zall" : (x - mean of all mice) / sd of all mice (conventional row scaling)
#     "log2" : log2(mean + 1) unscaled (colour then mostly shows gene abundance)
# Output: DEG_DIR/heatmap_main/heatmap_main_<ct>.jpeg
MAIN_SCALE <- "zSD"
MAIN_CLIP  <- 3      # colour limit for the z panels
MAIN_DIR   <- file.path(DEG_DIR, "heatmap_main")
dir.create(MAIN_DIR, showWarnings = FALSE, recursive = TRUE)

main_files <- list.files(DEG_DIR, pattern = "^mouse_means_hits_.*\\.csv$", full.names = TRUE)
main_files <- main_files[!basename(main_files) %in% paste0("mouse_means_hits_", safe(DROP_TYPES), ".csv")]
for (fm in main_files) {
  ct_safe <- sub("^mouse_means_hits_(.*)\\.csv$", "\\1", basename(fm))
  fd <- file.path(DEG_DIR, paste0("DEG_percell_HFDvsSD_", ct_safe, ".csv"))
  if (!file.exists(fd)) next
  deg <- fread(fd)
  hits <- deg[p_bonferroni < BONF_CUT & abs(avg_log2FC) > LFC_CUT][order(p_val)]
  mm   <- fread(fm)
  hits <- head(hits[gene %in% mm$gene], HEAT_MAX)
  if (nrow(hits) < 2) next
  gene_lv <- hits[order(avg_log2FC), gene]               # bottom to top: most up at the top
  mm <- mm[gene %in% hits$gene]
  mm[, `:=`(diet = factor(diet, levels = c("SD", "HFD")), sex = factor(sex, levels = c("F", "M")))]
  mm[, x := log2(mean_expr + 1)]
  mm[, value := switch(MAIN_SCALE,
    zSD  = { s <- sd(x[diet == "SD"]);  if (is.na(s) || s == 0) sign(x - mean(x[diet == "SD"])) * MAIN_CLIP else (x - mean(x[diet == "SD"])) / s },
    zall = { s <- sd(x);                if (is.na(s) || s == 0) 0 else (x - mean(x)) / s },
    log2 = x), by = gene]
  mice_ord <- unique(mm[, .(sample, diet, sex)])[order(diet, sex, sample)]
  mm[, `:=`(gene = factor(gene, levels = gene_lv),
            col  = factor(sample, levels = mice_ord$sample))]
  xlab_of <- setNames(paste0(mice_ord$sample, " (", mice_ord$sex, ")"), mice_ord$sample)

  expr_fill <- if (MAIN_SCALE == "log2")
    scale_fill_viridis_c(name = "log2(mean + 1)") else
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", limits = c(-MAIN_CLIP, MAIN_CLIP),
                         oob = scales::squish,
                         name = paste0(if (MAIN_SCALE == "zSD") "z vs SD mice" else "z across mice",
                                       "\nlog2 mean expr\n(clipped +/-", MAIN_CLIP, ")"))
  base_th <- theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "top",
          legend.title = element_text(size = 7), legend.key.height = unit(6, "pt"),
          plot.margin = margin(2, 1, 2, 1))
  panel_expr <- function(dt, title, show_y) {
    ggplot(mm[diet == dt], aes(col, gene, fill = value)) + geom_tile() +
      scale_x_discrete(labels = xlab_of) + expr_fill + labs(title = title) + base_th +
      theme(axis.text.y = if (show_y) element_text(size = 6) else element_blank())
  }
  lfc_lim <- max(abs(hits$avg_log2FC))
  p_lfc <- ggplot(hits[, .(gene = factor(gene, levels = gene_lv), avg_log2FC, col = "log2FC")],
                  aes(col, gene, fill = avg_log2FC)) +
    geom_tile() + geom_text(aes(label = sprintf("%.2f", avg_log2FC)), size = 1.8) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         limits = c(-lfc_lim, lfc_lim), name = "log2FC\nHFD / SD") +
    labs(title = "log2FC") + base_th + theme(axis.text.y = element_blank())
  pg <- lapply(list(panel_expr("SD", "SD", TRUE), p_lfc, panel_expr("HFD", "HFD", FALSE)),
               ggplotGrob)
  n_col <- c(sum(mice_ord$diet == "SD"), 1.6, sum(mice_ord$diet == "HFD"))   # panel widths, in tiles
  for (i in 1:3) {
    pc <- unique(pg[[i]]$layout[grepl("^panel", pg[[i]]$layout$name), "l"])
    pg[[i]]$widths[pc] <- unit(n_col[i], "null")
  }
  g <- cbind(pg[[1]], pg[[2]], pg[[3]], size = "max")
  ttl <- sprintf("%s: HFD vs SD, %d per-cell hits (p_bonferroni < %g, |log2FC| > %g)",
                 gsub("_", " ", ct_safe), nrow(hits), BONF_CUT, LFC_CUT)
  g <- gtable::gtable_add_rows(g, unit(14, "pt"), pos = 0)
  g <- gtable::gtable_add_grob(g, grid::textGrob(ttl, gp = grid::gpar(fontsize = 9, fontface = "bold")),
                               t = 1, l = 1, r = ncol(g))
  jpeg(file.path(MAIN_DIR, paste0("heatmap_main_", ct_safe, ".jpeg")),
       width = 3 + 0.3 * (nrow(mice_ord) + 1.6) + 1.5, height = 2 + 0.13 * nrow(hits),
       units = "in", res = 300)
  grid::grid.newpage(); grid::grid.draw(g); invisible(dev.off())
  msg("main heatmap: ", ct_safe, " (", nrow(hits), " genes)")
}
