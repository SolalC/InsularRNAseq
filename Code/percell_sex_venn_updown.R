# Venn of female vs male HFD-vs-SD DEG per cell type, split into up and down in HFD, from DEG_percell_sexsplit_<ct>.csv.
# Sex DEG: within-sex p_bonferroni < 0.05 (Bonferroni over the cell type's pooled hits) and
# |within-sex avg_log2FC| > 0.5. Universe: pooled hits (p_bonferroni < 0.05, |log2FC| > 0.5).
library(data.table); library(ggplot2)
DEG_DIR  <- "/QRISdata/Q9732/scratch/Results/percell_wilcox_diet/all_celltypes"
BONF_CUT <- 0.05; LFC_CUT <- 0.5

files <- list.files(DEG_DIR, pattern = "^DEG_percell_sexsplit_.*\\.csv$", full.names = TRUE)
files <- files[!grepl("_Unidentified\\.csv$", files)]
vn <- rbindlist(lapply(files, function(f) {
  d <- fread(f)
  d <- d[p_bonferroni < BONF_CUT & abs(avg_log2FC) > LFC_CUT]
  d[, .(annot = sub("^DEG_percell_sexsplit_(.*)\\.csv$", "\\1", basename(f)),
        gene, F_avg_log2FC, M_avg_log2FC, F_p_bonferroni, M_p_bonferroni)]
}))
vn[, `:=`(F_deg = !is.na(F_p_bonferroni) & F_p_bonferroni < BONF_CUT & abs(F_avg_log2FC) > LFC_CUT,
          M_deg = !is.na(M_p_bonferroni) & M_p_bonferroni < BONF_CUT & abs(M_avg_log2FC) > LFC_CUT)]
vn[, region := fifelse(F_deg & M_deg, "both", fifelse(F_deg, "F only", fifelse(M_deg, "M only", "neither")))]

# Four sets per cell type: F up, M up, F down, M down (up = log2FC > 0, higher in HFD).
# F up / F down (and M up / M down) cannot share a gene, so regions containing both are 0.
SETS <- c("F up", "M up", "F down", "M down")
vn[, `:=`(`F up`   = F_deg & F_avg_log2FC > 0, `F down` = F_deg & F_avg_log2FC < 0,
          `M up`   = M_deg & M_avg_log2FC > 0, `M down` = M_deg & M_avg_log2FC < 0)]

# Classic 4-ellipse Venn layout; every one of the 15 regions exists.
ell <- data.table(set = SETS, x0 = c(0.38, 0.47, 0.53, 0.62), y0 = c(0.45, 0.55, 0.55, 0.45),
                  a = 0.35, b = 0.2, rot = c(-45, -45, 45, 45) * pi / 180)
inside <- function(x, y, e) {   # scaled radius: < 1 inside the ellipse
  dx <- x - e$x0; dy <- y - e$y0
  u <- dx * cos(e$rot) + dy * sin(e$rot); v <- -dx * sin(e$rot) + dy * cos(e$rot)
  sqrt((u / e$a)^2 + (v / e$b)^2)
}
# Label position of each region: the grid point of that region farthest from any outline.
grid <- CJ(x = seq(0, 1, by = 0.004), y = seq(0, 1.05, by = 0.004))
r <- sapply(seq_len(nrow(ell)), function(i) inside(grid$x, grid$y, ell[i]))
grid[, key := apply(r < 1, 1, function(z) paste(as.integer(z), collapse = ""))]
grid[, depth := apply(abs(r - 1), 1, min)]
pos <- grid[key != "0000"][order(-depth), .SD[1], by = key][, .(key, x, y)]

vn[, key := paste0(as.integer(`F up`), as.integer(`M up`), as.integer(`F down`), as.integer(`M down`))]
cnt <- CJ(annot = unique(vn$annot), key = pos$key)
cnt <- merge(cnt, vn[key != "0000", .N, by = .(annot, key)], by = c("annot", "key"), all.x = TRUE)
cnt[is.na(N), N := 0L]
cnt <- merge(cnt, pos, by = "key")
strip <- vn[, .(lab = sprintf("%s\n%d pooled hits; %d DEG in neither sex", gsub("_", " ", annot[1]),
                              .N, sum(key == "0000"))), by = annot]
cnt[, facet := strip$lab[match(annot, strip$annot)]]

th <- seq(0, 2 * pi, length.out = 200)
outl <- ell[, {
  u <- a * cos(th); v <- b * sin(th)
  .(x = x0 + u * cos(rot) - v * sin(rot), y = y0 + u * sin(rot) + v * cos(rot))
}, by = set]
outl[, set := factor(set, levels = SETS)]
# set names just beyond each ellipse's upper end
lab_sets <- ell[, .(set, x = x0 - 1.1 * a * cos(rot) * sign(rot) * -1, y = y0 + 1.1 * a * abs(sin(rot)) + 0.04)]
gp <- ggplot() +
  geom_polygon(data = outl, aes(x, y, fill = set, group = set), alpha = 0.2, colour = "grey30", linewidth = 0.3) +
  geom_text(data = lab_sets, aes(x, y, label = set, colour = set), size = 2.8, fontface = "bold") +
  geom_text(data = cnt, aes(x, y, label = N), size = 2.6) +
  facet_wrap(~ facet, ncol = 4) +
  scale_fill_manual(values = c("F up" = "#d95f02", "M up" = "#1b9e77", "F down" = "#7570b3", "M down" = "#e7298a"),
                    guide = "none", aesthetics = c("fill", "colour")) +
  coord_equal(xlim = c(-0.02, 1.02), ylim = c(0.02, 1.02)) +
  labs(title = "HFD vs SD DEG within each sex, per cell type, by direction",
       subtitle = paste0("Sex DEG: within-sex Wilcoxon p_bonferroni < ", BONF_CUT,
                         " (Bonferroni over the cell type's pooled hits) and |log2FC| > ", LFC_CUT,
                         "; up = higher in HFD.\nOnly pooled hits were tested within sex. A gene cannot be up and ",
                         "down in the same sex, so those regions are 0. One sex only does not imply a ",
                         "sex-specific effect: each sex has half the nuclei.")) +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(size = 7, margin = margin(2, 0, 2, 0)), plot.subtitle = element_text(size = 6.5))
ggsave(file.path(DEG_DIR, "sex_venn_F_vs_M_updown.jpeg"), gp, width = 12, dpi = 300, bg = "white",
       height = 1.2 + 3 * ceiling(uniqueN(vn$annot) / 4), limitsize = FALSE)
print(dcast(cnt, annot ~ key, value.var = "N"))   # key digits: F up, M up, F down, M down
