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

# Direction split: one Venn for genes up in HFD (F up vs M up), one for genes down.
# F up and F down are disjoint, so a single 4-set Venn would be mostly impossible regions.
# A gene up in one sex and down in the other is in "F only" of one Venn and "M only" of the
# other; their number is given in the strip.
vn[, `:=`(F_dir = fifelse(F_deg, fifelse(F_avg_log2FC > 0, "up", "down"), NA_character_),
          M_dir = fifelse(M_deg, fifelse(M_avg_log2FC > 0, "up", "down"), NA_character_))]
cnt <- rbindlist(lapply(c("up", "down"), function(dr) vn[, {
  f <- !is.na(F_dir) & F_dir == dr; m <- !is.na(M_dir) & M_dir == dr
  data.table(dir = dr, x = c(-1.05, 0, 1.05), n = c(sum(f & !m), sum(f & m), sum(m & !f)),
             n_opp = sum(F_deg & M_deg & F_dir != M_dir))
}, by = annot]))
cnt[, facet := factor(sprintf("%s: %s in HFD\n(%d opposite sign in F vs M)", gsub("_", " ", annot),
                              dir, n_opp))]
cnt[, facet := factor(facet, levels = unique(cnt[order(annot, -xtfrm(dir)), facet]))]  # up, then down

th <- seq(0, 2 * pi, length.out = 200)
circ <- rbind(data.table(sex = "F", x = -0.6 + cos(th), y = sin(th)),
              data.table(sex = "M", x =  0.6 + cos(th), y = sin(th)))
gp <- ggplot() +
  geom_polygon(data = circ, aes(x, y, fill = sex, group = sex), alpha = 0.25, colour = "grey30", linewidth = 0.3) +
  geom_text(data = data.table(x = c(-1.2, 1.2), y = 1.2, lab = c("Female", "Male")),
            aes(x, y, label = lab), size = 3, fontface = "bold") +
  geom_text(data = cnt, aes(x, 0, label = n), size = 4.5) +
  facet_wrap(~ facet, ncol = 4) +
  scale_fill_manual(values = c(F = "#d95f02", M = "#1b9e77"), guide = "none") +
  coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.1, 1.4)) +
  labs(title = "HFD vs SD DEG within each sex, per cell type, split by direction",
       subtitle = paste0("Sex DEG: within-sex Wilcoxon p_bonferroni < ", BONF_CUT,
                         " (Bonferroni over the cell type's pooled hits) and |log2FC| > ", LFC_CUT,
                         "; up = log2FC > 0 (higher in HFD).\nOnly pooled hits were tested within sex. ",
                         "One circle only does not imply a sex-specific effect: each sex has half the nuclei.")) +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(size = 7, margin = margin(2, 0, 2, 0)), plot.subtitle = element_text(size = 6.5))
ggsave(file.path(DEG_DIR, "sex_venn_F_vs_M_updown.jpeg"), gp, width = 12, dpi = 300, bg = "white",
       height = 1.2 + 2.6 * ceiling(2 * uniqueN(vn$annot) / 4), limitsize = FALSE)
print(dcast(cnt[, .(annot, dir, region = c("F only", "both", "M only")[match(x, c(-1.05, 0, 1.05))], n)],
            annot + dir ~ region, value.var = "n"))
