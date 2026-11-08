# Venn of female vs male HFD-vs-SD DEG per cell type, from DEG_percell_sexsplit_<ct>.csv.
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
fwrite(vn[order(annot, region, gene)], file.path(DEG_DIR, "sex_venn_membership.csv"))

updown <- function(x) sprintf("%d up / %d down", sum(x > 0), sum(x < 0))
cnt <- vn[, {
  fo <- region == "F only"; mo <- region == "M only"; bo <- region == "both"
  opp <- bo & sign(F_avg_log2FC) != sign(M_avg_log2FC)
  data.table(x = c(-1.05, 0, 1.05), n = c(sum(fo), sum(bo), sum(mo)),
             sub = c(updown(F_avg_log2FC[fo]),
                     paste0(updown(F_avg_log2FC[bo & !opp]), if (any(opp)) sprintf("\n%d opposite sign", sum(opp)) else ""),
                     updown(M_avg_log2FC[mo])),
             n_hit = .N, n_neither = sum(region == "neither"))
}, by = annot]
cnt[, facet := sprintf("%s\n%d pooled hits; %d DEG in neither sex", gsub("_", " ", annot), n_hit, n_neither)]

th <- seq(0, 2 * pi, length.out = 200)
circ <- rbind(data.table(sex = "F", x = -0.6 + cos(th), y = sin(th)),
              data.table(sex = "M", x =  0.6 + cos(th), y = sin(th)))
gp <- ggplot() +
  geom_polygon(data = circ, aes(x, y, fill = sex, group = sex), alpha = 0.25, colour = "grey30", linewidth = 0.3) +
  geom_text(data = data.table(x = c(-1.2, 1.2), y = 1.2, lab = c("Female", "Male")),
            aes(x, y, label = lab), size = 3, fontface = "bold") +
  geom_text(data = cnt, aes(x, 0.12, label = n), size = 4.5) +
  geom_text(data = cnt, aes(x, -0.25, label = sub), size = 1.9, lineheight = 0.9) +
  facet_wrap(~ facet, ncol = 4) +
  scale_fill_manual(values = c(F = "#d95f02", M = "#1b9e77"), guide = "none") +
  coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.1, 1.4)) +
  labs(title = "HFD vs SD DEG within each sex, per cell type",
       subtitle = paste0("Sex DEG: within-sex Wilcoxon p_bonferroni < ", BONF_CUT,
                         " (Bonferroni over the cell type's pooled hits) and |log2FC| > ", LFC_CUT,
                         ".\nOnly pooled hits were tested within sex. One circle only does not imply a ",
                         "sex-specific effect: each sex has half the nuclei.")) +
  theme_void(base_size = 8) +
  theme(strip.text = element_text(size = 7, margin = margin(2, 0, 2, 0)), plot.subtitle = element_text(size = 6.5))
ggsave(file.path(DEG_DIR, "sex_venn_F_vs_M.pdf"), gp, width = 12,
       height = 1.2 + 2.6 * ceiling(uniqueN(vn$annot) / 4), limitsize = FALSE)
print(dcast(vn, annot ~ region, fun.aggregate = length))
