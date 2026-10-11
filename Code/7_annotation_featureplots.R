#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
})

PROJECT <- "/scratch/user/uqschauq/Insular"
IN      <- file.path(PROJECT, "Data", "Seurat", "HFDseq_v2_umap.rds")
ANNO    <- file.path(PROJECT, "Results", "annotation_markers", "cluster_annotation.csv")
OUTDIR  <- file.path(PROJECT, "Results", "7_annotation_featureplots")

MIN_GENES  <- 4
PT         <- 0.35
RAST_DPI   <- c(1400, 1400)
COL_LO     <- "grey92"
COL_HI     <- "#08306b"
Q_HI       <- 0.99
CUT_FLOOR  <- 1
SEED       <- 42

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
stopifnot(file.exists(IN), file.exists(ANNO))

SETS <- list(
  "Astrocyte" =
    c("Aqp4", "Gja1", "Slc1a2", "Slc1a3", "Aldh1l1", "Gfap", "Slc7a10", "Grin2c"),
  "Astrocyte Gfap-high" =
    c("Gfap", "Myoc", "Aqp4", "Aldh1l1", "Gja1", "Slc1a2", "Tnc", "C4b"),
  "Oligodendrocyte" =
    c("Plp1", "Mbp", "Mog", "Mal", "St18", "Opalin", "Mag", "Sox10"),
  "OPC" =
    c("Pdgfra", "Cspg4", "Olig1", "Olig2", "Sox10", "Vcan", "Lhfpl3", "Gpr17"),
  "Microglia" =
    c("Cx3cr1", "P2ry12", "Csf1r", "Ctss", "Hexb", "Tmem119", "Selplg", "Trem2"),
  "Perivascular macrophage" =
    c("Mrc1", "Cd163", "F13a1", "Ms4a7", "Dab2", "Stab1", "Lyz2", "Csf1r"),
  "Endothelial" =
    c("Cldn5", "Flt1", "Pecam1", "Slc2a1", "Adgrl4", "Vwf", "Mecom", "Ptprb"),
  "Pericyte" =
    c("Pdgfrb", "Rgs5", "Abcc9", "Kcnj8", "Vtn", "Notch3", "Atp13a5", "Ano1"),
  "Smooth muscle cell" =
    c("Acta2", "Myh11", "Tagln", "Myl9", "Mylk", "Lmod1", "Carmn", "Pln"),
  "VLMC/fibroblast" =
    c("Dcn", "Col1a1", "Slc6a13", "Ranbp3l", "Aldh1a2", "Lama1", "Bmp6", "Slc13a4"),
  "Arachnoid barrier/meningeal" =
    c("Slc47a1", "Prg4", "Cdh1", "Bnc2", "Slc26a7", "Piezo2", "Ahnak2", "Aldh1a2"),
  "L2/3 IT" =
    c("Cux2", "Rasgrf2", "Lamp5", "Calb1", "Trpc6", "Slc17a7", "Satb2", "Snap25"),
  "L2/3 IT Otof" =
    c("Otof", "Cux2", "Rasgrf2", "Ccbe1", "Smoc2", "Cpne6", "Slc17a7", "Satb2"),
  "L4/5 IT" =
    c("Rorb", "Whrn", "Ntng1", "Cpne4", "Thsd7a", "Kcnh5", "Slc17a7", "Cux2"),
  "L5 IT" =
    c("Rspo2", "Cdh9", "C1ql3", "Grm8", "Deptor", "Il1rapl2", "Scube1", "Slc17a7"),
  "L6 IT" =
    c("Cpne4", "Zfp804b", "Pou6f2", "Rmst", "Hgf", "Scn7a", "Slc17a7", "Neurod6"),
  "L5 ET" =
    c("Fezf2", "Bcl11b", "Pou3f1", "Npsr1", "Stac", "Tcerg1l", "Trpc7", "Kcng1"),
  "L5/6 NP" =
    c("Tshz2", "Trhr", "Nxph3", "Etv1", "Vwc2l", "Fezf2", "Htr2c", "Abi3bp"),
  "L6 CT" =
    c("Foxp2", "Syt6", "Tle4", "Bcl11b", "Hs3st4", "Thsd7b", "Zfpm2", "Grik3"),
  "L6b/subplate" =
    c("Ccn2", "Nr4a2", "Drd1", "Clic5", "Tle4", "Prss12", "Pcsk5", "Cdh24"),
  "Claustrum/endopiriform" =
    c("Car3", "Nr4a2", "Gnb4", "Synpr", "Oprk1", "Col11a1", "Ntng2", "Slc17a6"),
  "Pvalb interneuron" =
    c("Pvalb", "Gad1", "Gad2", "Slc32a1", "Kcnmb2", "Kcns3", "Syt2", "Lhx6"),
  "Sst interneuron" =
    c("Sst", "Crhbp", "Elfn1", "Reln", "Grin3a", "Gad1", "Slc32a1", "Lhx6"),
  "Vip interneuron" =
    c("Vip", "Adarb2", "Crh", "Grpr", "Cnr1", "Igf1", "Gad1", "Slc32a1"),
  "Lamp5/Ndnf interneuron" =
    c("Ndnf", "Lamp5", "Npy", "Kit", "Sv2c", "Reln", "Gad1", "Gad2"),
  "Striatal MSN" =
    c("Drd1", "Drd2", "Ppp1r1b", "Rgs9", "Gpr88", "Penk", "Rxrg", "Bcl11b"),
  "Excitatory unresolved" =
    c("Ndst4", "Lmo3", "Ntng1", "Scube1", "Rasgrf2", "Cux2", "Slc17a7", "Snap25"),
  "Low-complexity/mixed" =
    c("Snap25", "Slc17a7", "Aqp4", "Plp1", "Cx3cr1", "Csf1r", "Flt1", "Pdgfra"),
  "Oligodendrocyte-microglia doublet" =
    c("Plp1", "Mbp", "St18", "Mog", "Cx3cr1", "P2ry12", "Csf1r", "Hexb"),
  "Pan-neuronal reference" =
    c("Snap25", "Syt1", "Rbfox3", "Slc17a7", "Satb2", "Gad1", "Gad2", "Slc32a1")
)

obj <- readRDS(IN)
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "seurat_clusters"
msg("loaded ", ncol(obj), " nuclei | ", nlevels(obj$seurat_clusters), " clusters")

anno <- fread(ANNO, colClasses = list(character = "cluster"))
stopifnot(!anyDuplicated(anno$cluster),
          setequal(anno$cluster, levels(obj$seurat_clusters)))
m <- match(as.character(obj$seurat_clusters), anno$cluster)
obj$cluster_label <- anno$label[m]
obj$cell_class    <- anno$cell_class[m]

anno[, cluster_int := as.integer(cluster)]
setorder(anno, cluster_int)
byLabel <- anno[, .(clusters = paste(cluster, collapse = ", "),
                    n_clusters = .N,
                    conf = paste(unique(confidence), collapse = "/")),
                by = label]

missing_sets <- setdiff(byLabel$label, names(SETS))
if (length(missing_sets))
  stop("no marker set defined for: ", paste(missing_sets, collapse = " | "))

present <- rownames(obj[["RNA"]])
for (nm in names(SETS)) {
  keep <- SETS[[nm]][SETS[[nm]] %in% present]
  drop <- setdiff(SETS[[nm]], present)
  if (length(drop))
    msg("WARNING: '", nm, "' dropping absent gene(s): ", paste(drop, collapse = ", "))
  if (length(keep) < MIN_GENES)
    stop("'", nm, "' left with ", length(keep), " gene(s), below MIN_GENES=", MIN_GENES)
  SETS[[nm]] <- keep
}

sizes <- as.data.table(obj@meta.data)[, .N, by = .(cluster_label, cell_class)]
byLabel <- merge(byLabel, sizes, by.x = "label", by.y = "cluster_label", all.x = TRUE)
setorder(byLabel, -N)
fwrite(data.table(label = names(SETS),
                  genes = vapply(SETS, paste, "", collapse = " "))[
         byLabel, on = "label"],
       file.path(OUTDIR, "featureplot_gene_sets.csv"))
print(byLabel)

slug <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))

expr <- GetAssayData(obj, assay = "RNA", layer = "data")
cutoffs <- function(g) {
  v <- vapply(g, function(x) {
    e <- expr[x, ]
    e <- e[e > 0]
    if (!length(e)) return(CUT_FLOOR)
    max(as.numeric(quantile(e, Q_HI)), CUT_FLOOR)
  }, numeric(1))
  fwrite(data.table(gene = g, max_cutoff = round(v, 3)),
         file.path(OUTDIR, "featureplot_scale_cutoffs.csv"), append = TRUE)
  unname(v)
}

ref <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
               label = TRUE, label.size = 4, repel = TRUE,
               raster = TRUE, raster.dpi = RAST_DPI, pt.size = PT) +
  coord_equal() + NoLegend() +
  labs(title = sprintf("%d clusters | %s nuclei", nlevels(obj$seurat_clusters),
                       format(ncol(obj), big.mark = ",")))
ggsave(file.path(OUTDIR, "umap_reference_clusters.png"), ref,
       width = 11, height = 10, dpi = 200)

lab <- DimPlot(obj, reduction = "umap", group.by = "cluster_label",
               label = TRUE, label.size = 3, repel = TRUE,
               raster = TRUE, raster.dpi = RAST_DPI, pt.size = PT) +
  coord_equal() + NoLegend() + labs(title = "Proposed labels")
ggsave(file.path(OUTDIR, "umap_reference_labels.png"), lab,
       width = 11, height = 10, dpi = 200)

cls <- DimPlot(obj, reduction = "umap", group.by = "cell_class",
               label = TRUE, label.size = 3.4, repel = TRUE,
               raster = TRUE, raster.dpi = RAST_DPI, pt.size = PT) +
  coord_equal() + labs(title = "Proposed classes")
ggsave(file.path(OUTDIR, "umap_reference_classes.png"), cls,
       width = 13, height = 10, dpi = 200)
msg("reference UMAPs written")

panel_theme <- theme(plot.title = element_text(size = 13, face = "bold"),
                     axis.title = element_blank(),
                     axis.text  = element_text(size = 6),
                     legend.key.width = grid::unit(0.25, "cm"))

for (i in seq_len(nrow(byLabel))) {
  nm   <- byLabel$label[i]
  cl   <- strsplit(byLabel$clusters[i], ", ")[[1]]
  gs   <- SETS[[nm]]
  cells <- colnames(obj)[as.character(obj$seurat_clusters) %in% cl]

  hl <- DimPlot(obj, reduction = "umap", raster = TRUE, raster.dpi = RAST_DPI,
                pt.size = PT, cells.highlight = list(x = cells),
                cols.highlight = "#d1341f", sizes.highlight = PT) +
    coord_equal() + NoLegend() +
    labs(title = sprintf("cluster %s  (n = %s)", byLabel$clusters[i],
                         format(byLabel$N[i], big.mark = ","))) +
    panel_theme

  fp <- FeaturePlot(obj, features = gs, reduction = "umap", combine = FALSE,
                    order = TRUE, raster = TRUE, raster.dpi = RAST_DPI,
                    pt.size = PT, cols = c(COL_LO, COL_HI),
                    max.cutoff = cutoffs(gs))
  fp <- lapply(fp, function(p) p + coord_equal() + panel_theme)

  p <- wrap_plots(c(list(hl), fp), ncol = 3) +
    plot_annotation(
      title = nm,
      subtitle = sprintf("proposed label for cluster %s | confidence %s | %s markers",
                         byLabel$clusters[i], byLabel$conf[i], length(gs)),
      theme = theme(plot.title = element_text(size = 20, face = "bold"),
                    plot.subtitle = element_text(size = 12, colour = "grey35")))

  nr <- ceiling((length(gs) + 1) / 3)
  f  <- file.path(OUTDIR, sprintf("featureplot_%02d_%s.png", i, slug(nm)))
  ggsave(f, p, width = 16, height = 5.2 * nr + 0.8, dpi = 150, limitsize = FALSE)
  msg("[", i, "/", nrow(byLabel), "] ", basename(f), " | ", paste(gs, collapse = " "))
}

gs <- SETS[["Pan-neuronal reference"]]
fp <- FeaturePlot(obj, features = gs, reduction = "umap", combine = FALSE,
                  order = TRUE, raster = TRUE, raster.dpi = RAST_DPI,
                  pt.size = PT, cols = c(COL_LO, COL_HI),
                  max.cutoff = cutoffs(gs))
fp <- lapply(fp, function(p) p + coord_equal() + panel_theme)
p <- wrap_plots(fp, ncol = 3) +
  plot_annotation(title = "Pan-neuronal reference",
                  subtitle = "not a cluster assignment -- orientation panel for the neuronal half of the UMAP",
                  theme = theme(plot.title = element_text(size = 20, face = "bold"),
                                plot.subtitle = element_text(size = 12, colour = "grey35")))
ggsave(file.path(OUTDIR, "featureplot_00_pan_neuronal_reference.png"), p,
       width = 16, height = 5.2 * ceiling(length(gs) / 3) + 0.8, dpi = 150,
       limitsize = FALSE)
msg("pan-neuronal reference written | ", paste(gs, collapse = " "))

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
msg("done: ", OUTDIR)


  library(Nebulosa)
  library(cowplot)
  p0 <- UMAPPlot(obj, label = T, raster = F) + theme(legend.position= 'none')
  # Cluster of interest: 13
  # Clusters:
  cl <- as.character(obj$seurat_clusters)
  # Excitatory neurons:
  p1 <- plot_density(obj, 'Slc17a7')
  p2 <- plot_density(obj, 'Slc17a6')
  p3 <- plot_density(obj, 'Slc17a8')
  p4 <- plot_density(obj, 'Camk2a')
  p5 <- plot_density(obj, 'Grin2b')
  pExcitatory <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pExcitatory, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/ExcitatoryNeurons.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(2,3,4,32,21,26,8,16,19,6,27,30,5, 17, 12, 41, 13))] <- 'Excitatory Neurons'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Inhibitory neurons:
  p1 <- plot_density(obj, 'Gad1')
  p2 <- plot_density(obj, 'Gad2')
  p3 <- plot_density(obj, 'Slc32a1')
  p4 <- plot_density(obj, 'Pvalb')
  p5 <- plot_density(obj, 'Sst')
  pInhibitory <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pInhibitory, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/InhibitoryNeurons.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(9,11,35,20,14, 25))] <- 'Inhibitory Neurons'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Microglia
  p1 <- plot_density(obj, 'Aif1')
  p2 <- plot_density(obj, 'Itgam')
  p3 <- plot_density(obj, 'Ptprc')
  p4 <- plot_density(obj, 'Cx3cr1')
  p5 <- plot_density(obj, 'Tmem119')
  pMicroglia <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pMicroglia, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/Microglia.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(7, 29))] <- 'Microglia'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Astrocytes
  p1 <- plot_density(obj, 'Gfap')
  p2 <- plot_density(obj, 'Aldh1l1')
  p3 <- plot_density(obj, 'S100b')
  p4 <- plot_density(obj, 'Sox9')
  p5 <- plot_density(obj, 'Slc1a2')
  pAstrocytes <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pAstrocytes, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/Astrocytes.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(0,23))] <- 'Astrocytes'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Oligodendrocytes:
  p1 <- plot_density(obj, 'Olig2')
  p2 <- plot_density(obj, 'Sox10')
  p3 <- plot_density(obj, 'Mbp')
  p4 <- plot_density(obj, 'Plp1')
  p5 <- plot_density(obj, 'Myrf')
  pOligo <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pOligo, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/Olygodendrocytes.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(1))] <- 'Oligodendrocytes'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # OPC:
  p1 <- plot_density(obj, 'Pdgfra')
  p2 <- plot_density(obj, 'Olig2')
  p3 <- plot_density(obj, 'Cspg4')
  p4 <- plot_density(obj, 'Sox10')
  pOPC <- plot_grid(p0,p1,p2,p3,p4, nrow = 3)
  ggsave(pOPC, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/OPC.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(10, 40, 36))] <- 'OPC'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # VLMC:
  p1 <- plot_density(obj, 'Dcn')
  p2 <- plot_density(obj, 'Bmp6')
  p3 <- plot_density(obj, 'Cped1')
  p4 <- plot_density(obj, 'Cfh')
  p5 <- plot_density(obj, 'Slc6a13')
  pVLCM <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pVLCM, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/VLCM.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(15, 38))] <- 'VLCM'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Endothelial cells:
  p1 <- plot_density(obj, 'Flt1')
  p2 <- plot_density(obj, 'Ly6c1')
  p3 <- plot_density(obj, 'Pecam1')
  p4 <- plot_density(obj, 'Cdh5')
  p5 <- plot_density(obj, 'Cldn5')
  pEndo <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pEndo, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/Endothelial.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(18, 31, 37))] <- 'Endothelial Cell'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Meningeal cells:
  p1 <- plot_density(obj, 'Bnc2')
  p2 <- plot_density(obj, 'Slc47a1')
  p3 <- plot_density(obj, 'Slc6a13')
  pMeninge <- plot_grid(p0,p1,p2,p3, nrow = 2)
  ggsave(pMeninge, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/MeningealCells.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(22))] <- 'Meningeal Cell'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Pericytes:
  p1 <- plot_density(obj, 'Kcnj8')
  p2 <- plot_density(obj, 'Vtn')
  p3 <- plot_density(obj, 'Abcc9')
  p4 <- plot_density(obj, 'Higd1b')
  p5 <- plot_density(obj, 'Cd248')
  pPericytes <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pPericytes, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/Pericytes.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(24, 39))] <- 'Pericytes'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Unidentified
  # Assign clusters:
  cl[which(cl %in% c(28))] <- 'Unidentified'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Smooth muscle cell:
  p1 <- plot_density(obj, 'Acta2')
  p2 <- plot_density(obj, 'Myh11')
  p3 <- plot_density(obj, 'Tagln')
  p4 <- plot_density(obj, 'Myocd')
  p5 <- plot_density(obj, 'Lmod1')
  pSmoothMuscle <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pSmoothMuscle, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/SmoothMuscleCell.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(33))] <- 'Smooth Muscle'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T)
  # Macrophage:
  p1 <- plot_density(obj, 'Mrc1')
  p2 <- plot_density(obj, 'Pf4')
  p3 <- plot_density(obj, 'Lyz2')
  p4 <- plot_density(obj, 'Ms4a7')
  p4 <- plot_density(obj, 'F13a1')
  p5 <- plot_density(obj, 'Lyve1')
  pMacrophage <- plot_grid(p0,p1,p2,p3,p4,p5, nrow = 3)
  ggsave(pMacrophage, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/Macrophage.png', 
        width = 14, height = 14, dpi = 300)
  # Assign clusters:
  cl[which(cl %in% c(34))] <- 'Macrophage'
  obj$annot <- cl
  UMAPPlot(obj, group.by = 'annot', label = T, raster= F)
  # Write the annotated UMAP:
  source('/scratch/user/uqschauq/Insular/Code/lib_celltype_palette.R')
  (x <- UMAPPlot(obj, group.by = 'annot', label = T, raster = F,
                 cols = celltype_cols(levels(factor(obj$annot)))))
  ggsave(x, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/AnnotatedUmap.png',
          width = 14, height = 14, dpi = 300)
  # Write the dotplot with all the genes
markers <- list(
  'Excitatory Neurons' = c('Slc17a7','Slc17a6','Slc17a8','Camk2a','Grin2b', 'Tubb3', 'Eno2', 'Syp'),
  'Inhibitory Neurons' = c('Gad1','Gad2','Slc32a1','Pvalb','Sst'),
  'Microglia' = c('Aif1','Itgam','Ptprc','Cx3cr1','Tmem119'),
  'Astrocytes' = c('Gfap','Aldh1l1','S100b','Sox9','Slc1a2'),
  'Oligodendrocytes' = c('Olig2','Sox10','Mbp','Plp1','Myrf'),
  'OPC' = c('Pdgfra','Cspg4','Olig2','Sox10'),
  'VLCM' = c('Dcn','Bmp6','Cped1','Cfh','Slc6a13'),
  'Endothelial Cell' = c('Flt1','Ly6c1','Pecam1','Cdh5','Cldn5'),
  'Meningeal Cell' = c('Bnc2','Slc47a1','Slc6a13'),
  'Pericytes' = c('Kcnj8','Vtn','Abcc9','Higd1b','Cd248'),
  'Smooth Muscle' = c('Acta2','Myh11','Tagln','Myocd','Lmod1'),
  'Macrophage' = c('Mrc1','Pf4','Lyz2','Ms4a7','F13a1','Lyve1'))
markers <- setNames(lapply(seq_along(markers), function(i)
  setdiff(markers[[i]], unlist(markers[seq_len(i - 1)]))), names(markers))
genes <- unlist(markers, use.names = F)
celltypes <- names(markers)
nmk <- lengths(markers)
bands <- data.frame(xmin = cumsum(c(0, head(nmk, -1))) + 0.5, xmax = cumsum(nmk) + 0.5)
obj$annot <- factor(obj$annot, levels = rev(c(celltypes, setdiff(unique(obj$annot), celltypes))))
Idents(obj) <- obj$annot
p <- DotPlot(obj, features = genes, group.by = 'annot')
p$layers <- c(annotate('rect', xmin = bands$xmin, xmax = bands$xmax, ymin = -Inf, ymax = Inf,
                       fill = rep(c('grey90','white'), length.out = nrow(bands)), alpha = 0.5),
              p$layers)
(pDot <- p +
  annotate('text', x = (bands$xmin + bands$xmax) / 2, y = Inf, label = celltypes,
           vjust = rep(c(-2.4, -0.7), length.out = length(celltypes)), size = 3) +
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.margin = margin(t = 48, r = 6, b = 6, l = 6)))
ggsave(pDot, filename = '/scratch/user/uqschauq/Insular/Results/7_annotation_featureplots/manual/AnnotatedDotPlot.jpeg',
        width = 16, height = 8, dpi = 300)
saveRDS(obj, '/scratch/user/uqschauq/Insular/Data/Seurat/HFDseq_v2_annotated.rds', compress = FALSE)
