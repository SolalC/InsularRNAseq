# One palette for the cell types, shared by every figure that colours by
# `annot` -- the composition bars, the annotated UMAP, the reporting UMAPs.
# Sourced, not copied, so a colour changed here changes everywhere at once.
#
# The colours are MetBrewer (Signac, Renoir, Redon, Isfahan1, Lakota, Egypt).
# They were not taken from one palette in order: 13 categories is more than any
# MetBrewer palette carries, and the ones that are long enough spend several
# slots on near-identical pinks. These were picked by farthest-point sampling in
# CIE Lab over the pooled MetBrewer colours, restricted to L 38-78 and chroma
# >= 22 so nothing lands on near-white or near-black, then assigned by lineage:
# neurons warm, macroglia green/blue, vascular and mural muted, Unidentified a
# deliberate neutral grey so a non-category does not read as a cell type.
#
# Worst pair is dE 30.7 (VLCM vs Pericytes), against 29.7 for the Seurat
# polychrome palette this replaced -- so it is no less separable, and the hues
# are organised rather than random.

CELLTYPE_COLS <- c(
  "Excitatory Neurons" = "#bb292c",
  "Inhibitory Neurons" = "#da7901",
  "Astrocytes"         = "#31c7ba",
  "Oligodendrocytes"   = "#376795",
  "Microglia"          = "#c969a1",
  "OPC"                = "#92c051",
  "VLCM"               = "#df8d71",
  "Endothelial Cell"   = "#d8b847",
  "Meningeal Cell"     = "#a1c2ed",
  "Pericytes"          = "#845d29",
  "Smooth Muscle"      = "#633372",
  "Macrophage"         = "#4060c8",
  "Unidentified"       = "#7a7a7a")

# Colours for whatever levels an object actually has, in that order. Anything
# not in the map above -- a renamed label, or a cluster the manual pass never
# named -- is filled from Signac rather than silently dropped, and reported, so
# a new cell type shows up as a colour AND as a line in the log.
celltype_cols <- function(levels, quiet = FALSE) {
  levels <- as.character(levels)
  out <- CELLTYPE_COLS[levels]
  names(out) <- levels
  unknown <- levels[is.na(out)]
  if (length(unknown)) {
    spare <- setdiff(as.character(MetBrewer::met.brewer("Signac", type = "discrete")),
                     CELLTYPE_COLS)
    if (length(unknown) > length(spare))
      spare <- grDevices::rainbow(length(unknown))
    out[unknown] <- spare[seq_along(unknown)]
    if (!quiet)
      message("celltype_cols(): no colour defined for ", paste(unknown, collapse = ", "),
              " -- filled from Signac; add them to lib_celltype_palette.R")
  }
  out
}
