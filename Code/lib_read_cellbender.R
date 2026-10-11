# Reader for CellBender remove-background output h5.
#
# Seurat's Read10X_h5 CANNOT read these files. It treats every top-level HDF5
# group as a genome and tries to read "<group>/data", but CellBender writes the
# CellRanger-format matrix under `matrix/` alongside its own `droplet_latents`,
# `global_latents` and `metadata` groups, so Read10X_h5 dies with
#   An object with name droplet_latents/data does not exist in this group
# Measured on cellbender 0.4.0 output, 2026-08-31.
#
# This reads the `matrix/` group directly -- standard CellRanger v3 CSC layout
# -- and returns a dgCMatrix with gene names as rownames.

read_cellbender_h5 <- function(path, group = "matrix", use = c("name", "id")) {
  stopifnot(requireNamespace("hdf5r", quietly = TRUE))
  use <- match.arg(use)
  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  if (!(group %in% names(h5))) {
    stop("no '", group, "' group in ", path,
         " -- top-level groups are: ", paste(names(h5), collapse = ", "))
  }
  g <- h5[[group]]

  shape <- g[["shape"]]$read()            # c(n_features, n_barcodes)
  m <- Matrix::sparseMatrix(
    i = g[["indices"]]$read() + 1L,       # CellRanger stores 0-based row indices
    p = g[["indptr"]]$read(),
    x = as.numeric(g[["data"]]$read()),
    dims = as.integer(shape),
    repr = "C")

  feat <- g[["features"]]
  rownames(m) <- make.unique(feat[[use]]$read())
  colnames(m) <- g[["barcodes"]]$read()
  m
}
