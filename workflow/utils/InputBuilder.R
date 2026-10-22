suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(Rarr)
  library(reticulate)
})


# --------------------------------------------------
# Multiple argument for one flag helper
# --------------------------------------------------
get_multi_arg <- function(flag, args) {
    idx <- which(args == flag | startsWith(args, paste0(flag, "=")))
    if (length(idx) == 0) return(NULL)

    # Handle --flag=value case
    if (grepl("=", args[idx])) {
        return(sub("^[^=]+=", "", args[idx]))
    }

    start <- idx + 1
    end <- which(startsWith(args, "-") & seq_along(args) > idx)[1] - 1
    if (is.na(end)) end <- length(args)

    args[start:end]
}


# --------------------------------------------------
# Seurat builder helper 
# --------------------------------------------------
detect_input_type <- function(path) {
  if (grepl("\\.h5ad$", path)) return("h5ad")
  if (dir.exists(path)) return("zarr")
  stop("Unsupported input format: ", path)
}

load_npz_matrix <- function(path, transpose = TRUE) {
    np <- import("numpy", delay_load = TRUE)
    npz <- np$load(path, allow_pickle = TRUE)

    # Extract CSR components
    data    <- as.numeric(npz$f[["data"]])
    indices <- as.integer(npz$f[["indices"]])
    indptr  <- as.integer(npz$f[["indptr"]])
    shape   <- as.integer(npz$f[["shape"]])

    n_rows <- shape[1]
    n_cols <- shape[2]

    # Reconstruct CSR → dgCMatrix
    rows <- rep(seq_len(n_rows), diff(indptr))
    cols <- indices + 1  # 0-based → 1-based

    mat <- Matrix::sparseMatrix(
        i = rows,
        j = cols,
        x = data,
        dims = c(n_rows, n_cols)
    )

    if (transpose) {
        mat <- Matrix::t(mat)
    }

    return(mat)
}

list_zarr_layers <- function(zarr_path) {
  layers_path <- file.path(zarr_path, "layers")
  if (!dir.exists(layers_path)) return(character())
  list.dirs(layers_path, full.names = FALSE, recursive = FALSE)
}

read_zarr_csr <- function(base_path, cells, features) {
  data    <- read_zarr_array(file.path(base_path, "data"))
  indices <- read_zarr_array(file.path(base_path, "indices"))
  indptr  <- read_zarr_array(file.path(base_path, "indptr"))

  rows <- rep(seq_along(cells), diff(indptr))
  cols <- indices + 1

  counts <- Matrix::sparseMatrix(
    i = rows,
    j = cols,
    x = as.numeric(data),
    dims = c(length(cells), length(features))
  )

  counts <- t(counts)
  colnames(counts) <- as.character(cells)
  rownames(counts) <- as.character(features)

  return(counts)
}

load_seurat_from_zarr <- function(zarr_path,
                                  batch_key = NULL,
                                  base_layer = "raw_counts",
                                  norm_layer = "norm_counts") {

  cells    <- read_zarr_array(file.path(zarr_path, "obs/_index"))
  features <- read_zarr_array(file.path(zarr_path, "var/_index"))

  # ---- discover layers ----
  layers <- list_zarr_layers(zarr_path)
  stopifnot(base_layer %in% layers)

  message("→ Found layers: ", paste(layers, collapse = ", "))

  # ---- load base layer ----
  message("   → Loading raw counts: ", base_layer, " as counts.")
  mat <- read_zarr_csr(file.path(zarr_path, "layers", base_layer), cells, features)
  seu <- CreateSeuratObject(
    counts = mat,
    assay  = "RNA"
  )
  rm(mat)
  gc()

  # ---- load norm layer ----
  if (norm_layer %in% layers) {
    message("   → Loading normalized counts: ", norm_layer, " as data.")
    mat <- read_zarr_csr(file.path(zarr_path, "layers", norm_layer), cells, features)
    seu <- SetAssayData(
      seu,
      assay = "RNA",
      layer = "data", 
      new.data = mat
    )
    rm(mat)
    gc()
  }

  # ---- load remaining layers as assays ----
  other_layers <- setdiff(layers, c(base_layer, norm_layer))

  for (layer in other_layers) {
    message("   → Loading layer: ", layer, " as ",layer,".")
    mat <- read_zarr_csr(file.path(zarr_path, "layers", layer), cells, features)
    seu <- SetAssayData(
      seu,
      assay = "RNA",
      layer = layer, 
      new.data = mat
    )
    rm(mat)
    gc()

  }

  # ---- batch key ----
  if (!is.null(batch_key)) {
    cats  <- read_zarr_array(file.path(zarr_path, "obs", batch_key, "categories"))
    codes <- read_zarr_array(file.path(zarr_path, "obs", batch_key, "codes"))
    seu[[batch_key]] <- factor(cats[codes + 1])
  }

  # ---- HVGs ----
  if (file.exists(file.path(zarr_path, "var/highly_variable"))) {
    message("→ Loading highly_variable from var")
    hv <- read_zarr_array(file.path(zarr_path, "var/highly_variable"))
    VariableFeatures(seu) <- rownames(seu)[hv]
  }

  # ---- PCA ----
  if (dir.exists(file.path(zarr_path, "obsm/X_pca"))) {
    message("→ Loading PCA from obsm/X_pca")

    X_pca <- read_zarr_array(file.path(zarr_path, "obsm/X_pca"))
    rownames(X_pca) <- colnames(seu)

    loadings <- read_zarr_array(file.path(zarr_path, "varm/X_pca"))
    rownames(loadings) <- rownames(seu)

    seu[["pca"]] <- CreateDimReducObject(
      embeddings = X_pca,
      loadings   = loadings,
      assay      = DefaultAssay(seu),
      key        = "PC_"
    )
  }

  # ---- PCA variance ----
  if (dir.exists(file.path(zarr_path, "uns/X_pca"))) {
    var <- read_zarr_array(file.path(zarr_path, "uns/X_pca/variance"))
    vr  <- read_zarr_array(file.path(zarr_path, "uns/X_pca/variance_ratio"))
    seu[["pca"]]@misc <- list(
      variance = as.numeric(var),
      variance_ratio = as.numeric(vr)
    )
  }

  return(seu)
}



# --------------------------------------------------
# Seurat builder functions
# --------------------------------------------------
build_seurat <- function(
  raw_path,
  layer = "norm_counts",
  batch_key = NULL
) {
  input_type <- detect_input_type(raw_path)

  if (input_type == "h5ad") {
    message("\n→ No way to load h5ad directly yet, please input as zarr. See convert_h5ad_to_zarr rule.")

  } else {
    message("\n→ Loading ZARR via Rarr")

    seu <- load_seurat_from_zarr(
      zarr_path = raw_path,
      batch_key = batch_key
    )
  }
  return(seu)
}
