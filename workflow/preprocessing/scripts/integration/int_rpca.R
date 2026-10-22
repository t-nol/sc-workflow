suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(future)
})

# Increase allowed memory for parallel exports
options(future.globals.maxSize = 8 * 1024^3) #here 8Gb
plan("multisession", workers = 1)

# Load prepro.R functions
script_path <- normalizePath(sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]))
project_root <- dirname(dirname(dirname(dirname(script_path))))
source(file.path(project_root, "utils", "InputBuilder.R"))

# ----------------------------
# Argument parsing
# ----------------------------
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Input path"),
    make_option(c("-o", "--output"), type = "character", help = "Output path"),
    make_option(c("-n", "--nfeatures"), type = "integer", help = "Number of highly variable features to identidy"),
    make_option(c("-f", "--feature_space"), type = "character", help = "Wether to use full or hvg feature subset for integration"),
    make_option(c("-b", "--batch_key"), type = "character", help = "Batch key (cell metadata column)")
)
parser <- OptionParser(option_list = option_list)
raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(parser, args = raw_args)

# ----------------------------
# Load input as seurat object
# ----------------------------
obj <- build_seurat(
    raw_path = args$input, 
    batch_key = args$batch_key
    )

# ----------------------------
# RPCA integration
# ----------------------------
## RPCA anc CCA integration are seurat specific methods that require
## using Seurat's own preprocessing workflow.
# Preprocessing per batch
obj[["RNA"]] <- split(obj[["RNA"]], f = obj@meta.data[[args$batch_key]])

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = args$nfeatures)

if (args$feature_space == "full") {
    features_to_use <- rownames(obj)
} else if (args$feature_space == "hvg") {
    features_to_use <- Seurat::VariableFeatures(obj)
} else {
    stop("feature_subset must be 'full' or 'hvg'")
}

obj <- ScaleData(obj, features = features_to_use, verbose = FALSE)
obj <- RunPCA(obj, features = features_to_use, verbose = FALSE)

# Integration
obj <- IntegrateLayers(
    object = obj, 
    features = features_to_use,
    method = RPCAIntegration,
    orig.reduction = "pca", 
    new.reduction = "rpca",
    verbose = FALSE
    )

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

embedding <- obj[["rpca"]]@cell.embeddings

# ----------------------------
# Write output
# ----------------------------
np <- reticulate::import("numpy", delay_load = TRUE)

np$savez_compressed(
    args$output,
    X_corrected = embedding
)
