suppressPackageStartupMessages({
    library(optparse)
    library(Seurat)
    library(simspec)
    library(Rarr)
})

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
    make_option(c("-f", "--feature_space"), type = "character", help = "Wether to use full or hvg feature subset for integration"),
    make_option(c("-b", "--batch_key"), type = "character", help = "Batch key (cell metadata column)"),
    make_option(c("-c", "--correlation"), type = "character", default = "spearman",help = "CSS correlation method")
)
parser <- OptionParser(option_list = option_list)
raw_args <- commandArgs(trailingOnly = TRUE)
args <- parse_args(parser, args = raw_args)

# Enforce required arguments
if (is.null(args$input) || is.null(args$output) || is.null(args$batch)) {
    print_help(parser)
    stop("Missing required arguments", call. = FALSE)
}

# define correlation method and spectrum type
if(args$correlation == "pearson"){
    spectrum <- "corr_kernel"
}else if(args$correlation == "spearman"){
    spectrum <- "corr_ztransform"
}else{
    stop("Correlation method must be pearson or spearman.", call. = FALSE)
}

# ----------------------------
# Load input as seurat object
# ----------------------------
obj <- build_seurat(
    raw_path = args$input, 
    batch_key = args$batch_key
    )

if (args$feature_space == "full") {
    var_genes <- rownames(obj)
} else if (args$feature_space == "hvg") {
    var_genes <- Seurat::VariableFeatures(obj)
} else {
    stop("feature_space must be 'full' or 'hvg'")
}

# ----------------------------
# CSS integration
# ----------------------------
obj <- cluster_sim_spectrum(
    object = obj, 
    label_tag = args$batch_key,
    corr_method = args$correlation,
    spectrum_type = spectrum,
    use_dr = "pca",
    dims_use = seq_len(ncol(obj[["pca"]]@cell.embeddings)),
    var_genes = var_genes
    )

embedding <- obj[["css"]]@cell.embeddings

# ----------------------------
# Write output 
# ----------------------------
np <- reticulate::import("numpy", delay_load = TRUE)

np$savez_compressed(
    args$output,
    X_corrected = embedding
)
