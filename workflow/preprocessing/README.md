# PREPROCESSING & INTEGRATION MODULE

## Goal

Preprocessing on a single-cell AnnData object, optionally batch corrected by batch (as defined by `batch_key` in config). 

## Steps

- **Normalization** : Normalizes gene expression per cell (e.g., library-size scaling and log-transformation).

- **Feature space** : Either takes all features in dataset or identifies highly variable genes (hvg) per batch or across the full dataset based on the configuration. These genes are flagged in `var/highly_variable` and used in downstream dimensionality reduction and integration.

- **Embedding space** : Performs PCA after scaling data. Optional regression of covariates (cell_cycle) can be applied to reduce unwanted variation if specified in config in `pca/regress_cell_cycle`. 

- **Graph space** : Neighbors graph built by knn or integration method. 

- **Projection** : Performs UMAP (2D and/or 3D) as well as diffmap dimension reduction methods.

- **Integration** : Batch effect correction happens at different steps of the preprocessing workflow depending on the method used.

- **Report** : Generates an interactive html report for each feature_space + integration method.

## Config parameters
