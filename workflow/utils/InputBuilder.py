from pathlib import Path
import re
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix, load_npz


# ----------------------------
# Adata builder functions
# ----------------------------
def read_adata(path: Path) -> ad.AnnData:
    if path.suffix == ".zarr":
        return ad.read_zarr(path)
    elif path.suffix == ".h5ad":
        return sc.read_h5ad(path)
    else:
        raise ValueError(f"Unsupported input type: {path}")


def build_adata(
    raw_path: Path,
    norm_path: Path | None = None,
    features_path: Path | None = None,
    pca_path: Path | None = None,
    neighbors_path: Path | None = None,
) -> ad.AnnData:
    """
    Incrementally rebuild AnnData depending on which inputs are provided.
    Was used in previous workflow to save disk space for large datasets.
    """

    # ---- base object ----
    adata = read_adata(raw_path)

    if "raw_counts" not in adata.layers:
        adata.layers["raw_counts"] = adata.X.copy()

    # ---- normalized counts ----
    if norm_path is not None:
        adata.layers["norm_counts"] = load_npz(norm_path)

    # ---- feature annotations ----
    if features_path is not None:
        features_df = pd.read_csv(features_path, index_col=0)
        adata.var = adata.var.join(features_df)

    # ---- PCA ----
    if pca_path is not None:
        pca_data = np.load(pca_path)
        adata.obsm["X_pca"] = pca_data["X_pca"]

        if "loadings" in pca_data.files:
            adata.varm["PCs"] = pca_data["loadings"]

        if {"variance", "variance_ratio"} <= set(pca_data.files):
            adata.uns["pca"] = {
                "variance": pca_data["variance"],
                "variance_ratio": pca_data["variance_ratio"],
            }

    # ---- neighbors ----
    if neighbors_path is not None:
        neighbors_data = np.load(neighbors_path, allow_pickle=True)
        n = adata.n_obs

        adata.obsp["connectivities"] = csr_matrix(
            (
                neighbors_data["connectivities"],
                neighbors_data["connectivities_indices"],
                neighbors_data["connectivities_indptr"],
            ),
            shape=(n, n),
        )

        adata.obsp["distances"] = csr_matrix(
            (
                neighbors_data["distances"],
                neighbors_data["distances_indices"],
                neighbors_data["distances_indptr"],
            ),
            shape=(n, n),
        )

        adata.uns["neighbors"] = neighbors_data["neighbors_params"].item()

    return adata


def adata_join_without_overlaps(adata, df, join_on='obs'):
    """
    Join a DataFrame to an AnnData object without overwriting existing columns.

    Parameters
    ----------
    adata : AnnData
        AnnData object to modify.
    df : pd.DataFrame
        DataFrame to join to adata.obs or adata.var.
    join_on : str, default 'obs'
        Whether to join on 'obs' (cells) or 'var' (genes).

    Returns
    -------
    adata : AnnData
        AnnData object with the DataFrame joined, excluding overlapping columns.
    """
    existing_columns = getattr(adata, join_on).columns
    df_filtered = df.loc[:, ~df.columns.isin(existing_columns)]

    if join_on == 'obs':
        adata.obs = adata.obs.join(df_filtered)
    elif join_on == 'var':
        adata.var = adata.var.join(df_filtered)
    
    return adata


def validate_adata(adata, split_by):
    """Check that the AnnData object has all required fields."""

    # Check var columns
    required_var_cols = ["gene_id", "gene_name"]
    missing_cols = [c for c in required_var_cols if c not in adata.var.columns]
    if missing_cols:
        raise ValueError(f"Input AnnData is missing required columns in .var: {missing_cols}")
    
    # Required obs columns
    required_obs_cols = []
    if split_by is not None:
        required_obs_cols.append(split_by)
    missing_cols = [c for c in required_obs_cols if c not in adata.obs.columns]
    if missing_cols:
        raise ValueError(f"Input AnnData is missing required columns in .obs: {missing_cols}")
    
    print("Validation passed: all required fields are present.")


def sanitize_batch_name(batch):
    """
    Convert a batch name to a safe string for filenames:
    - Replace spaces and any non-alphanumeric/dash/underscore characters with underscore
    - Keeps letters, numbers, dash (-), and underscore (_)
    """
    batch_str = str(batch)  # in case batch is numeric or other type
    safe_batch = re.sub(r"[^A-Za-z0-9_-]", "_", batch_str)
    return safe_batch

