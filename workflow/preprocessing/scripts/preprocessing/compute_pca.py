#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import numpy as np
import anndata as ad
import scanpy as sc

project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-i", "--input", required=True, type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
p.add_argument("-n", "--n_comps", required=True, type=int)
p.add_argument("-l", "--layer_key", required=True, type=str)
p.add_argument("-f", "--feature_space", required=True, type=str)
p.add_argument("-r", "--regress_cc", required=False, type=bool)
p.add_argument("-b", "--batch_key", required=False, type=str, default=None)
args = p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Read inputs
    adata = build_adata(raw_path=args.input)
    adata.X = adata.layers[args.layer_key].copy()

    # Cell Cycle impact regression
    if args.regress_cc:
        if "S_score" in adata.obs.columns and "G2M_score" in adata.obs.columns:
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    
    # Scaling
    sc.pp.scale(adata)

    # PCA
    if args.batch_key is None:
        # Standard PCA
        sc.pp.pca(
            adata,
            n_comps=args.n_comps,
            mask_var=None if args.feature_space == "full" else "highly_variable",
            svd_solver="covariance_eigh",
            key_added="X_pca",
        )
    
    else:
        # Batch-wise PCA
        X_pca = np.zeros((adata.n_obs, args.n_comps), dtype=float)

        for batch in adata.obs[args.batch_key].unique():
            batch_idx = adata.obs[args.batch_key] == batch
            n_cells = int(batch_idx.sum())

            if n_cells <= 1:
                continue

            n_batch_comps = min(args.n_comps, n_cells - 1)
            solver = "arpack" if n_cells > args.n_comps else "full"

            ad_batch = adata[batch_idx].copy()

            sc.pp.pca(
                ad_batch,
                n_comps=n_batch_comps,
                mask_var=None if args.feature_space == "full" else "highly_variable",
                svd_solver=solver,
                key_added="X_batchwise-pca",
            )

            X_pca[batch_idx, :n_batch_comps] = ad_batch.obsm["X_batchwise-pca"]
        
        adata.obsm["X_batchwise-pca"] = X_pca

    # Write output
    adata.write_zarr(args.output)
