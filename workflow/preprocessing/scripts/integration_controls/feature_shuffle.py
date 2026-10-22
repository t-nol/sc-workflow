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
p.add_argument("-l", "--layer_key", required=True, type=str)
p.add_argument("-s", "--shuffle_by", required=False, type=str, default=None)
args = p.parse_args()

# ----------------------------
# Helpers
# ----------------------------
def shuffle_features(X, partition):
    X_out = X.copy()
    if partition is None:
        partition = np.full(X.shape[0], 0)
    else:
        partition = np.asarray(partition)
    for p in np.unique(partition):
        idx = np.argwhere(partition == p).flatten()
        X_out[idx] = X[np.random.permutation(idx)]
    return X_out

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Read inputs
    adata = build_adata(raw_path=args.input)

    # Shuffle features
    if args.shuffle_by is None:
        partition = None
    else:
        partition = adata.obs[args.shuffle_by]
    
    X = adata.layers["norm_counts"]
    corrected = shuffle_features(X, partition=partition)
    adata.layers[args.layer_key] = corrected

    # Write output
    adata.write_zarr(args.output)
