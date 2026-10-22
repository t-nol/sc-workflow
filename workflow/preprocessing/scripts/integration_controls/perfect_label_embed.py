#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import scanpy as sc
import numpy as np
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-i", "--input", required=True, type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
p.add_argument("-f", "--feature_space", required=True, type=str)
p.add_argument("-l", "--layer_key", required=True, type=str)
p.add_argument("-k", "--label_key", required=True, type=str)
p.add_argument("-j", "--jitter", required=False, type=float, default=0)
args = p.parse_args()

# ----------------------------
# Utilities
# ----------------------------
def compute_perfect_embedding(labels, jitter):
    """
    Create a one-hot (perfect) embedding from categorical labels,
    with optional uniform jitter.
    """
    encoded = LabelEncoder().fit_transform(labels)[:, None]
    embedding = OneHotEncoder(sparse_output=False).fit_transform(encoded)
    embedding = embedding.astype(float)

    if jitter > 0:
        embedding += np.random.uniform(-jitter, jitter, embedding.shape)

    return embedding

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":

    # Read inputs
    adata = build_adata(raw_path=args.input)
    adata.X = adata.layers[args.layer_key]

    # Feature space selection
    if args.feature_space == "hvg":
        ad_to_use = adata[:, adata.var["highly_variable"]].copy()
    else:
        ad_to_use = adata.copy()

    # Build perfect latent embedding
    X_perfect = compute_perfect_embedding(
        ad_to_use.obs[args.label_key],
        jitter=args.jitter
    )
    del ad_to_use

    # Save into object
    if args.jitter > 0:
        adata.obsm[f"X_perfect-embed-label-jitter"] = X_perfect
    else:
        adata.obsm[f"X_perfect-embed-label"] = X_perfect
    
    # Write output
    adata.write_zarr(args.output)
