#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy as sc
import numpy as np

project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-i", "--input", required=True, type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
p.add_argument("-n", "--knn_neighbors", required=True, type=int)
p.add_argument("-e", "--embed_key", required=False, type=str, default="X_pca")
args = p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Read inputs
    adata = build_adata(raw_path=args.input)
    print(adata)
    
    # KNN graph computation
    sc.pp.neighbors(adata, n_neighbors=args.knn_neighbors, use_rep=args.embed_key, key_added=None)

    # Write outputs
    adata.write_zarr(args.output)
    