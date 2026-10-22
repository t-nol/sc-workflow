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
p.add_argument("-m", "--methods", required=True, nargs="+", type=str)
p.add_argument("-d", "--do_3d", required=True, type=bool)
p.add_argument("-k", "--neighbors_key", required=False, type=str, default=None)
args = p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     

    # Read inputs
    adata = build_adata(raw_path=args.input)

    # UMAP embeddings
    if "umap" in args.methods:
        # 2D UMAP
        sc.tl.umap(adata, n_components=2, key_added="X_umap_2D")
        # 3D UMAP
        if args.do_3d:
            sc.tl.umap(adata, n_components=3, key_added="X_umap_3D")
    
    # UMAP embeddings
    if "diffmap" in args.methods:
        sc.tl.diffmap(adata, n_comps=15, neighbors_key=args.neighbors_key) 
        #resulting embedding in adata.obsm['X_diffmap']

    # Write output
    adata.write_zarr(args.output)
