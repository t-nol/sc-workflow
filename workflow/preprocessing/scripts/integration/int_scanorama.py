#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy.external as sce

project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True, type=Path)
    p.add_argument("-o", "--output", required=True, type=Path)
    p.add_argument("-b", "--batch_key", required=True, type=str)
    p.add_argument("-e", "--embed_key", required=False, type=str, default="X_pca")
    p.add_argument("-n", "--n_neighbors", required=True, type=int)
    return p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    # Scanorama
    print(f"    Making sure batches are contiguous. Ordering cells based on {args.batch_key}...")
    idx = adata.obs.sort_values(args.batch_key).index
    adata = adata[idx,]

    print(f"Running scanorama on {args.embed_key}...")
    sce.pp.scanorama_integrate(
        adata,
        key=args.batch_key,
        basis=args.embed_key,
        adjusted_basis="X_scanorama",
        knn=args.n_neighbors,
        batch_size=5000 #set as large value that runs into memory
    )

    # Write output
    print(f"Saving results to {args.output}")
    adata.write_zarr(f"{args.output}")
