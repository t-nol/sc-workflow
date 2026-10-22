#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy as sc
import bbknn

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
    p.add_argument("-n", "--n_comps", required=True, type=int)
    p.add_argument("-b", "--batch_key", required=True, type=str)
    p.add_argument("-e", "--embed_key", required=False, type=str, default="X_pca")
    p.add_argument("-l", "--label_key", required=False, type=str, default=None)
    return p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    # BBKNN
    if args.label_key:
        print(f"    Running bbknn with ridge regression of {args.label_key} on {args.embed_key}...")
        bbknn.ridge_regression(adata, batch_key=args.batch_key, confounder_key=args.label_key)
        sc.tl.pca(adata, key_added=f"X_bbknn")
        bbknn.bbknn(adata, batch_key=args.batch_key, use_rep=f"X_bbknn", n_pcs=args.n_comps)
    else :
        print(f"    Running bbknn on {args.embed_key}...")
        bbknn.bbknn(adata, batch_key=args.batch_key, use_rep=args.embed_key, n_pcs=args.n_comps)
 
    # Write output
    print(f"    Saving results to {args.output}")
    adata.write_zarr(f"{args.output}")
