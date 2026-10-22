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
    return p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    # Harmony
    print(f"Running harmony on {args.embed_key}...")
    sce.pp.harmony_integrate(
        adata,
        key=args.batch_key,
        basis=args.embed_key,
        adjusted_basis="X_harmony"
    )

    # Write output
    print(f"Saving results to {args.output}")
    adata.write_zarr(f"{args.output}")
