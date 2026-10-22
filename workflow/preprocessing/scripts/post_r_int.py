#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import numpy as np

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-i", "--input", required=True, type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
p.add_argument("-u", "--uncorrected", required=True, type=Path)
p.add_argument("-s", "--space", required=True, type=str)
p.add_argument("-m", "--int_method", required=True, type=str)
args = p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Read inputs
    adata = build_adata(raw_path=args.uncorrected)

    if args.space == "embedding_space":
        embed_data = np.load(args.input)
        adata.obsm[f"X_{args.int_method}"] = embed_data["X_corrected"]

    # Write outputs
    adata.write_zarr(args.output)
    