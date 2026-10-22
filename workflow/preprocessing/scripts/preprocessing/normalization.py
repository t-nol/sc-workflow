#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
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
p.add_argument("-m", "--method", required=True, type=str)
args = p.parse_args()

# ----------------------------
# Functions
# ----------------------------
def shifted_log(adata):
    sc.pp.normalize_total(adata, target_sum=10000, inplace=True)
    sc.pp.log1p(adata)
    return adata

def pearson_residuals(adata):
    sc.experimental.pp.normalize_pearson_residuals(adata, inplace=True)
    return adata

def normalize(adata, norm_method):
    adata.uns['normalization_method'] = norm_method
    if norm_method == "shifted_logarithm":
        return shifted_log(adata)
    if norm_method == "pearson_residuals":
        return pearson_residuals(adata)
    raise ValueError(f"Unknown method: {norm_method}")

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Read inputs
    adata = build_adata(raw_path=args.input)

    # Normalize
    adata.X = adata.layers["raw_counts"].copy()
    adata = normalize(adata, norm_method=args.method)
    adata.layers["norm_counts"] = adata.X

    # Write output
    adata.write_zarr(args.output)
