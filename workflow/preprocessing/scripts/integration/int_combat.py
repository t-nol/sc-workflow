
#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix

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
    p.add_argument("-f", "--feature_space", required=True, type=str)
    p.add_argument("-l", "--label_key", required=False, type=str)
    p.add_argument("-c", "--covariates", nargs="+", required=False, type=str, default=None)
    return p.parse_args()

# ----------------------------
# Functions
# ----------------------------
def run_integration(adata, batch_key, covariates):
    sc.pp.combat(
        adata,
        key=batch_key,
        covariates=covariates
    )

    corrected_counts = csr_matrix(adata.X)

    return corrected_counts

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    # Combat
    if not args.covariates:
        args.covariates = []
    if args.label_key:
        args.covariates.append(args.label_key)
    
    if args.feature_space == "hvg":
        hvgs = adata.var.query("highly_variable").index
        adata_hvg = adata[:, hvgs].copy()
        corrected_hvg = run_integration(adata_hvg, args.batch_key, args.covariates)
        # When using only hvg, corrected matrix needs to be of the same shape as the full one for loading into adata.
        # After ComBat, create a matrix of same shape as full gene set:
        # 1. Take original normalized counts, they won't be used anyway for PCA for non hvg genes afterwards,
        # 2. Replace HVGs with ComBat-corrected values.
        full_matrix = adata.X.copy()
        full_matrix[:, adata.var_names.isin(hvgs)] = corrected_hvg
        corrected_counts = full_matrix
    
    elif args.feature_space == "full":
        corrected_counts = run_integration(adata, args.batch_key, args.covariates)
    
    adata.layers["combat_counts"] = corrected_counts
    
    # Write output
    print(f"    Saving results to {args.output}")
    adata.write_zarr(f"{args.output}")
    