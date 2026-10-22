#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy as sc
import torch
import concord as ccd

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
    return p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    if args.feature_space == "hvg":
        feature_list = adata.var.query("highly_variable").index.to_list()
        preload_dense = True
    else:
        feature_list = adata.var_names.to_list()
        preload_dense = False
    
    # CONCORD
    ccd_obj = ccd.Concord(
        adata=adata, 
        input_feature=feature_list, 
        domain_key=args.batch_key, 
        latent_dim=50,
        seed=0, 
        preload_dense=preload_dense,
        verbose=True, 
        device='cpu'
    ) 

    # Encode data, saving the latent embedding in adata.obsm
    ccd_obj.fit_transform(output_key="X_concord")
    
    # Write output
    print(f"\nSaving results to {args.output}")
    ccd_obj.adata.write_zarr(f"{args.output}")
