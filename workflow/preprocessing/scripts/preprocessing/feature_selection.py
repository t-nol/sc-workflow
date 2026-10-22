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
p.add_argument("-n", "--n_top_genes", required=True, type=int)
p.add_argument("-b", "--batch_key", required=False, default=None, type=str)
args = p.parse_args()

# ----------------------------
# Functions
# ----------------------------
def feature_selection(adata, selection_method, n_genes, batch_key):
    
    if selection_method == "highly_variable_genes":
        selected_genes = sc.pp.highly_variable_genes(
            adata, 
            layer='norm_counts',
            flavor='seurat',
            n_top_genes=n_genes, 
            batch_key=batch_key,
            check_values=True,
            inplace=False,
            )
    
    if selection_method == "highly_deviant_genes":
        selected_genes = highly_deviant_genes(
            adata=adata,
            layer="norm_counts", 
            top_pct=None, 
            top_n=None, 
            family='binomial', 
            gene_names="gene_name",
            )
        
    return selected_genes[['means','dispersions','dispersions_norm','highly_variable']] #need the first 2 for plottng in report


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     

    # Read inputs
    adata = build_adata(args.input)

    # Feature selection
    genes_df = feature_selection(adata, args.method, args.n_top_genes, args.batch_key)
    adata.var = adata.var.join(genes_df)
    
    # Write output
    adata.write_zarr(args.output)
