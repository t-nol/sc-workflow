#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd

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
p.add_argument("-l", "--gene_lists", nargs="+", type=Path, help="Gene lists: S phase and G2M phase CSV files")
args = p.parse_args()

# ----------------------------
# Functions
# ----------------------------
def scanpy_scoring(adata, s_path, g2m_paths):
    # Load cell cycle genes
    df_s_genes = pd.read_csv(s_path)
    s_genes = df_s_genes['gene_id'].tolist()

    df_g2m_genes = pd.read_csv(g2m_paths)
    g2m_genes = df_g2m_genes['gene_id'].tolist()

    # Keep only genes present in the dataset
    s_genes = [gene for gene in s_genes if gene in adata.var_names]
    g2m_genes = [gene for gene in g2m_genes if gene in adata.var_names]

    # Compute scores
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    return adata

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Read inputs
    adata = build_adata(raw_path=args.input)

    # Compute scores
    if args.method == "scanpy":
        adata = scanpy_scoring(adata, args.gene_lists[0], args.gene_lists[1])

    # Write output
    adata.write_zarr(args.output)