
#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy.external as sce
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
    return p.parse_args()

# ----------------------------
# Function
# ----------------------------
def run_integration(adata, batch_key, feature_space):
    if feature_space == "hvg":
        vargenes = adata.var["highly_variable"].values
    elif feature_space == "full":
        vargenes = None

    batches = adata.obs[batch_key].cat.categories.tolist()
    
    corrected_data = sce.pp.mnn_correct(
        adata, 
        svd_dim = 50, 
        batch_key = batch_key, 
        batch_categories = batches,
        var_subset = vargenes,
        do_concatenate = True,
        )
    
    corrected_adata = corrected_data[0][0]
    corrected_counts = csr_matrix(corrected_adata.X)
    
    return corrected_counts

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    # MNN
    corrected_counts = run_integration(adata, args.batch_key, args.feature_space)
    adata.layers["mnn_counts"] = corrected_counts

   # Write output
    print(f"    Saving results to {args.output}")
    adata.write_zarr(f"{args.output}")
    