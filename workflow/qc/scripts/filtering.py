#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata


# ----------------------------
# Helper Functions
# ----------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Filter AnnData based on QC metrics")
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("-f", "--db_filter", required=True, type=str)
    return parser.parse_args()


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()

    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    # Filter out cells / features
    kept_cells = (adata.obs['qc_pass_global'] == True).sum()
    removed_cells = (adata.obs['qc_pass_global'] == False).sum()
    adata = adata[adata.obs['qc_pass_global'] == True].copy()
    print(f"Filtering out {removed_cells} cells. {kept_cells} cells passed QC.")

    kept_features = (adata.var['qc_pass_global'] == True).sum()
    removed_features = (adata.var['qc_pass_global'] == False).sum()
    adata = adata[:, adata.var['qc_pass_global'] == True].copy()
    print(f"Filtering out {removed_features} features. {kept_features} features passed QC.")

    # Filter out doublets
    if args.db_filter == "all":
        adata = adata[~adata.obs["is_doublet_all"]].copy()
        print(f"Filtered doublets identified with ALL tools: {adata.n_obs} cells remain.")
    elif args.db_filter == "any":
        adata = adata[~adata.obs["is_doublet_any"]].copy()
        print(f"Filtered doublets identified with ANY (at least one) tool: {adata.n_obs} cells remain.")
    else:
        print("No doublet filtering applied. Keeping all cells.")

    # Save output
    print(f"Saving prepared AnnData to {args.output}...")
    adata.write_zarr(args.output)
    print("Done.")