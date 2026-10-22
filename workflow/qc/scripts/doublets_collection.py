#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import pandas as pd

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata, adata_join_without_overlaps


# ----------------------------
# Helper Functions
# ----------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("-c", "--csvs", required=True, type=Path, nargs="+")
    return parser.parse_args()

def extract_tool(path: Path) -> str:
    """
    Extract tool name from a path like: <filename>.<tool>/<batch>.csv
    """
    file = path.name
    return file.split(".", 1)[0]


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()

    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    for csv in args.csvs:
        tool = extract_tool(csv)
        print(f"    Loading doublets from {tool}...")
        doublet_df = pd.read_csv(csv, header=0, index_col=0)
        adata = adata_join_without_overlaps(adata, doublet_df, join_on='obs')

        # Count how many cells each tool flagged as doublet
        col_name = f"{tool}_doublet_prediction"
        print(f"      -> found {int(adata.obs[col_name].sum())} doublets.")
        
    # Combine results
    outlier_cols = [col for col in adata.obs.columns if col.endswith("_doublet_prediction")]
    mask = adata.obs[outlier_cols].fillna(False) #in .all, NAs from cells not passing qc are treated as True (which we don't want)
    adata.obs["is_doublet_all"] = mask.all(axis=1)
    adata.obs["is_doublet_any"] = mask.any(axis=1)

    # Track how many cells would be removed by each filter
    n_cells_before = adata.n_obs
    n_removed_all = int(adata.obs["is_doublet_all"].sum())
    n_removed_any = int(adata.obs["is_doublet_any"].sum())
    
    print(f"{n_removed_all} doublets identified with ALL tools.")
    print(f"{n_removed_any} doublets identified with at least one (ANY) tool.")

    # Save output
    print(f"Saving prepared AnnData to {args.output}...")
    adata.write_zarr(args.output)
    print("Done.")
