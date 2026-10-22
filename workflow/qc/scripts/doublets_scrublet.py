#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata


# ----------------------------
# Helper Functions
# ----------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("-o", "--output", required=True, type=Path)
    parser.add_argument("-b", "--batch_key", default=None, type=str)
    return parser.parse_args()


def run_scrublet(adata):
    max_comps = min(adata.n_obs, adata.n_vars) - 2

    sc.pp.scrublet(
        adata,
        expected_doublet_rate=0.05,
        stdev_doublet_rate=0.02,
        n_prin_comps=min(30, max_comps)
    )

    return (
        adata.obs[["doublet_score", "predicted_doublet"]]
        .rename(
            columns={
                "doublet_score": "scrublet_doublet_score",
                "predicted_doublet": "scrublet_doublet_prediction",
            }
        )
    )


def qc_filtering(adata, batch):
    col = f"qc_pass_{batch}"
    # force proper boolean dtype
    adata.obs[col] = adata.obs[col].to_numpy(dtype=bool)
    adata.var[col] = adata.var[col].to_numpy(dtype=bool)
    # filter
    adata = adata[adata.obs[col], adata.var[col]].copy()
    print("Post filtering :", adata)
    return adata


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()

    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    # ----------------------------
    # No batch splitting
    # ----------------------------
    if args.batch_key == "None":
        print("\n\nRunning Scrublet without batch splitting.")
        adata = qc_filtering(adata, batch="global")
        
        # Skip empty batches
        if adata.n_obs == 0 or adata.n_vars == 0:
            print(f"    -> Skipping: no cells or genes left after QC.")
            doublet_df = pd.DataFrame(columns=["scrublet_doublet_score","scrublet_doublet_prediction"])
        
        # Skip very small batches
        elif adata.n_obs < 50:
            print(f"    -> Skipping scrublet: too small, less than 50 cells left after QC.")
            doublet_df = pd.DataFrame(columns=["scrublet_doublet_score","scrublet_doublet_prediction"])
        
        else :
            doublet_df = run_scrublet(adata)
            print(f"    -> Found {(doublet_df.scrublet_doublet_prediction).sum()} doublets.")

    # ----------------------------
    # Batch-aware Scrublet
    # ----------------------------
    else:
        print(f"\n\nRunning Scrublet per batch using key: {args.batch_key}")
        batches = adata.obs[args.batch_key].astype("category")

        dfs = []

        for batch in batches.cat.categories:
            print(f"\n\n  Processing batch: {batch}")
            adata_b = adata[batches == batch].copy()
            adata_b = qc_filtering(adata_b, batch)
            
            # Skip empty batches
            if adata_b.n_obs == 0 or adata_b.n_vars == 0:
                print(f"    -> Skipping batch {batch}: no cells or genes left after QC.")
                continue
            
            # Skip very small batches
            if adata_b.n_obs < 50:
                print(f"    -> Skipping batch {batch}: too small, less than 50 cells left after QC.")
                continue

            df_b = run_scrublet(adata_b)
            dfs.append(df_b)
            print(f"\n    -> Found {(df_b.scrublet_doublet_prediction).sum()} doublets.")

        # row-wise merge (cells are unique)
        doublet_df = pd.concat(dfs, axis=0)

    # ----------------------------
    # Save output
    # ----------------------------
    print(f"Saving results to {args.output}...")
    doublet_df.to_csv(args.output, header=True, index=True)
    print("Done.")
