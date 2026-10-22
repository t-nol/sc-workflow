#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.stats import median_abs_deviation as mad

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata


# ----------------------------
# Helper Functions
# ----------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True, type=Path)
    p.add_argument("-o", "--output", required=True, type=Path)
    p.add_argument("--n_cells_min", required=True, type=int)
    p.add_argument("--n_counts_max", required=True, type=int)
    p.add_argument("--n_counts_min", required=True, type=int)
    p.add_argument("--n_genes_min", required=True, type=int)
    p.add_argument("--percent_mito_max", required=True, type=float)
    p.add_argument("--use_mad", required=True, type=str)
    p.add_argument("-b", "--batch_key", default=None, type=str)
    return p.parse_args()


def calculate_qc_metrics(adata):
    adata.var["mt"] = adata.var["gene_name"].str.startswith("MT-")
    adata.var["ribo"] = adata.var["gene_name"].str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var["gene_name"].str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        layer="raw_counts",
        inplace=True,
        percent_top=[20, 50, 100],
        log1p=True,
    )

    return adata


def filter_features(adata, batch: str, n_cells_min: int):
    gene_subset, _ = sc.pp.filter_genes(
        adata, min_cells=n_cells_min, inplace=False
    )
    adata.var[f"qc_pass_{batch}"] = gene_subset
    return adata


def is_outlier(adata, metric: str, nmads: int, upper_only=False):
    M = adata.obs[metric]
    median_val = np.median(M)
    mad_val = mad(M)

    upper_cutoff = median_val + nmads * mad_val
    lower_cutoff = median_val - nmads * mad_val

    print(f"\nMAD lower threshold for {metric}: {lower_cutoff}")
    print(f"MAD upper threshold for {metric}: {upper_cutoff}")

    if upper_only:
        return M > upper_cutoff
    return (M < lower_cutoff) | (M > upper_cutoff)


def filter_cells(
    adata,
    batch: str,
    n_counts_max: int,
    n_counts_min: int,
    n_genes_min: int,
    percent_mito_max: int,
    use_mad: str,
):
    # Step 1: initial filtering on genes per cell
    cell_subset, _ = sc.pp.filter_cells(
        adata, min_genes=n_genes_min, inplace=False
    )
    adata.obs["initial_filter"] = cell_subset

    # Step 2: outlier detection
    if use_mad == "True":
        outlier_metrics = [
            "log1p_total_counts",
            "log1p_n_genes_by_counts",
            "pct_counts_in_top_20_genes",
            "pct_counts_mt",
        ]

        for metric in outlier_metrics:
            if metric == "pct_counts_mt":
                adata.obs[f"{metric}_outlier"] = (
                    is_outlier(adata, metric, nmads=3, upper_only=True)
                    | (adata.obs[metric] > percent_mito_max)
                )
            else:
                adata.obs[f"{metric}_outlier"] = is_outlier(
                    adata, metric, nmads=5
                )
    else:
        adata.obs["total_counts_outlier"] = (
            (adata.obs.total_counts < n_counts_min)
            | (adata.obs.total_counts > n_counts_max)
        )
        adata.obs["pct_counts_mt_outlier"] = (
            adata.obs.pct_counts_mt > percent_mito_max
        )

    # Step 3: final QC pass
    outlier_cols = [c for c in adata.obs.columns if c.endswith("_outlier")]
    adata.obs[f"qc_pass_{batch}"] = (
        adata.obs["initial_filter"]
        & ~(adata.obs[outlier_cols].any(axis=1))
    )

    return adata


def process_adata(adata, batch, args):
    adata = calculate_qc_metrics(adata)
    adata = filter_features(adata, batch=batch, n_cells_min=args.n_cells_min)
    adata = filter_cells(
        adata,
        batch=batch,
        n_counts_max=args.n_counts_max,
        n_counts_min=args.n_counts_min,
        n_genes_min=args.n_genes_min,
        percent_mito_max=args.percent_mito_max,
        use_mad=args.use_mad,
    )
    return adata


def process_adata_per_batch(adata, batches, args):
    if batches==["global"]:
        return process_adata(adata, "global", args)
    
    else:
        adata_list = []
        filtered_cells = set() # union of cells across filtered batches
        filtered_features = set()  # union of features across filtered batches

        for batch in batches:
            print(f"\n\nProcessing batch: {batch}")
            adata_b = adata[adata.obs[args.batch_key] == batch].copy()
            adata_b = process_adata(adata_b, batch, args)
            adata_list.append(adata_b)
            filtered_cells.update(adata_b.obs.loc[adata_b.obs[f"qc_pass_{batch}"], :].index)
            filtered_features.update(adata_b.var.loc[adata_b.var[f"qc_pass_{batch}"], :].index)
            print(f"\n{batch} anndata object :")
            print(adata_b)

        adata = ad.concat(
            adata_list,
            join="outer",
            merge="unique",
            index_unique=None,
        )
        
        adata.obs["qc_pass_global"] = adata.obs_names.isin(filtered_cells)
        adata.var["qc_pass_global"] = adata.var_names.isin(filtered_features)

        print("\n\nFull anndata object :")
        print(adata)

        return adata


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()

    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    print(args.batch_key, type(args.batch_key))

    if args.batch_key == "None":
        print("\n\nRunning QC without batch splitting.")
        batches = ["global"]
    else:
        print(f"\n\nRunning QC per batch using key: {args.batch_key}")
        batches = adata.obs[args.batch_key].astype("category").cat.categories.to_list()
        print("    -> batches :",batches)
    
    adata = process_adata_per_batch(adata, batches, args)

    print(f"Saving prepared AnnData to {args.output}...")
    adata.write_zarr(args.output)
    print("Done.")
