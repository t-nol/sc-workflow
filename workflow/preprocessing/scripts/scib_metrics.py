#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scib.metrics as metrics

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

import warnings
warnings.filterwarnings("ignore")

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-u", "--unintegrated", required=True, type=Path)
p.add_argument("-i", "--integrated", required=True, type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
p.add_argument("-f", "--feature_space", required=True, type=str)
p.add_argument("-b", "--batch_key", required=True, type=str)
p.add_argument("-l", "--label_key", required=False, type=str)
p.add_argument("-e", "--embed_key", required=False, type=str, default="X_pca")
p.add_argument("-m", "--integration_method", required=True, type=str)
p.add_argument("-p", "--integration_product", required=True, type=str)
args = p.parse_args()


# ----------------------------
# Helper functions
# ----------------------------
def safe_metric(func, *fargs, label_required=False, batch_required=False, **fkwargs):
    """Compute metric safely, returning np.nan if conditions not met or metric fails."""
    if label_required and not label_exists:
        return np.nan
    if batch_required and not batch_exists:
        return np.nan
    try:
        return func(*fargs, **fkwargs)
    except Exception as e:
        print(f"Warning: {func.__name__} failed: {e}")
        return np.nan


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Load AnnDatas
    adata_pre = build_adata(raw_path=args.unintegrated)
    adata_post = build_adata(raw_path=args.integrated)
    print(f"\nLoaded unintegrated: {adata_pre}")
    print(f"\nLoaded integrated: {adata_post}")

    # Check keys
    label_exists = args.label_key is not None and args.label_key in adata_post.obs
    batch_exists = args.batch_key in adata_post.obs

    # Cluster optimization if label exists
    if label_exists:
        print("\n\nLooking for optimal clustering using scib...")
        metrics.cluster_optimal_resolution(
            adata_post,
            label_key=args.label_key,
            cluster_key="scib_optimal_clustering",
            resolutions=[0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0],
            force=True,
            flavor="igraph",
            verbose=True
        )
        print(adata_post.obs["scib_optimal_clustering"])

    results = {}

    # ----------------------------
    # Feature-space metrics
    # ----------------------------
    print(f"\n\n\nComputing Feature Space metrics...")
    print(f"\n\n  hvg_overlap")
    results["hvg_overlap"] = safe_metric(
        metrics.hvg_overlap,
        adata_pre,
        adata_post,
        batch_key=args.batch_key,
        n_hvg=500,
        verbose=True,
        label_required=False,
        batch_required=batch_exists
    ) if args.feature_space == "full" and args.integration_product == "feature_space" else np.nan

    # ----------------------------
    # Embedding-space metrics
    # ----------------------------
    compute_embedding = args.embed_key != "X_pca" or args.integration_product == "feature_space" or args.integration_product == "None"
    cc_available = {"S_score", "G2M_score"}.issubset(adata_pre.obs.columns)

    embedding_metrics = {
        "cell_cycle": (metrics.cell_cycle, {"adata_pre": adata_pre, "adata_post": adata_post,
                                            "batch_key": args.batch_key, "embed": args.embed_key,
                                            "recompute_cc": False, "linreg_method": "numpy"},
                               cc_available),
        "pcr_comparison": (metrics.pcr_comparison, {"adata_pre": adata_pre, "adata_post": adata_post,
                                                    "covariate": args.batch_key, "embed": args.embed_key,
                                                    "n_comps":50, "linreg_method":"numpy", "scale":True, "n_threads":1},
                           batch_exists),
        "silhouette": (metrics.silhouette, {"adata": adata_post, "label_key": args.label_key, "embed": args.embed_key, "metric":"euclidean", "scale":True},
                       label_exists),
        "isolated_labels_asw": (metrics.isolated_labels_asw, {"adata": adata_post, "label_key": args.label_key,
                                                               "batch_key": args.batch_key, "embed": args.embed_key,
                                                               "scale":True, "verbose":True},
                                label_exists),
        "silhouette_batch": (metrics.silhouette_batch, {"adata": adata_post, "batch_key": args.batch_key,
                                                        "label_key": args.label_key, "embed": args.embed_key,
                                                        "metric":"euclidean", "scale":True, "return_all":False, "verbose":True},
                             label_exists)
    }

    print(f"\n\n\nComputing Embedding Space metrics...")
    for name, (func, fargs, cond) in embedding_metrics.items():
        print(f"\n\n  {name}")
        results[name] = safe_metric(func, **fargs) if compute_embedding and cond else np.nan

    # ----------------------------
    # Graph-space metrics
    # ----------------------------
    graph_metrics = {
        "clisi_graph": (metrics.clisi_graph, {"adata": adata_post, "label_key": args.label_key, "type_":"knn",
                                              "k0":90, "scale":True, "n_cores":1, "verbose":True}, label_exists),
        "ilisi_graph": (metrics.ilisi_graph, {"adata": adata_post, "batch_key": args.batch_key, "type_":"knn",
                                              "k0":90, "scale":True, "n_cores":1, "verbose":True}, batch_exists),
        "kBET": (metrics.kBET, {"adata": adata_post, "batch_key": args.batch_key, "label_key": args.label_key,
                                "type_":"knn", "scaled":True, "return_df":False, "verbose":True}, label_exists),
        "graph_connectivity": (metrics.graph_connectivity, {"adata": adata_post, "label_key": args.label_key}, label_exists),
        "isolated_labels_f1": (metrics.isolated_labels_f1, {"adata": adata_post, "label_key": args.label_key,
                                                             "batch_key": args.batch_key, "embed": None,
                                                             "cluster_key":"scib_optimal_clustering",
                                                             "verbose":True}, label_exists),
        "ari": (metrics.ari, {"adata": adata_post, "cluster_key":"scib_optimal_clustering",
                              "label_key": args.label_key}, label_exists),
        "nmi": (metrics.nmi, {"adata": adata_post, "cluster_key":"scib_optimal_clustering",
                              "label_key": args.label_key}, label_exists)
    }

    print(f"\n\n\nComputing Graph Space metrics...")
    for name, (func, fargs, cond) in graph_metrics.items():
        print(f"\n\n  {name}")
        results[name] = safe_metric(func, **fargs) if cond else np.nan

    # ----------------------------
    # Save results as CSV
    # ----------------------------
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results_df = pd.DataFrame([results], index=[f"{args.feature_space}_{args.integration_method}"])
    results_df.index.name = "integration"
    results_df.to_csv(args.output)
    print(f"scIB metrics saved to {args.output}")
