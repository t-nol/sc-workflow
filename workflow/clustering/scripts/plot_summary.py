import sys
from pathlib import Path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root))

import argparse
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from plot_utils import create_cluster_tree


# -----------------------------------------------------------
# Parse args
# -----------------------------------------------------------
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input-zarr", required=True, type=Path)
parser.add_argument("-m", "--main-csv", required=True, nargs="+", type=Path)
parser.add_argument("-a", "--algorithm", type=str, help="Algorithm used for clustering.")
parser.add_argument("-s", "--sub-csv", nargs="*", default=[], type=Path)
parser.add_argument("--main-res", required=True, nargs="+",
                    help="Main resolutions (one for sub-summary, many for main-summary).")
parser.add_argument("--sub-res", nargs="*", default=[],
                    help="Subclustering resolutions (optional).")
parser.add_argument("-o", "--output", required=True, type=Path)

args = parser.parse_args()


# -----------------------------------------------------------
# Convenience
# -----------------------------------------------------------
zarr_path = args.input_zarr
output_path = args.output
main_csvs = args.main_csv
sub_csvs = args.sub_csv
algorithm = args.algorithm
main_resolutions = args.main_res
sub_resolutions = args.sub_res

has_sub = len(sub_csvs) > 0
batch_key = "sample"


# -----------------------------------------------------------
# Load AnnData
# -----------------------------------------------------------
adata = ad.read_zarr(zarr_path)


# -----------------------------------------------------------
# Load ALL cluster annotations first (main + sub)
# -----------------------------------------------------------
def load_cluster_column(csv_file: Path, colname: str):
    """Add a new column in adata.obs with categorical cluster assignments."""
    df = pd.read_csv(csv_file, index_col=0)
    adata.obs[colname] = df[colname].astype("category")


# MAIN RESOLUTIONS -----------------------------------------
for res, csv_file in zip(main_resolutions, main_csvs):
    cluster_key = f"{algorithm}_{res}"
    load_cluster_column(csv_file, cluster_key)

# SUB RESOLUTIONS ------------------------------------------
if has_sub:
    main_res = main_resolutions[0]
    for sub_r, csv_file in zip(sub_resolutions, sub_csvs):
        cluster_key = f"{algorithm}_{main_res}_sub_{sub_r}"
        load_cluster_column(csv_file, cluster_key)


# -----------------------------------------------------------
# Plotting function (no more CSV loading inside!)
# -----------------------------------------------------------
def plot_one_resolution(adata, cluster_key, pdf):
    """Generate UMAP + barplots for a single resolution."""

    # -----------------------
    # 1) UMAP page
    # -----------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    ax1, ax2 = axes

    sc.pl.embedding(
        adata,
        basis="X_umap_2D",
        color=[cluster_key],
        legend_loc='on data',
        show=False,
        ax=ax1,
        title=f'UMAP: {cluster_key}'
    )

    sc.pl.embedding(
        adata,
        basis="X_umap_2D",
        color=[batch_key],
        show=False,
        ax=ax2,
        title=f'UMAP: {batch_key}'
    )

    pdf.savefig()
    plt.close()

    # -----------------------
    # 2) Proportion page
    # -----------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    ax1, ax2 = axes
    plt.subplots_adjust(wspace=0.35)

    counts = (
        adata.obs
        .groupby([cluster_key, batch_key], observed=False)
        .size()
        .unstack(fill_value=0)
    )

    # Colors
    cluster_names = counts.index.tolist()
    cluster_colors = adata.uns.get(f"{cluster_key}_colors",
                                   ["gray"] * len(cluster_names))

    batch_names = counts.columns.tolist()
    batch_colors = adata.uns.get(f"{batch_key}_colors",
                                 ["lightgray"] * len(batch_names))

    # Left: normalized by batch
    props_sample = counts.divide(counts.sum(axis=0), axis=1) * 100
    props_sample.T.plot(kind="barh", stacked=True, ax=ax1, color=cluster_colors)

    ax1.set_title("Cluster Proportions per Batch", fontweight="bold")
    ax1.set_xlabel("Percent (%)")
    ax1.set_ylabel(batch_key)
    ax1.legend(title="Cluster", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax1.grid(False)

    # Right: normalized by cluster
    props_cluster = counts.divide(counts.sum(axis=1), axis=0) * 100
    props_cluster.plot(kind="bar", stacked=True, ax=ax2, color=batch_colors)

    ax2.set_title("Batch Composition per Cluster", fontweight="bold")
    ax2.set_xlabel("Cluster")
    ax2.set_ylabel("Percent (%)")
    ax2.legend(title=batch_key, bbox_to_anchor=(1.05, 1), loc="upper left")
    ax2.grid(False)

    pdf.savefig(fig)
    plt.close()


# -----------------------------------------------------------
# Main PDF generation
# -----------------------------------------------------------
output_path.parent.mkdir(parents=True, exist_ok=True)

with PdfPages(output_path) as pdf:

    # ======================================================
    # MAIN SUMMARY
    # ======================================================
    if not has_sub:
        # cluster keys for clustree
        levels = [f"{algorithm}_{r}" for r in main_resolutions]

        clustree_fig = create_cluster_tree(adata, levels, title="Main Clustering Tree")
        clustree_fig.tight_layout()
        pdf.savefig(clustree_fig)
        plt.close(clustree_fig)

        # plot each resolution
        for res in main_resolutions:
            cluster_key = f"{algorithm}_{res}"
            plot_one_resolution(adata, cluster_key, pdf)

    else:
        # ======================================================
        # MAIN + SUB SUMMARY
        # ======================================================
        main_res = main_resolutions[0]
        main_key = f"{algorithm}_{main_res}"

        # clustree ordering
        levels = [main_key] + [f"{main_key}_sub_{sub}" for sub in sub_resolutions]

        clustree_fig = create_cluster_tree(adata, levels, title="Sub Clustering Tree")
        clustree_fig.tight_layout()
        pdf.savefig(clustree_fig)
        plt.close(clustree_fig)

        # main resolution
        plot_one_resolution(adata, main_key, pdf)

        # sub resolutions
        for sub in sub_resolutions:
            sub_key = f"{main_key}_sub_{sub}"
            plot_one_resolution(adata, sub_key, pdf)
