import argparse
import yaml
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, type=Path)
parser.add_argument("-c", required=True, type=Path)
parser.add_argument("-m", required=True, type=Path)
parser.add_argument("-d", required=True, type=Path)
parser.add_argument("-k", required=True, type=str)
parser.add_argument("-o", required=True, type=Path)
parser.add_argument("-n", required=True, type=int)
args = parser.parse_args()

input_path = args.i
output_path = args.o
markers_path = args.m
dendro_path = args.d
cluster_key = args.k
clusters = args.c
n_markers = args.n

# ----------------------------
# Load data
# ----------------------------
adata = ad.read_zarr(input_path)
clusters = pd.read_csv(clusters, index_col=0)
adata.obs[cluster_key] = clusters[cluster_key].astype("category")

# Load markers
with open(markers_path) as f:
    markers = yaml.safe_load(f)

valid_genes = set(adata.var["gene_name"]) # keep only markers defined in the adata
markers = {
    ct: [g for g in genes if g in valid_genes]
    for ct, genes in markers.items()
}


# Load dendrogram
with open(dendro_path, "rb") as f:
    dendro = pickle.load(f)

adata.uns[f'dendrogram_{cluster_key}'] = {
    'linkage': dendro.get('linkage'),
    'reordered_obs_names': dendro.get(
        'reordered_obs_names',
        list(adata.obs[cluster_key].cat.categories)
    )
}

# ----------------------------
# Create output file
# ----------------------------
output_path.parent.mkdir(parents=True, exist_ok=True)

with PdfPages(output_path) as pdf:

    # Dotplot -------------------
    dp = sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby=cluster_key,
        gene_symbols="gene_name",
        layer="norm_counts",
        standard_scale="var",
        cmap="Reds",
        return_fig=True,
        show=False,
    )
    dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5)
    fig = dp.make_figure()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

