import argparse
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
parser.add_argument("-s", required=False, type=bool, default=False)
args = parser.parse_args()

input_path = args.i
output_path = args.o
markers_path = args.m
dendro_path = args.d
cluster_key = args.k
clusters = args.c
n_markers = args.n
subclustering = args.s

# ----------------------------
# Load data
# ----------------------------
adata = ad.read_zarr(input_path)
clusters = pd.read_csv(clusters, index_col=0)
adata.obs[cluster_key] = clusters[cluster_key].astype("category")

# Load markers & pick top N per group
markers = pd.read_csv(markers_path)
markers_top = (
    markers.groupby('group')
    .apply(lambda x: x.nlargest(n_markers, 'scores'))
    .reset_index(drop=True)
)

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

    # --------------------------------------------------------------------
    # either main clustering : single figure
    # or subclustering : make one figure per parent cluster
    # --------------------------------------------------------------------

    if not subclustering:
        # MAIN CLUSTERING
        to_process = {"ALL": markers_top}

    else:
        # SUB-CLUSTERING
        parent_clusters = sorted(
            {lab.split(".")[0] for lab in adata.obs[cluster_key].astype(str)}
        )
        
        to_process = {}
        for parent in parent_clusters:
            # Extract subcluster labels belonging to this parent
            sub_groups = [
                g for g in adata.obs[cluster_key].cat.categories.tolist()
                if str(g).startswith(f"{parent}.")
            ]

            mk = markers_top[markers_top["group"].isin(sub_groups)]
            if len(mk) == 0:
                continue

            to_process[parent] = mk


    # --------------------------------------------------------------------
    # Plotting
    # --------------------------------------------------------------------
    for parent_label, mk in to_process.items():

        title = (
            f"{cluster_key} (all clusters)"
            if parent_label == "ALL"
            else f"{cluster_key} – parent {parent_label}"
        )

        genes_by_cluster = {
            str(k): v
            for k, v in mk.groupby("group")["gene_name"].apply(list).to_dict().items()
        }

        # Heatmap -------------------
        # only for main markers, not very visible on children cluster specific markers
        if parent_label == "ALL":
            sc.pl.heatmap(
                adata,
                var_names=genes_by_cluster,
                groupby=cluster_key,
                gene_symbols="gene_name",
                layer="norm_counts",
                standard_scale="var",
                cmap="RdBu_r",
                swap_axes=True,
                show=False,
            )
            pdf.savefig(bbox_inches="tight")
            plt.close()

        # Dotplot -------------------
        dp = sc.pl.dotplot(
            adata,
            var_names=genes_by_cluster,
            groupby=cluster_key,
            gene_symbols="gene_name",
            layer="norm_counts",
            standard_scale="var",
            cmap="Reds",
            #dendrogram=cluster_key if parent_label == "ALL" else None,
            title=f"{title}",
            return_fig=True,
            show=False,
        )
        dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5)
        fig = dp.make_figure()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
        del dp
