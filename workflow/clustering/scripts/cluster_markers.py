import argparse
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd
import pickle
from gprofiler import GProfiler


# ----------------------------
# Argument parsing
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=Path, help="Raw zarr input")
parser.add_argument("-c", "--clusters", required=True, type=Path, help="Clustering CSV")
parser.add_argument("-k", "--key", required=True, type=str, help="Clustering key name")
parser.add_argument("-o", "--output", required=True, type=Path, nargs=2, help="Output CSV and pkl paths")
parser.add_argument("-s", "--subclustering", default=False, required=False, type=bool, help="Wether to compute markers between subclusters or not")
args = parser.parse_args()

# ----------------------------
# Load AnnData and clustering
# ----------------------------
adata = ad.read_zarr(args.input)

df = pd.read_csv(args.clusters, index_col=0)
df.index = df.index.astype(str)

adata.obs[args.key] = df[args.key].astype("category")

is_subcluster = "_sub_" in args.key

# ----------------------------
# Dendrogram
# ----------------------------
sc.tl.dendrogram(adata, groupby=args.key)

# ----------------------------
# Marker detection
# ----------------------------
all_markers = []

if not args.subclustering:
    # --- Main cluster: compute markers one-vs-rest ---
    # Use all clusters and compare each to all others
    sc.tl.rank_genes_groups(
        adata,
        groupby=args.key,
        method="wilcoxon",
        use_raw=False,
        layer="norm_counts",
        key_added="dea"
    )
    adata.uns['dea'] = adata.uns['dea'].copy()  # defragment

    markers = sc.get.rank_genes_groups_df(
        adata,
        group=None,  # None means all groups
        gene_symbols="gene_name",
        key="dea",
        pval_cutoff=0.05,  # TODO: pass as config parameter
    )
    markers.rename(columns={'names': 'gene_id'}, inplace=True)
    markers['parent_cluster'] = None
    markers['subcluster'] = markers['group'].astype(str)

    all_markers.append(markers)

else:
    # --- Subcluster: compute markers within each parent cluster ---
    def split_parent_sub(label):
        parts = label.rsplit(".", 1)
        return parts[0], parts[1]

    adata.obs['parent_cluster'] = adata.obs[args.key].apply(lambda x: split_parent_sub(str(x))[0])
    adata.obs['subcluster_label'] = adata.obs[args.key].apply(lambda x: split_parent_sub(str(x))[1])

    # Compute markers within each parent cluster
    for parent in adata.obs['parent_cluster'].unique():
        sub_adata = adata[adata.obs['parent_cluster'] == parent].copy()
        subclusters = sub_adata.obs[args.key].cat.categories
        if len(subclusters) < 2:
            print(f"Skipping subcluster marker detection for parent cluster '{parent}' (only {len(subclusters)} subcluster(s))")
            continue

        sc.tl.rank_genes_groups(
            sub_adata,
            groupby=args.key,
            method="wilcoxon",
            use_raw=False,
            layer="norm_counts",
            key_added="dea"
        )
        sub_adata.uns['dea'] = sub_adata.uns['dea'].copy()

        markers_sub = sc.get.rank_genes_groups_df(
            sub_adata,
            group=None,
            gene_symbols="gene_name",
            key="dea"
        )
        markers_sub.rename(columns={'names': 'gene_id'}, inplace=True)
        markers_sub['parent_cluster'] = parent
        markers_sub['subcluster'] = markers_sub['group'].astype(str)

        all_markers.append(markers_sub)

# Concatenate all markers into one DataFrame
# Or create empty dataframe and exit if no subclusters
if len(all_markers) == 0:
    print("No marker tables were produced — skipping concatenation.")
    markers = pd.DataFrame(columns=[
        "gene_id", "group", "logfoldchanges", "pvals", "pvals_adj",
        "scores", "parent_cluster", "subcluster"
    ])
    args.output[0].parent.mkdir(parents=True, exist_ok=True)
    markers.to_csv(args.output[0], index=False)

    with open(args.output[1], "wb") as f:
        pickle.dump(adata.uns[f"dendrogram_{args.key}"], f)

    exit(0)
else:
    markers = pd.concat(all_markers, ignore_index=True)


# ----------------------------
# Marker annotation using gProfiler
# ----------------------------
gp = GProfiler(return_dataframe=True)

ensembl_ids = markers["gene_id"].dropna().unique().tolist()
if ensembl_ids:
    anno = gp.convert(
        organism="hsapiens",
        query=ensembl_ids,
        target_namespace="ENTREZGENE"
    )
else:
    anno = pd.DataFrame()


# Clean up dataframe :
# gProfiler returns rows with:
# - incoming: original ID
# - name: HGNC symbol
# - description: full name
# - converted: entrez id
# - protein: uniprot if available
anno["full_name"] = (
    anno["description"]
    .str.replace(r"\s*\[.*\]$", "", regex=True)
)
anno["full_name_source"] = (
    anno["description"]
    .str.extract(r"\[(.*)\]")
)

cols_to_drop = [
    "converted",
    "n_incoming",
    "n_converted",
    "name",
    "namespaces",
    "query",
    "description",
]
anno = anno.drop(columns=cols_to_drop, errors="ignore")

anno = anno.rename(columns={
    "incoming": "gene_id",
})

markers = markers.merge(anno, on="gene_id", how="left")

# ----------------------------
# Save results
# ----------------------------
args.output[0].parent.mkdir(parents=True, exist_ok=True)
markers.to_csv(args.output[0], index=False)

# Save dendrogram for this clustering key
with open(args.output[1], "wb") as f:
    pickle.dump(adata.uns[f"dendrogram_{args.key}"], f)
