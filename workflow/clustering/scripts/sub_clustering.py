import argparse
from pathlib import Path
import scanpy as sc
import anndata as ad
import pandas as pd

# ----------------------------
# Arguments
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, type=Path)
parser.add_argument("-a", required=True, type=str)
parser.add_argument("-m", required=True, type=Path)
parser.add_argument("-s", required=True, type=float)
parser.add_argument("-o", required=True, type=Path)
args = parser.parse_args()

input_path = args.i
algorithm = args.a
main_csv_path = args.m
sub_res = args.s
output_csv = args.o

# ----------------------------
# Load input
# ----------------------------
if input_path.suffix == ".zarr":
    adata = ad.read_zarr(input_path)
elif input_path.suffix == ".h5ad":
    adata = ad.read_h5ad(input_path)
else:
    raise ValueError(f"Unsupported format: {input_path}")

main_clusters = pd.read_csv(main_csv_path, index_col=0, dtype=str)
main_key = main_clusters.columns[0]
adata.obs[main_key] = main_clusters[main_key].astype("category")

# ----------------------------
# Prepare output column
# ----------------------------
sub_key = f"{main_key}_sub_{sub_res}"
adata.obs[sub_key] = adata.obs[main_key].astype(str)

# Temporary column used during restrict_to clustering
tmp_key = f"{sub_key}_tmp"
adata.obs[tmp_key] = pd.Categorical([""] * adata.n_obs)

# ----------------------------
# Run subclustering per parent cluster
# ----------------------------
for cluster in adata.obs[main_key].cat.categories:
    mask = adata.obs[main_key] == cluster

    # Subclustering only in the subset
    if algorithm == "leiden":
        sc.tl.leiden(
            adata,
            restrict_to=[main_key, [cluster]],
            resolution=float(sub_res),
            key_added=tmp_key,
            flavor="igraph",
        )
    elif algorithm == "louvain":
        sc.tl.louvain(
            adata,
            restrict_to=[main_key, [cluster]],
            resolution=float(sub_res),
            key_added=tmp_key,
        )
    else:
        raise ValueError(f"Unsupported algorithm: {algorithm}")

    # convert "1,0" → "1.0" and save to sub_key
    adata.obs.loc[mask, sub_key] = adata.obs.loc[mask, tmp_key].astype(str).str.replace(",", ".")

del adata.obs[tmp_key]

# ----------------------------
# Save CSV
# ----------------------------
output_csv.parent.mkdir(parents=True, exist_ok=True)
adata.obs[sub_key].to_csv(output_csv, header=True)
