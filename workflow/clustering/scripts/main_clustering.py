import argparse
from pathlib import Path
import anndata as ad
import scanpy as sc

# ----------------------------
# Arguments
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, type=Path, help="Input .zarr or .h5ad")
parser.add_argument("-o", required=True, type=Path, help="Output CSV")
parser.add_argument("-a", required=True, type=str, help="Clustering algorithm")
parser.add_argument("-r", required=True, type=float, help="Resolution")
args = parser.parse_args()

input_path = args.i
output_csv = args.o
alg = args.a.lower()
res = args.r

# ----------------------------
# Load input AnnData
# ----------------------------
if input_path.suffix == ".zarr":
    adata = ad.read_zarr(input_path)
elif input_path.suffix == ".h5ad":
    adata = ad.read_h5ad(input_path)
else:
    raise ValueError(f"Unsupported format: {input_path}")

# ----------------------------
# Run clustering
# ----------------------------
main_key = f"{alg}_{res}"

if alg == "leiden":
    sc.tl.leiden(adata, resolution=res, key_added=main_key, flavor="igraph")
elif alg == "louvain":
    sc.tl.louvain(adata, resolution=res, key_added=main_key)
else:
    raise ValueError(f"Unsupported clustering algorithm: {alg}")

# ----------------------------
# Save output CSV
# ----------------------------
output_csv.parent.mkdir(parents=True, exist_ok=True)
adata.obs[main_key].to_csv(output_csv)
