import argparse
from pathlib import Path
import anndata as ad
import pandas as pd
from gprofiler import GProfiler

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=Path)
parser.add_argument("-m", "--markers", required=True, type=Path)
parser.add_argument("-o", "--output", required=True, type=Path)
parser.add_argument("-n", "--top_n", type=int, default=300)
parser.add_argument("-p", "--adj_pval", type=float, default=0.05)

args = parser.parse_args()

organism = "hsapiens"
sources = [
    'GO:BP', 'GO:MF', 'GO:CC', 
    'KEGG', 'REAC', 'WP',
    'HPA', 'HP', 'TF',
    ]

# ----------------------------
# Load inputs
# ----------------------------
if args.input.suffix == ".zarr":
    adata = ad.read_zarr(args.input)
elif args.input.suffix == ".h5ad":
    adata = ad.read_h5ad(args.input)
else:
    raise ValueError(f"Unsupported format: {args.input}")

markers = pd.read_csv(args.markers)

# Check dataframe columns
if "pvals_adj" not in markers.columns:
    raise ValueError("Markers CSV must contain a 'pvals_adj' column to select genes.")

if "logfoldchanges" not in markers.columns:
    raise ValueError("Markers CSV must contain a 'logfoldchanges' column to split UP/DOWN genes.")

# select on adj_pvalues
markers = markers[markers["pvals_adj"] <= args.adj_pval]

# ----------------------------
# Run gProfiler on UP and DOWN separately
# ----------------------------
gp = GProfiler(return_dataframe=True)

results = []
for clust, sub in markers.groupby("group"):
    
    # UP-regulated genes: logFC > 0, ordered by adj p-value
    up_genes = (
        sub[sub["logfoldchanges"] > 0]
        .sort_values("pvals_adj")
        .gene_name.unique()[:args.top_n]
    )

    # DOWN-regulated genes: logFC < 0, ordered by adj p-value
    down_genes = (
        sub[sub["logfoldchanges"] < 0]
        .sort_values("pvals_adj")
        .gene_name.unique()[:args.top_n]
    )

    for direction, genes in [("UP", up_genes), ("DOWN", down_genes)]:
        if len(genes) == 0:
            continue

        enr = gp.profile(
            query=list(genes),
            organism=organism,
            sources=sources,
            user_threshold=0.05,
            significance_threshold_method="fdr",
            background=adata.var["gene_name"].tolist(),
            no_evidences=False,
        )

        if enr is None or enr.empty:
            continue
        
        # add cluster information and save
        enr["group"] = str(clust)
        enr["direction"] = direction
        enr = enr.sort_values(["source","direction", "p_value"])

        results.append(enr)

out = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
args.output.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(args.output, index=False)
