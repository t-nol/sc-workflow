import argparse
from pathlib import Path
import scanpy as sc
import celltypist
import sys

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=Path)
parser.add_argument("-o", "--output", required=True, type=Path)
parser.add_argument("--gene-id-column", required=True)
parser.add_argument("--cell-type-column", required=True)
parser.add_argument("--feature-selection", action="store_true")
parser.add_argument("--top-genes", type=int, required=True)
parser.add_argument("--check-expression", action="store_true")
parser.add_argument("-n", "--n-jobs", type=int)
args = parser.parse_args()


def train(adata):
    adata.var_names = adata.var[args.gene_id_column].astype(str)

    if args.cell_type_column not in adata.obs:
      raise ValueError("Cell type column not found")

    if not args.check_expression:
      sc.pp.normalize_total(adata, target_sum=1e4)
      sc.pp.log1p(adata)
      
    model = celltypist.train(
      X=adata,
      labels=args.cell_type_column,
      feature_selection=args.feature_selection,
      top_genes=args.top_genes,
      check_expression=args.check_expression,
      n_jobs=args.n_jobs,
    )

    return model

if __name__ == "__main__":

    adata = build_adata(raw_path=args.input)
    model = train(adata)
    model.write(str(args.output))
