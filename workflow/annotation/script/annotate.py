import argparse
from pathlib import Path
import sys
import scanpy as sc
import celltypist
from celltypist.models import Model
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# -----------------------------------
# Argument Parsing
# -----------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=Path)
parser.add_argument("-m", "--models", required=True, type=Path, nargs="+")
parser.add_argument("-o", "--output", required=True, type=Path)
parser.add_argument("--norm-layer", required=True)
parser.add_argument("--gene-name-column", required=True)
parser.add_argument("--majority-voting", action="store_true")
parser.add_argument("--cluster_key", required=False, default=None)
parser.add_argument("--dotplot", action="store_true", default=False, help="Generate a dotplot")
parser.add_argument("--dotplot-output", required=False, type=Path, default=None)
parser.add_argument("--use-as-reference", required=False, default=None)
parser.add_argument("--use-as-prediction", default="majority_voting")
parser.add_argument("--filter-prediction", type=float, default=0.0)
parser.add_argument("--swap-axes", action="store_true")
parser.add_argument("--title", default="CellTypist label transfer")
args = parser.parse_args()

if args.cluster_key == "None":
    args.cluster_key = None

# -----------------------------------
# Functions
# -----------------------------------

def get_celltype_colname(model_path: Path):
    name = model_path.stem
    return name.replace("model.", "")


def annotation_celltypist(adata, model):
    if args.norm_layer not in adata.layers:
        raise ValueError(f"Layer '{args.norm_layer}' not found in adata.layers")
    adata.X = adata.layers[args.norm_layer].copy()
    
    model = Model.load(str(model))
    predicted = celltypist.annotate(
        adata,
        model=model,
        majority_voting=args.majority_voting,
        over_clustering=args.cluster_key,
    )

    return predicted

# -----------------------------------
# Main
# -----------------------------------

if __name__ == "__main__":
    adata = build_adata(raw_path=args.input)

    # Open PDF once if dotplot is requested
    pdf = None
    if args.dotplot and args.use_as_reference is not None and args.dotplot_output is not None:
        pdf = PdfPages(args.dotplot_output)

    for model in args.models:
        print(f"\n\nAnnotating dataset using Celltypist model : {model}\n", flush=True)

        predicted = annotation_celltypist(adata, model)
        
        celltype_colname = get_celltype_colname(model)
        if args.majority_voting:
            adata.obs[f"celltypist_labels_{celltype_colname}"] = predicted.predicted_labels["majority_voting"]
            adata.obs[f"celltypist_labels_no_voting_{celltype_colname}"] = predicted.predicted_labels["predicted_labels"]
        else:
            adata.obs[f"celltypist_labels_{celltype_colname}"] = predicted.predicted_labels["predicted_labels"]
        
        adata.obsm[f"celltypist_probability_{celltype_colname}"] = predicted.probability_matrix
        adata.obs[f"celltypist_conf_score_{celltype_colname}"] = predicted.probability_matrix.max(axis=1)

        # Dotplot
        if pdf is not None:
            print(f"\nPrinting dotplot...", flush=True)
            fig_result = celltypist.dotplot(
                predicted,
                use_as_reference=args.use_as_reference,
                use_as_prediction=args.use_as_prediction,
                filter_prediction=args.filter_prediction,
                swap_axes=args.swap_axes,
                title=f"{args.title} ({celltype_colname})",  # add annotation level to title
                show=False,
                return_fig=True,
            )
            fig_result.make_figure()
            pdf.savefig(fig_result.fig, dpi=300, bbox_inches="tight")
            plt.close(fig_result.fig)

    # Close PDF if it was opened
    if pdf is not None:
        pdf.close()

    # Save annotated adata
    print(f"\n\nSaving annotated anndata to {args.output}.", flush=True)
    adata.write_h5ad(args.output)