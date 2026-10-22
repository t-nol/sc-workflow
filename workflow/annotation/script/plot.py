import argparse
from pathlib import Path
import sys
import scanpy as sc
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import warnings
warnings.filterwarnings('ignore')

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

parser = argparse.ArgumentParser(description="Plot CellTypist annotation using existing UMAP or create fallback")
parser.add_argument("-i", "--input", required=True, type=Path, help="Input h5ad file")
parser.add_argument("-p", "--plot", required=True, type=Path, help="Output plot file")
parser.add_argument("--use-existing-umap", action="store_true", help="Use existing UMAP if present")
parser.add_argument("--umap-key", type=str, required=True, help="Key of UMAP in adata.obsm")

args = parser.parse_args()

# ------------------------------------
# Main
# ------------------------------------

if __name__ == "__main__":
    adata = build_adata(raw_path=args.input)
    
    # Get umap key
    if args.use_existing_umap and args.umap_key in adata.obsm:
        umap_key = args.umap_key
    else:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
        sc.tl.umap(adata)
        umap_key = "X_umap"

    # Find all annotation columns
    annotation_cols = [col for col in adata.obs.columns if col.startswith("celltypist_labels_")]
    if not annotation_cols:
        raise ValueError("No columns starting with 'celltypist_labels_' found in adata.obs")
    
    # One plot for each annotation levels
    with PdfPages(args.plot) as pdf:
        for col in annotation_cols:
            level = col.replace("celltypist_labels_", "")  # extract annotation level
            sc.pl.embedding(
                adata,
                basis=umap_key.replace("X_", ""),
                color=col,
                legend_loc="right margin",
                show=False,
                title=f"CellTypist Annotation ({level})"
            )
            plt.tight_layout()
            pdf.savefig(dpi=300, bbox_inches="tight")
            plt.close()
