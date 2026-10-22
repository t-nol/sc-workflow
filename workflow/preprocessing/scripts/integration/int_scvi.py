#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import scvi

project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True, type=Path)
    p.add_argument("-o", "--output", required=True, type=Path)
    p.add_argument("-b", "--batch_key", required=True, type=str)
    p.add_argument("-l", "--label_key", required=False, type=str, default=None)
    p.add_argument("-f", "--feature_space", required=True, type=str)
    return p.parse_args()

# ----------------------------
# Function
# ----------------------------
def run_integration(adata, batch_key, label_key, latent_key):  
    print(f"    Running scVI...")    
    scvi.model.SCVI.setup_anndata(adata, layer="raw_counts", batch_key=batch_key)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    model.train(
        max_epochs=400,
        early_stopping=True,
        check_val_every_n_epoch=10,
    )
    
    if latent_key=="X_scvi":
        adata.obsm[latent_key] = model.get_latent_representation(adata)
    
    elif latent_key=="X_scanvi":
        print(f"    Running scANVI...")    
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            model,
            adata=adata,
            labels_key=label_key,
            unlabeled_category="Unknown",
        )
        scanvi_model.train(
            max_epochs=200,
            early_stopping=True,
            check_val_every_n_epoch=10,
            n_samples_per_label=100
            )
        adata.obsm[latent_key] = scanvi_model.get_latent_representation(adata)

    return(adata)

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # Load input
    print(f"Loading AnnData from {args.input}...")
    adata = build_adata(raw_path=args.input)

    if args.feature_space == "hvg":
        adata_int = adata[:, adata.var["highly_variable"]].copy()
    else:
        adata_int = adata

    # scVI/scANVI
    if args.label_key:
        latent_key = "X_scanvi"
        adata_int.obs[args.label_key] = adata_int.obs[args.label_key].astype("category")
        if adata_int.obs[args.label_key].isna().sum() > 0 :
            adata_int.obs[args.label_key] = adata_int.obs[args.label_key].cat.add_categories("Unknown")
            adata_int.obs[args.label_key] = adata_int.obs[args.label_key].fillna("Unknown")
    else:  
        latent_key = "X_scvi"
    
    adata_int = run_integration(adata_int, args.batch_key, args.label_key, latent_key)

    adata.obsm[latent_key] = adata_int.obsm[latent_key]

    # Write output
    print(f"    Saving results to {args.output}")
    adata.write_zarr(f"{args.output}")
