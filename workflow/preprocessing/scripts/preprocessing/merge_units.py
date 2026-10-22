#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad

project_root = Path(__file__).resolve().parents[3]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Generate interactive QC HTML report")
    p.add_argument("-i", "--input", nargs="+", required=True, type=Path)
    p.add_argument("-o", "--output", required=True, type=Path)
    p.add_argument("-b", "--batch_key", required=False, type=str)
    p.add_argument("-p", "--per_batch", required=False, action='store_true')
    return p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    args = parse_args()
    
    # load and save input if only one file
    if len(args.input) == 1 :
        adata = build_adata(raw_path=args.input[0])
        print("Saving : ")
        print(adata)
        adata.write_zarr(args.output)
    
    # load and merge inputs if more
    else:
        print("Merging batches :")
        adatas = []
        for f in args.input:
            adata = build_adata(raw_path=f)
            print(adata, "\n")
            adatas.append(adata)
        
        if args.per_batch: 
            #if all batches have been processed separately, keep as much information per batch
            #note that it is better to provide several separated objects in this case to avoid
            #loosing information.
            adata_merged = ad.concat(
                adatas, 
                label=args.batch_key,
                #keys=batch_names,
                join="outer",
                merge="same",
                index_unique=None
            )
        
        print("\nInto :")
        print(adata_merged)

        adata_merged.write_zarr(args.output)
