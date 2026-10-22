#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import anndata as ad
import numpy as np

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-i", "--input", required=True, type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
args = p.parse_args()

# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":     
    
    # Read inputs
    adata = ad.read_zarr(args.input)

    # Write outputs
    adata.write_h5ad(args.output)
