#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import sys
import numpy as np
import plotly.express as px

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-i", "--inputs", required=True, nargs="+", type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
args = p.parse_args()

# ----------------------------
# Metric groups (ORDER DEFINES COLUMN ORDER)
# ----------------------------
BIO_METRICS = [
    "hvg_overlap",
    "cell_cycle",
    "silhouette",
    "isolated_labels_asw",
    "clisi_graph",
    "isolated_labels_f1",
    "ari",
    "nmi",
    "trajectory_conservation",
]

BATCH_METRICS = [
    "pcr_comparison",
    "silhouette_batch",
    "ilisi_graph",
    "kBET",
    "graph_connectivity",
]

# ----------------------------
# Helpers
# ----------------------------
def load_metrics(csv_path: Path) -> pd.DataFrame | None:
    if not csv_path.exists():
        print(f"[WARN] Missing file: {csv_path}", file=sys.stderr)
        return None
    try:
        df = pd.read_csv(csv_path, index_col=0)
        if df.empty:
            print(f"[WARN] Empty CSV skipped: {csv_path}", file=sys.stderr)
            return None
        return df
    except Exception as e:
        print(f"[ERROR] Failed to read {csv_path}: {e}", file=sys.stderr)
        return None


def minmax_scale(df: pd.DataFrame) -> pd.DataFrame:
    """Min–max scale each metric independently, ignoring NaNs."""
    scaled = df.copy()
    for col in scaled.columns:
        values = scaled[col].dropna()
        if values.empty:
            continue
        vmin, vmax = values.min(), values.max()
        if vmin == vmax:
            scaled[col] = np.nan
        else:
            scaled[col] = (scaled[col] - vmin) / (vmax - vmin)
    return scaled


# ----------------------------
# Main
# ----------------------------
if __name__ == "__main__":
    dfs = []

    for csv_path in args.inputs:
        df = load_metrics(csv_path)
        if df is not None:
            dfs.append(df)

    if not dfs:
        raise RuntimeError("No valid scIB metric files found.")

    # ----------------------------
    # Combine metrics
    # ----------------------------
    combined = pd.concat(dfs, axis=0, sort=True)
    combined = combined.sort_index()

    for col in combined.columns:
        try:
            combined[col] = pd.to_numeric(combined[col])
        except (ValueError, TypeError):
            pass

    combined.index.name = "integration"

    # ----------------------------
    # Reorder columns by category
    # ----------------------------
    ordered_columns = (
        [c for c in BIO_METRICS if c in combined.columns]
        + [c for c in BATCH_METRICS if c in combined.columns]
    )

    combined = combined[ordered_columns]

    # ----------------------------
    # Write CSV
    # ----------------------------
    combined.to_csv(args.output)
    print(f"\n[INFO] Combined scIB metrics CSV written to {args.output}")
