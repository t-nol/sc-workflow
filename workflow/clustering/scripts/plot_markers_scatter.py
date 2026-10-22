#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
from matplotlib.lines import Line2D

from adjustText import adjust_text
adjust_text.debug = False

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--markers_csv", required=True, type=Path)
parser.add_argument("-k", "--clust_key", required=True, type=str)
parser.add_argument("-o", "--output", required=True, type=Path)

# thresholds
parser.add_argument("-p", "--adj_pval", type=float, default=0.05)
parser.add_argument("-l", "--logfc", type=float, default=1.5)

args = parser.parse_args()

markers = pd.read_csv(args.markers_csv)
args.output.parent.mkdir(parents=True, exist_ok=True)


def neglog10(x):
    return -np.log10(np.clip(x, 1e-300, 1))


# Fill missing parent_cluster with "rest_of_cells"
if "parent_cluster" not in markers.columns:
    markers["parent_cluster"] = "rest_of_cells"
else:
    markers["parent_cluster"] = markers["parent_cluster"].fillna("rest_of_cells")


with PdfPages(args.output) as pdf:
    for parent, parent_df in markers.groupby("parent_cluster"):

        child_clusters = parent_df["group"].unique()
        n_children = len(child_clusters)

        # Determine page title: only for subclusters
        page_title = f"Parent cluster: {parent}" if n_children > 1 else ""

        # Layout: max 3 columns
        n_cols = min(3, n_children)
        n_rows = math.ceil(n_children / n_cols)
        fig, axes = plt.subplots(
            nrows=n_rows, ncols=n_cols, figsize=(6*n_cols, 6*n_rows),
            squeeze=False
        )

        # Create legend elements once per page
        comparison_text = "siblings" if n_children > 1 else "rest of cells"
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=f'Upregulated vs {comparison_text}',
                   markerfacecolor='green', markersize=8),
            Line2D([0], [0], marker='o', color='w', label=f'Downregulated vs {comparison_text}',
                   markerfacecolor='red', markersize=8),
            Line2D([0], [0], marker='o', color='w', label=f'Significant but low FC',
                   markerfacecolor='lightgrey', markersize=8)
        ]

        for idx, cluster in enumerate(child_clusters):
            row_idx = idx // n_cols
            col_idx = idx % n_cols
            ax = axes[row_idx, col_idx]

            df = parent_df[parent_df["group"] == cluster].copy()
            df["neglog10p"] = neglog10(df["pvals_adj"])

            # Classify genes
            df["category"] = "grey"
            df.loc[df["logfoldchanges"] >= args.logfc, "category"] = "up"
            df.loc[df["logfoldchanges"] <= -args.logfc, "category"] = "down"
            color_map = {"up": "green", "down": "red", "grey": "lightgrey"}

            # Filter only significant genes
            sig_df = df[df["pvals_adj"] <= args.adj_pval]
            if sig_df.empty:
                ax.set_visible(False)
                continue

            # Top label candidates
            label_candidates = []
            for cat in ["up", "down"]:
                sub = sig_df[sig_df["category"] == cat]
                if len(sub) > 0:
                    label_candidates.append(sub.nsmallest(5, "pvals_adj"))
            label_df = pd.concat(label_candidates) if label_candidates else pd.DataFrame()

            # Scatter: horizontal volcano
            for cat, col in color_map.items():
                sub = sig_df[sig_df["category"] == cat]
                if len(sub) > 0:
                    ax.scatter(
                        sub["neglog10p"],     # X axis = -log10(p)
                        sub["logfoldchanges"], # Y axis = log2FC
                        s=35,
                        c=col,
                        alpha=0.8,
                        linewidth=0
                    )

            # Threshold lines (horizontal volcano)
            ax.axvline(-np.log10(args.adj_pval), color="black", ls="--", lw=1)  # vertical p-value threshold
            ax.axhline(args.logfc, color="black", ls="--", lw=1)               # horizontal FC
            ax.axhline(-args.logfc, color="black", ls="--", lw=1)

            # Symmetrically center 0 on y-axis
            max_fc = max(abs(sig_df["logfoldchanges"].min()), abs(sig_df["logfoldchanges"].max()))
            ax.set_ylim(-max_fc * 1.1, max_fc * 1.1)

            # Labels with offset
            texts = []
            for _, row in label_df.iterrows():
                x = row["neglog10p"]
                y = row["logfoldchanges"]
                label = row["gene_name"]

                dy = 0.3 if y >= 0 else -0.3
                dx = 0.3

                ax.plot([x, x + dx], [y, y + dy], color="grey", lw=0.6)
                texts.append(ax.text(
                    x + dx, y + dy, label,
                    fontsize=8,
                    ha="left",
                    va="bottom" if y >= 0 else "top"
                ))

            # Adjust labels
            if texts:
                adjust_text(
                    texts, 
                    ax=ax, 
                    arrowprops=dict(
                        arrowstyle="-", 
                        color='grey', 
                        lw=0.5,
                        shrinkA=5,
                        shrinkB=5,
                    ), 
                    debug=False)

            ax.set_xlabel("-log10(adj p-value)")
            ax.set_ylabel("log2 Fold Change")
            ax.set_title(f"{cluster} vs {comparison_text}")


        # Hide unused axes
        total_axes = n_rows * n_cols
        for idx in range(len(child_clusters), total_axes):
            row_idx = idx // n_cols
            col_idx = idx % n_cols
            axes[row_idx, col_idx].set_visible(False)

        # Page title only for subclusters
        if page_title:
            plt.suptitle(page_title, fontsize=12)

        # Add single legend above subplots
        fig.legend(
            handles=legend_elements,
            loc='upper right',
            bbox_to_anchor=(0.97, 0.97),
            fontsize=8
        )

        plt.tight_layout(rect=[0, 0, 0.95, 0.95])
        pdf.savefig(fig)
        plt.close()
