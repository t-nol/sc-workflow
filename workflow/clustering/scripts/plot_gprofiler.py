import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import textwrap
import seaborn as sb
from matplotlib import colors, rcParams
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.backends.backend_pdf import PdfPages


# ---------------------------------------
# Utility
# ---------------------------------------
def scale_data_dynamic(data, scale_min=0.05, scale_max=0.75):
    """
    Scale data to a dynamic range per dataset.
    Returns values between scale_min and scale_max.
    """
    mind = np.min(data)
    maxd = np.max(data)
    if maxd == mind:
        maxd += 1
        mind -= 1
    dr = maxd - mind
    scaled = ((data - mind) / dr * (scale_max - scale_min)) + scale_min
    return scaled * 100  # scale to marker sizes


def plot_enrich_subplot(fig, gs, df, n_terms, cmap, direction_title):
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])

    ax1.set_title(direction_title, fontsize=14, fontweight='bold')

    if df.empty:
        ax1.text(0.5, 0.5, "No terms", ha='center', va='center')
        ax1.set_axis_off()
        return

    dfp = df.sort_values("p_value").iloc[:n_terms].copy()
    dfp["scaled_overlap"] = scale_data_dynamic(dfp["intersection_size"])
    dfp["wrapped_name"] = ["\n".join(textwrap.wrap(n, 50)) for n in dfp["name"]]

    # color scaling
    min_pval = dfp['p_value'].min()
    max_pval = dfp['p_value'].max()
    norm = colors.LogNorm(min_pval, max_pval)

    sc = ax1.scatter(
        data=dfp,
        x="recall",
        y="wrapped_name",
        c="p_value",
        cmap=plt.get_cmap(cmap + "_r"),
        norm=norm,
        linewidth=1,
        edgecolor="grey",
        s=[(i + 10) ** 1.5 for i in dfp['scaled_overlap']],
    )

    ax1.invert_yaxis()
    ax1.set_xlabel('Gene ratio')
    ax1.set_ylabel('')
    ax1.grid(axis='y')

    # dynamic legend sizes based on this plot's min/max
    min_olap = dfp['intersection_size'].min()
    max_olap = dfp['intersection_size'].max()
    olap_range = max_olap - min_olap
    size_leg_vals = [np.round(i/5)*5 for i in 
                     [min_olap, min_olap + 0.45*olap_range, max_olap]]
    size_scaled = scale_data_dynamic(np.array(size_leg_vals))
    ax2.cla()  # clear axis for legend
    for s, raw in zip(size_scaled, size_leg_vals):
        ax2.scatter([], [], s=s, color='grey', alpha=0.6, label=str(int(raw)))
    ax2.legend(title='Overlap size', scatterpoints=1, frameon=False, loc='center')
    ax2.axis("off")

    # colorbar
    cbar = fig.colorbar(sc, cax=ax3)
    cbar.set_label("p-value")


def plot_enrich_pages(df, sources, params, n_terms=10):
    """
    Yield figures with max 3 sources per page.
    """
    max_sources_per_page = 3
    # Split sources into chunks of 3
    for i in range(0, len(sources), max_sources_per_page):
        page_sources = sources[i:i + max_sources_per_page]
        fig = plt.figure(figsize=(8 * len(page_sources), 10), constrained_layout=True)
        outer_gs = GridSpec(
            2, len(page_sources), figure=fig,
            height_ratios=[1, 1], width_ratios=[1] * len(page_sources),
            hspace=0.25, wspace=0.3,
        )

        for col, src in enumerate(page_sources):
            # prep data
            srcdata = df[df["source"] == src].copy()
            mask = (
                (srcdata["term_size"] >= params[src]["min_term_size"]) &
                (srcdata["term_size"] <= params[src]["max_term_size"])
            )
            srcdata = srcdata[mask]
            up_df = srcdata[srcdata["direction"].str.lower() == "up"]
            down_df = srcdata[srcdata["direction"].str.lower() == "down"]

            # plot column subplots
            inner_up_gs = outer_gs[0, col].subgridspec(2, 2, width_ratios=[1, 0.15], wspace=0)
            plot_enrich_subplot(fig, inner_up_gs, up_df, n_terms, "Greens", f"{source_params[src]['title']}\nUP")

            inner_down_gs = outer_gs[1, col].subgridspec(2, 2, width_ratios=[1, 0.15], wspace=0)
            plot_enrich_subplot(fig, inner_down_gs, down_df, n_terms, "Reds", f"{source_params[src]['title']}\nDOWN")

        yield fig


# ---------------------------------------
# Argument parsing
# ---------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=Path)
parser.add_argument("-o", "--output", required=True, type=Path)
parser.add_argument("-n", type=int, default=10)

args = parser.parse_args()

# ---------------------------------------
# Source filters
# ---------------------------------------
source_params = {
    'GO:BP': {'title': 'GO\nBiological Process', 'min_term_size': 10, 'max_term_size': 2000},
    'GO:MF': {'title': 'GO\nMolecular Function', 'min_term_size': 10, 'max_term_size': 2000},
    'GO:CC': {'title': 'GO\nCellular Component', 'min_term_size': 10, 'max_term_size': 2000},
    'KEGG': {'title': 'KEGG', 'min_term_size': 10, 'max_term_size': 500},
    'REAC': {'title': 'Reactome', 'min_term_size': 10, 'max_term_size': 500},
    'WP':{'title': 'WikiPathways', 'min_term_size': 10, 'max_term_size': 500},
    'TF':{'title': 'TRANSFAC', 'min_term_size': 1, 'max_term_size': 75},
    'HPA': {'title': 'Human Protein Atlas', 'min_term_size': 20, 'max_term_size': 1000},
    'HP': {'title': 'Human Phenotype Ontology', 'min_term_size': 10, 'max_term_size': 400},
    'CORUM': {'title': 'CORUM Complexes', 'min_term_size': 1, 'max_term_size': 75},
}

# ---------------------------------------
# Load data
# ---------------------------------------
df = pd.read_csv(args.input)
required = ["source", "name", "p_value", "group", "intersection_size",
            "direction", "recall", "term_size"]
if not all(c in df.columns for c in required):
    raise ValueError("Missing required columns.")

clusters = df["group"].unique()
args.output.parent.mkdir(parents=True, exist_ok=True)

# ---------------------------------------
# Generate PDF (with max 3 sources per page)
# ---------------------------------------
with PdfPages(args.output) as pdf:
    for clust in clusters:
        cluster_df = df[df["group"] == clust].copy()
        sources = list(cluster_df["source"].unique())

        for fig in plot_enrich_pages(cluster_df, sources, source_params, n_terms=args.n):
            fig.text(0.01, 0.99, f"Cluster {clust}", fontsize=18, fontweight='bold', va='top', ha='left')
            pdf.savefig(fig)
            plt.close(fig)
