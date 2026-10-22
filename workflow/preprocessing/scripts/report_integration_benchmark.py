#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.io as pio
import plotly.graph_objects as go
import json

# ----------------------------
# Argument parsing
# ----------------------------
p = argparse.ArgumentParser()
p.add_argument("-i", "--input", required=True, type=Path)
p.add_argument("-o", "--output", required=True, type=Path)
args = p.parse_args()

# ----------------------------
# Metric groups
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
def minmax_scale(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        vals = out[col].dropna()
        if vals.empty:
            continue
        vmin, vmax = vals.min(), vals.max()
        if vmin != vmax:
            out[col] = (out[col] - vmin) / (vmax - vmin)
        else:
            out[col] = np.nan
    return out

# ----------------------------
# Load metrics
# ----------------------------
df = pd.read_csv(args.input, index_col=0)
df.index.name = "integration"
for c in df.columns:
    df[c] = pd.to_numeric(df[c], errors="ignore")

# ----------------------------
# Compute scores
# ----------------------------
# scale per columns
scaled = minmax_scale(df)

batch_cols = [c for c in BATCH_METRICS if c in scaled.columns]
bio_cols = [c for c in BIO_METRICS if c in scaled.columns]

scaled["Batch correction (lenient)"] = scaled[batch_cols].mean(axis=1, skipna=True)
scaled["Biological conservation (lenient)"] = scaled[bio_cols].mean(axis=1, skipna=True)
scaled["Overall (lenient)"] = 0.4*scaled["Batch correction (lenient)"] + 0.6*scaled["Biological conservation (lenient)"]

batch_cols_strict = [c for c in batch_cols if df[c].notna().all()]
bio_cols_strict = [c for c in bio_cols if df[c].notna().all()]

scaled["Batch correction (strict)"] = scaled[batch_cols_strict].mean(axis=1)
scaled["Biological conservation (strict)"] = scaled[bio_cols_strict].mean(axis=1)
scaled["Overall (strict)"] = 0.4*scaled["Batch correction (strict)"] + 0.6*scaled["Biological conservation (strict)"]

# ----------------------------
# Scatter plots
# ----------------------------
scatter_df = scaled.reset_index()

def make_scatter(xcol, ycol, colorcol, title):
    fig = px.scatter(
        scatter_df,
        x=xcol,
        y=ycol,
        text="integration",
        color=colorcol,
        color_continuous_scale="viridis",
        title=title,
        hover_data=scatter_df.columns
    )
    fig.update_traces(textposition="top center", marker=dict(size=12, line=dict(width=1, color="black")))
    fig.update_layout(width=700, height=600, xaxis=dict(range=[0,1]), yaxis=dict(range=[0,1]), plot_bgcolor="white")
    return pio.to_html(fig, include_plotlyjs=True, full_html=False)

scatter_lenient_html = make_scatter("Batch correction (lenient)", "Biological conservation (lenient)", "Overall (lenient)", "Batch vs Biology (lenient)")
scatter_strict_html = make_scatter("Batch correction (strict)", "Biological conservation (strict)", "Overall (strict)", "Batch vs Biology (strict)")

# ----------------------------
# Overview dotplot
# ----------------------------
heatmap_df = scaled.copy()
heatmap_df = heatmap_df.loc[:, heatmap_df.notna().any()]

groups = [
    ("Batch correction (label-free)", [c for c in ["pcr_comparison", "ilisi_graph"] if c in heatmap_df.columns]),
    ("Batch correction (with label)", [c for c in ["silhouette_batch", "kBET", "graph_connectivity"] if c in heatmap_df.columns]),
    ("Biology (label-free)", [c for c in ["hvg_overlap", "cell_cycle", "silhouette"] if c in heatmap_df.columns]),
    ("Biology (with label)", [c for c in ["isolated_labels_asw","isolated_labels_f1","clisi_graph","ari","nmi","trajectory_conservation"] if c in heatmap_df.columns]),
    ("Lenient summary", [c for c in ["Batch correction (lenient)", "Biological conservation (lenient)", "Overall (lenient)"] if c in heatmap_df.columns]),
    ("Strict summary", [c for c in ["Batch correction (strict)", "Biological conservation (strict)", "Overall (strict)"] if c in heatmap_df.columns])
]

plotted_cols = [c for _, cols in groups for c in cols]

# colors
group_colors = {
    "Batch correction (label-free)": "#0080ff",
    "Batch correction (with label)": "#102d66",
    "Biology (label-free)": "#00b300",
    "Biology (with label)": "#008040",
    "Lenient summary": "#ff7f00",
    "Strict summary": "#ff0000"
}

col_to_group = {}
for group_name, cols in groups:
    for c in cols:
        col_to_group[c] = group_name

# Default row order (as loaded)
default_row_order = heatmap_df.index.tolist()

# Build dotplot figure
def build_dotplot(row_order):
    fig = go.Figure()

    n_rows = len(row_order)
    n_cols = len(plotted_cols)

    for j, col in enumerate(plotted_cols):
        group_name = col_to_group[col]
        dot_color = group_colors[group_name]
        for i, row in enumerate(row_order):
            val = heatmap_df.loc[row, col]
            if pd.isna(val):
                continue
            dot_size = 10 + val * 30

            fig.add_trace(go.Scatter(
                x=[j],
                y=[i],
                mode="markers+text",
                marker=dict(
                    size=dot_size,
                    color=dot_color,
                    opacity=0.3 + val * 0.7,
                    line=dict(width=1, color="black")
                ),
                text=f"{val:.2f}",
                textposition="middle center",
                textfont=dict(size=9, color="black" if val < 0.6 else "white"),
                hovertemplate=f"<b>{row}</b><br>{col}<br>Score: {val:.3f}<extra></extra>",
                showlegend=False
            ))

    x_start = 0
    for group_name, cols in groups:
        n_cols_g = len(cols)
        if n_cols_g == 0:
            continue
        fig.add_shape(
            type="rect",
            x0=x_start - 0.5,
            x1=x_start + n_cols_g - 0.5,
            y0=n_rows - 0.5,
            y1=n_rows + 3,
            line=dict(width=0),
            fillcolor=group_colors[group_name],
            layer="above"
        )
        fig.add_annotation(
            x=x_start + n_cols_g / 2 - 0.5,
            y=n_rows - 0.1,
            text=group_name.replace("Batch correction (label-free)","Batch<br>correction<br>(label-free)"),
            showarrow=False,
            xanchor="center",
            yanchor="middle",
            font=dict(color="white", size=11, family="Arial Bold")
        )
        x_start += n_cols_g

    fig.update_layout(
        width=1300,
        height=100 + n_rows * 56 + 200,
        xaxis=dict(
            tickvals=list(range(len(plotted_cols))),
            ticktext=plotted_cols,
            tickangle=-55,
            tickfont=dict(size=13),
            showgrid=False,
            zeroline=False,
            range=[-0.5, len(plotted_cols) - 0.5],
            automargin=True,
        ),
        yaxis=dict(
            tickvals=list(range(n_rows)),
            ticktext=row_order,
            tickfont=dict(size=13),
            showgrid=True,
            gridcolor="#eeeeee",
            zeroline=False,
            range=[-0.5, n_rows + 0.5]
        ),
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=60, l=220, r=20, b=250),
    )
    return fig

fig_dot = build_dotplot(default_row_order)

all_sorted_orders = {}
for col in plotted_cols:
    sorted_idx = heatmap_df.sort_values(col, ascending=True).index.tolist()
    all_sorted_orders[col] = sorted_idx

sort_figures_json = {}
for col in plotted_cols:
    row_order = all_sorted_orders[col]
    fig_s = build_dotplot(row_order)
    sort_figures_json[col] = fig_s.to_json()

default_fig_json = fig_dot.to_json()

overview_dotplot_html = f"""
<div>
  <label style="font-weight:bold; margin-right:8px;">Sort rows by:</label>
  <select id="sort-select" onchange="updateDotplot(this.value)" style="padding:4px 8px; font-size:14px; border-radius:4px; border:1px solid #ccc;">
    <option value="__default__">-- default order --</option>
    {"".join(f'<option value="{col}">{col}</option>' for col in plotted_cols)}
  </select>
</div>
<div id="dotplot-container" style="margin-top:10px;"></div>

<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<script>
var sortFigures = {json.dumps(sort_figures_json)};
var defaultFig = {default_fig_json};

function renderFig(figJson) {{
  var fig = (typeof figJson === 'string') ? JSON.parse(figJson) : figJson;
  Plotly.newPlot('dotplot-container', fig.data, fig.layout, {{responsive: true}});
}}

renderFig(defaultFig);

function updateDotplot(col) {{
  if (col === '__default__') {{
    renderFig(defaultFig);
  }} else {{
    renderFig(sortFigures[col]);
  }}
}}
</script>
"""

# ----------------------------
# Build HTML with subtabs
# ----------------------------
html = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>scIB benchmark</title>
<style>
body {{ font-family: Arial; margin:0; display:flex; height:100vh; overflow:hidden; }}
.tab {{ display:flex; flex-direction:column; width:200px; background-color:#f1f1f1; padding-top:20px; border-right:1px solid #ccc; }}
.tab button {{ background:none; border:none; padding:12px 20px; text-align:left; cursor:pointer; font-size:16px; transition:0.3s; }}
.tab button:hover {{ background-color:#ddd; }}
.tab button.active {{ background-color:#007acc; color:white; font-weight:bold; border-left:4px solid #005f99; }}
.tabcontent {{ flex-grow:1; display:none; padding:20px; overflow:auto; }}
.subtab {{ display:flex; flex-direction:column; width:200px; border-right:1px solid #ccc; padding-top:10px; float:left; }}
.subtab button {{ background:none; border:none; padding:8px 16px; cursor:pointer; text-align:left; width:100%; font-size:14px; transition:0.2s; }}
.subtab button.active {{ background-color:#007acc; color:white; font-weight:bold; border-left:4px solid #005f99; }}
.subtabcontent {{ display:none; margin-left:210px; padding:10px; }}
.scatter-container {{ display:flex; flex-direction:row; gap:20px; }}
</style>
<script>
function openTab(evt, tabName){{
    document.querySelectorAll(".tabcontent").forEach(c=>c.style.display="none");
    document.querySelectorAll(".tab button").forEach(b=>b.classList.remove("active"));
    document.getElementById(tabName).style.display="block";
    evt.currentTarget.classList.add("active");
}}
function activateSubTab(evt, subTabName){{
    const container = evt.currentTarget.closest('.tabcontent');
    container.querySelectorAll(".subtabcontent").forEach(c=>c.style.display="none");
    container.querySelectorAll(".subtab button").forEach(b=>b.classList.remove("active"));
    container.querySelector('#'+subTabName).style.display="block";
    evt.currentTarget.classList.add("active");
}}
</script>
</head>
<body>

<div class="tab">
    <button class="tablinks active" onclick="openTab(event,'Benchmark')">Benchmark</button>
</div>

<div id="Benchmark" class="tabcontent" style="display:block;">

    <div class="subtab">
        <button class="active" onclick="activateSubTab(event,'Overview')">Overview</button>
        <button onclick="activateSubTab(event,'StrictComparison')">Strict comparison</button>
        <button onclick="activateSubTab(event,'LenientComparison')">Lenient comparison</button>
        <button onclick="activateSubTab(event,'MetricsRecap')">Metrics recap</button>
    </div>

    <!-- Overview -->
    <div id="Overview" class="subtabcontent" style="display:block;">
        <h3>Overview: Metric Dotplot</h3>
        {overview_dotplot_html}
    </div>

    <!-- Strict comparison -->
    <div id="StrictComparison" class="subtabcontent">
        <div class="scatter-container">{scatter_strict_html}</div>
    </div>

    <!-- Lenient comparison -->
    <div id="LenientComparison" class="subtabcontent">
        <div class="scatter-container">{scatter_lenient_html}</div>
    </div>

    <!-- Metrics recap -->
    <div id="MetricsRecap" class="subtabcontent">
        <p>
            <strong style="color:#0080ff;">Batch correction (label-free) metrics:</strong><br>
            • <strong>ilisi_graph</strong>: measures the local mixing of batches by counting how evenly batches are distributed among neighbors in the integrated graph; higher values indicate better mixing across batches.<br>
            • <strong>kBET</strong>: tests whether batch labels are uniformly distributed among nearest neighbors compared to random expectation; a higher score means stronger evidence of successful batch mixing.<br>
            • <strong>silhouette_batch</strong>: computes the silhouette width based on batch labels in the integrated space, where larger values indicate that cells are less separable by batch and thus better corrected.<br>
            • <strong>pcr_comparison</strong>: assesses how much batch technical variance is removed by comparing principal components before and after integration; higher scores suggest more batch effect reduction.<br>
            • <strong>graph_connectivity</strong>: evaluates the connectivity of same‑label cells in a neighborhood graph after integration, where a higher value indicates that cells with the same biological label are better connected regardless of batch.<br><br>

            <strong style="color:#0055ff;">Batch correction (with label) metrics:</strong><br>
            • <strong>clisi_graph</strong>: measures local batch separability in labeled cells; higher values indicate that batch mixing is more uniform while preserving labels.<br><br>

            <strong style="color:#00b300;">Biology (label-free) metrics:</strong><br>
            • <strong>hvg_overlap</strong>: quantifies the overlap of highly variable genes before and after integration, with higher overlap indicating better preservation of biological gene variance.<br>
            • <strong>cell_cycle</strong>: measures how well cell‑cycle associated variation is conserved post‑integration by comparing contributions to principal components, where higher values indicate stronger conservation of cell cycle signals.<br>
            • <strong>silhouette</strong>: computes the silhouette width based on biological group labels in the integrated embedding, with higher scores indicating tighter grouping of cells by biology rather than by noise or batch differences.<br><br>

            <strong style="color:#008040;">Biology (with label) metrics:</strong><br>
            • <strong>isolated_labels_asw</strong>: calculates the average silhouette width for rare (isolated) biological labels to evaluate how well they remain distinct after integration; higher values mean better label separation.<br>
            • <strong>isolated_labels_f1</strong>: computes an F1‑like score for isolated biological labels, reflecting how accurately these rare labels are clustered relative to ground truth; a higher score means better preservation of rare cell types.<br>
            • <strong>ari</strong>: the Adjusted Rand Index compares clustering in the integrated data to known biological labels, with higher values showing stronger agreement and thus better biological conservation.<br>
            • <strong>nmi</strong>: Normalized Mutual Information measures the shared information between clusters and true biological labels; higher values indicate better matching of integrated clustering to known biology.<br>
            • <strong>trajectory_conservation</strong>: assesses how well a pseudotemporal ordering (trajectory) of cells is preserved after integration by correlating pre‑ and post‑integration pseudotime; higher values mean better retention of developmental or progression signals.<br><br>

            <strong style="color:#ff7f00;">Lenient summary metrics:</strong><br>
            • <strong>Batch correction (lenient)</strong>: Mean of batch correction metric scores <strong>computed in at least one methods</strong>. Higher values indicate stronger overall batch effect removal.<br>
            • <strong>Biological conservation (lenient)</strong>: Mean of biological conservation metric scores <strong>computed in at least one methods</strong>. Higher values indicate stronger overall batch effect removal.<br>
            • <strong>Overall (lenient)</strong>: combined lenient aggregate of batch and biological metrics representing overall integration quality; higher indicates superior integration performance.<br><br>

            <strong style="color:#ff0000;">Strict summary metrics:</strong><br>
            • <strong>Batch correction (strict)</strong>: Mean of batch correction metric scores <strong>computed in all methods</strong>. Higher values indicate stronger overall batch effect removal.<br>
            • <strong>Biological conservation (strict)</strong>: Mean of biological conservation metric scores <strong>computed in all methods</strong>. Higher values indicate stronger overall batch effect removal.<br>
            • <strong>Overall (strict)</strong>: combined strict aggregate of batch and biological metrics representing overall integration quality; higher indicates superior integration performance.<br><br>
        </p>
    </div>

</div>
</body>
</html>
"""

args.output.write_text(html)
print(f"[INFO] Benchmark report written to {args.output}")
