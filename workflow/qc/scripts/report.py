#!/usr/bin/env python3
import sys
from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
import plotly.io as pio
import plotly.express as px
import plotly.graph_objects as go

project_root = Path(__file__).resolve().parents[2]
sys.path.append(str(project_root))
from utils.InputBuilder import build_adata

# ----------------------------
# Argument parsing
# ----------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Generate interactive QC HTML report")
    p.add_argument("-u", "--unfiltered", required=True, type=Path)
    p.add_argument("-o", "--output",     required=True, type=Path,
                   help="Output directory (created if absent)")
    p.add_argument("-b", "--batch_key",  default=None,  type=str)
    return p.parse_args()

# ----------------------------
# Helper
# ----------------------------
def sanitize_id(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", s)

# ----------------------------
# Subsampling helper
# ----------------------------
SCATTER_MAX = 3000  # max points per scatter/violin trace

def _subsample(df: pd.DataFrame, n: int = SCATTER_MAX, seed: int = 42) -> pd.DataFrame:
    """Return at most n rows, preserving the original if already small enough."""
    if len(df) <= n:
        return df
    return df.sample(n=n, random_state=seed)


# ----------------------------
# Plot functions
# ----------------------------
def plot_qc_hist(adata, metric, outlier_col=None, kept_col=None, bins=50):
    data = adata.obs[[c for c in [metric, outlier_col, kept_col] if c and c in adata.obs.columns]].copy()

    if metric not in data.columns:
        return None

    if outlier_col is None or outlier_col not in data.columns:
        data['_outlier'] = False
        outlier_col = '_outlier'
    if kept_col is None or kept_col not in data.columns:
        data['_kept'] = True
        kept_col = '_kept'

    vals = data[metric].values
    _, bin_edges = np.histogram(vals, bins=bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_width = np.diff(bin_edges)

    stacking_order = [(False, True), (True, True), (False, False), (True, False)]
    colors_fill  = {True: "tomato",  False: "skyblue"}
    opacity_kept = {True: 0.7,       False: 0.4}

    fig = go.Figure()
    bottoms = np.zeros(len(bin_centers))
    for is_outlier, is_kept in stacking_order:
        mask = (data[outlier_col] == is_outlier) & (data[kept_col] == is_kept)
        hist, _ = np.histogram(data.loc[mask, metric], bins=bin_edges)
        fig.add_trace(go.Bar(
            x=bin_centers, y=hist, width=bin_width, base=bottoms,
            marker_color=colors_fill[is_outlier], opacity=opacity_kept[is_kept],
            name=f"{'Outlier' if is_outlier else 'Not Outlier'}, {'Kept' if is_kept else 'Discarded'}"
        ))
        bottoms += hist

    fig.update_layout(barmode='stack', title=metric, xaxis_title=metric, yaxis_title="Frequency")
    return fig


def plot_qc_violin(adata, metric, outlier_col=None, kept_col=None):
    cols = [c for c in [metric, outlier_col, kept_col] if c and c in adata.obs.columns]
    data = adata.obs[cols].copy()

    if metric not in data.columns:
        return None
    if kept_col is None or kept_col not in data.columns:
        data["_kept"] = True
        kept_col = "_kept"
    if outlier_col is not None and outlier_col not in data.columns:
        outlier_col = None

    outlier_palette = {True: "tomato",     False: "skyblue"}
    kept_palette    = {True: "mediumblue", False: "red"}
    fig = go.Figure()

    # Violin: ALL cells so the distribution shape is exact
    for kept_status in [True, False]:
        full = data[data[kept_col] == kept_status]
        fig.add_trace(go.Violin(
            x=[str(kept_status)] * len(full), y=full[metric],
            name="Kept" if kept_status else "Discarded",
            line_color=kept_palette[kept_status],
            fillcolor="rgba(0,0,0,0)", width=0.8, showlegend=False, points=False,
        ))

    # Scatter overlay: subsampled for display performance
    color_col    = outlier_col if outlier_col is not None else kept_col
    palette      = outlier_palette if outlier_col is not None else kept_palette
    legend_title = "Outlier" if outlier_col is not None else "Cell status"
    for value, color in palette.items():
        subset = _subsample(data[data[color_col] == value])
        n_total = (data[color_col] == value).sum()
        fig.add_trace(go.Scatter(
            x=subset[kept_col].astype(str), y=subset[metric],
            mode="markers", marker=dict(color=color, size=3, opacity=0.5),
            name=f"{legend_title}: {value} ({len(subset):,} / {n_total:,} shown)",
        ))

    fig.update_layout(title=metric, xaxis_title="Cells kept after full QC",
                      yaxis_title=metric, violingap=0.3, violinmode="group", height=500)
    return fig


def generate_plot(sub_label, subsub, adata, batch, bins=50):
    if sub_label == "Feature Filtering Overview":
        n_unf  = adata.shape[1]
        n_filt = int(adata.var[f"qc_pass_{batch}"].sum())
        df = pd.DataFrame({"Dataset": ["Unfiltered", "Post Filtering"],
                           "Features": [n_unf, n_filt]})
        return px.bar(df, x="Dataset", y="Features", title="Number of Features Post Filtering")

    if sub_label == "Cell QC Metrics":
        if subsub == f"qc_pass_{batch}" and subsub in adata.obs.columns:
            df = adata.obs[[subsub]].copy()
            df[subsub] = df[subsub].astype(str)
            return px.histogram(df, x=subsub, color=subsub, barmode="group",
                                title="Cells Passing All Filters")
        outlier_col = f"{subsub}_outlier" if f"{subsub}_outlier" in adata.obs.columns else None
        kept_col    = "qc_pass"           if "qc_pass"           in adata.obs.columns else None
        if subsub in {"pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"}:
            return plot_qc_violin(adata, subsub, outlier_col=outlier_col, kept_col=kept_col)
        return plot_qc_hist(adata, subsub, outlier_col=outlier_col, kept_col=kept_col, bins=bins)

    if sub_label == "Doublet Detection":
        n_total = adata[adata.obs[f"qc_pass_{batch}"]].n_obs
        if subsub in ("is_doublet_all", "is_doublet_any"):
            title      = "Doublets detected by all methods" if subsub == "is_doublet_all" \
                         else "Doublets detected by at least one method"
            n_doublets = int(adata.obs[subsub].sum())
        else:
            title      = f"Doublets detected by tool: {subsub.split('_')[0]}"
            n_doublets = adata[adata.obs[subsub] == "True"].n_obs
        df = pd.DataFrame({"status": ["Non-doublet", "Doublet"],
                           "count":  [n_total - n_doublets, n_doublets]})
        return px.bar(df, x="status", y="count", text="count", title=title, color="status")

    return None


# ----------------------------
# Per-batch HTML generation
# ----------------------------
def build_batch_html(filename: str, batch: str, adata, title: str) -> str:
    """Generate a complete standalone HTML for one batch."""

    subtabs = {
        "Feature Filtering Overview": [],
        "Cell QC Metrics": [
            f"qc_pass_{batch}",
            "total_counts", "log1p_total_counts",
            "n_genes_by_counts", "log1p_n_genes_by_counts",
            "pct_counts_in_top_20_genes",
            "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb",
        ],
        "Doublet Detection": ["is_doublet_all", "is_doublet_any"] + [
            col for col in adata.obs.columns if col.endswith("_doublet_prediction")
        ],
    }

    subtab_buttons_html = ""
    subtab_contents_html = ""
    # Plotly.js included only once — on the first figure
    plotlyjs_included = False
    js_plot_scripts = ""
    first_subtab = True

    for sub_label, children in subtabs.items():
        subtab_id  = sanitize_id(sub_label)
        active_cls = "active" if first_subtab else ""
        first_subtab = False

        subtab_buttons_html += (
            f'<button class="tablinks_subsub {active_cls}" '
            f'onclick="openSubTab(event,\'{subtab_id}\')">{sub_label}</button>\n'
        )

        if children:
            subsub_buttons = subsub_contents = ""
            first_child = True
            for child in children:
                child_id     = sanitize_id(f"{subtab_id}_{child}")
                active_child = "active" if first_child else ""
                first_child  = False

                subsub_buttons += (
                    f'<button class="tablinks_subsubsub {active_child}" '
                    f'onclick="openSubSubTab(event,\'{child_id}\')">{child}</button>'
                )
                subsub_contents += f'<div id="{child_id}" class="tabcontent_subsubsub"></div>'

                fig = generate_plot(sub_label, child, adata, batch)
                if fig is not None:
                    include_js = "cdn" if not plotlyjs_included else False
                    plotlyjs_included = True
                    js_plot_scripts += pio.to_html(
                        fig, full_html=False, include_plotlyjs=include_js, div_id=child_id
                    )

            inner = (f'<div class="subsub-buttons">{subsub_buttons}</div>'
                     f'<div class="subsub-content">{subsub_contents}</div>')
        else:
            inner = f'<div id="{subtab_id}_plot" class="tabcontent_subsub_plot"></div>'
            fig = generate_plot(sub_label, sub_label, adata, batch)
            if fig is not None:
                include_js = "cdn" if not plotlyjs_included else False
                plotlyjs_included = True
                js_plot_scripts += pio.to_html(
                    fig, full_html=False, include_plotlyjs=include_js, div_id=f"{subtab_id}_plot"
                )

        subtab_contents_html += (
            f'<div id="{subtab_id}" class="tabcontent_subsub">{inner}</div>'
        )

    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>{title}</title>
  <style>
    body {{font-family: Arial; margin: 0;}}
    .topbar {{background: #333; color: white; padding: 8px 16px; display: flex;
              align-items: center; gap: 16px; font-size: 14px;}}
    .topbar a {{color: #adf; text-decoration: none; font-size: 13px;}}
    .topbar a:hover {{text-decoration: underline;}}
    .container {{display: flex; height: calc(100vh - 37px);}}
    .subtab-buttons {{flex: 0 0 160px; display: flex; flex-direction: column;
                      background: #f1f1f1; padding-top: 10px; overflow-y: auto;}}
    .tablinks_subsub {{display: block; background: inherit; border: none;
                       padding: 8px 12px; margin-bottom: 2px; cursor: pointer; text-align: left;}}
    .tablinks_subsub.active {{background: #ccc;}}
    .subtab-content {{flex: 1; overflow: auto; padding: 10px;}}
    .tabcontent_subsub {{display: none; flex-direction: row;}}
    .subsub-buttons {{flex: 0 0 180px; display: flex; flex-direction: column; margin-right: 10px;}}
    .tablinks_subsubsub {{display: block; background: #ddd; border: none;
                          padding: 5px 10px; margin-bottom: 2px; cursor: pointer; text-align: left;}}
    .tablinks_subsubsub.active {{background: #bbb;}}
    .subsub-content {{flex: 1; overflow: auto;}}
    .tabcontent_subsubsub {{display: none; flex: 1;}}
    .tabcontent_subsub_plot {{flex: 1; min-height: 400px;}}
  </style>
</head>
<body>
  <div class="topbar">
    <span>{title}</span>
    <a href="index.html">&#8592; Back to index</a>
  </div>
  <div class="container">
    <div class="subtab-buttons">{subtab_buttons_html}</div>
    <div class="subtab-content">{subtab_contents_html}</div>
  </div>

  {js_plot_scripts}

  <script>
  function openSubTab(evt, tabName) {{
    document.querySelectorAll(".tabcontent_subsub")
            .forEach(el => el.style.display = "none");
    document.querySelectorAll(".tablinks_subsub")
            .forEach(el => el.classList.remove("active"));
    document.getElementById(tabName).style.display = "flex";
    evt.currentTarget.classList.add("active");
    const first = document.getElementById(tabName).querySelector(".tablinks_subsubsub");
    if (first) first.click();
  }}
  function openSubSubTab(evt, tabName) {{
    document.querySelectorAll(".tabcontent_subsubsub")
            .forEach(el => el.style.display = "none");
    document.querySelectorAll(".tablinks_subsubsub")
            .forEach(el => el.classList.remove("active"));
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.classList.add("active");
  }}
  document.querySelector(".tablinks_subsub").click();
  </script>
</body>
</html>"""


# ----------------------------
# Index page
# ----------------------------
def build_index_html(filename: str, batch_files: list[tuple[str, str]]) -> str:
    """batch_files: list of (batch_label, filename.html)"""
    links = ""
    for label, fname in batch_files:
        links += (
            f'<li><a href="{fname}">{label}</a></li>\n'
        )
    return f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>QC Report — {filename}</title>
  <style>
    body {{font-family: Arial; margin: 40px; color: #333;}}
    h1 {{font-size: 22px; margin-bottom: 24px;}}
    ul {{list-style: none; padding: 0;}}
    li {{margin-bottom: 10px;}}
    a {{color: #1a6fb5; font-size: 15px; text-decoration: none;}}
    a:hover {{text-decoration: underline;}}
  </style>
</head>
<body>
  <h1>QC Report — {filename}</h1>
  <ul>{links}</ul>
</body>
</html>"""


# ----------------------------
# Main
# ----------------------------
def generate_report(unfiltered_file: Path, output_dir: Path, batch_key) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    adata_unf = build_adata(unfiltered_file)
    filename  = unfiltered_file.name.split(".")[0]

    batches = (
        ["global"] + adata_unf.obs[batch_key].cat.categories.tolist()
        if batch_key and batch_key != "None"
        else ["global"]
    )
    print("Batches:", batches)

    batch_files = []

    for batch in batches:
        print(f"  Building {batch}…")
        adata  = adata_unf if batch == "global" \
                 else adata_unf[adata_unf.obs[batch_key] == batch]
        label  = f"{filename} / {batch}"
        fname  = f"report_{sanitize_id(batch)}.html"
        html   = build_batch_html(filename, batch, adata, label)
        (output_dir / fname).write_text(html, encoding="utf-8")
        batch_files.append((label, fname))
        print(f"    → {fname}")

    index_html = build_index_html(filename, batch_files)
    (output_dir / "index.html").write_text(index_html, encoding="utf-8")

    print(f"\nDone. Open {output_dir / 'index.html'} to start.")


if __name__ == "__main__":
    args = parse_args()
    generate_report(args.unfiltered, args.output, args.batch_key)