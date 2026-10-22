import argparse
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse
import scanpy as sc
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import rgb2hex
import io
import base64
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
import json

# ----------------------------
# Argument parsing
# ----------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, type=Path)
parser.add_argument("-o", required=True, type=Path)
parser.add_argument("-n", required=True, type=str, help="Normalization method")
parser.add_argument("-f", required=True, type=str, help="Feature selection method")
parser.add_argument("-p", required=True, type=int, help="Number of PCA components")
parser.add_argument("--integration", required=True, type=str, help="Integration method")
parser.add_argument("--featspace", required=True, type=str, help="Feature space")
parser.add_argument("-b", required=False, default=None, type=str, help="Batch key")
parser.add_argument("-l", required=False, default=None, type=str, help="Label key")
parser.add_argument("-k", required=False, nargs="+", default=None, help="Additional obs keys")
parser.add_argument("-e", required=False, type=str, default="X_pca", help="X_method key for embedding space (pca, integration...)")
args = parser.parse_args()

# ----------------------------
# Load AnnData
# ----------------------------
adata = ad.read_zarr(args.i)

# ----------------------------
# Helper functions
# ----------------------------
def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")

# ----------------------------
# Cache all gene expressions at once
# ----------------------------
keys_to_plot = [k for k in [args.b, args.l] + (args.k or []) if k]
all_gene_keys = set()

for key in keys_to_plot:
    if key in adata.var_names:
        all_gene_keys.add(key)
    elif key in adata.var['gene_name'].values:
        idxs = np.where(adata.var['gene_name'].values == key)[0]
        all_gene_keys.update(adata.var_names[idxs])

# Convert sparse to dense once
if issparse(adata.X):
    X_dense = adata.X.toarray()
else:
    X_dense = adata.X

# Build cache
gene_expr_cache = {gene: X_dense[:, adata.var_names.get_loc(gene)] for gene in all_gene_keys}

def get_expr(key):
    """Return expression values from cache or obs"""
    if key in adata.obs.columns:
        return adata.obs[key].values
    if key in adata.var_names:
        return gene_expr_cache[key]
    if key in adata.var['gene_name'].values:
        idxs = np.where(adata.var['gene_name'].values == key)[0]
        return {adata.var_names[i]: gene_expr_cache[adata.var_names[i]] for i in idxs}
    return None

# ----------------------------
# 2D PCA/UMAP plotting
# ----------------------------
def make_2d_plot(key, basis, images_list, size_scale=500_000):
    if key in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(8,8))
        sc.pl.embedding(adata, basis=basis, color=[key], ax=ax, size=size_scale/adata.n_obs, show=False)
        images_list.append((key, fig_to_base64(fig)))
    else:
        expr_dict = get_expr(key)
        if expr_dict is None:
            return
        if isinstance(expr_dict, dict):  # multiple genes with same gene_name
            n = len(expr_dict)
            n_cols = 3
            n_rows = int(np.ceil(n / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*8, n_rows*8))
            axes = axes.flatten()
            for i, (var_id, expr) in enumerate(expr_dict.items()):
                sc.pl.embedding(adata, basis=basis, color=[var_id], cmap=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True),
                                ax=axes[i], size=size_scale/adata.n_obs, 
                                title = f"{adata.var.loc[var_id,'gene_name']} : {var_id}" if var_id in adata.var_names else var_id,
                                show=False)
            for ax in axes[i+1:]:
                fig.delaxes(ax)
            images_list.append((key, fig_to_base64(fig)))
        else:  # single var
            var_id = list(expr_dict.keys())[0] if isinstance(expr_dict, dict) else key
            fig, ax = plt.subplots(figsize=(8,8))
            sc.pl.embedding(adata, basis=basis, color=[var_id], cmap=sns.cubehelix_palette(dark=0, light=.9, as_cmap=True),
                            ax=ax, size=size_scale/adata.n_obs, 
                            title = f"{adata.var.loc[var_id,'gene_name']} : {var_id}" if var_id in adata.var_names else var_id,
                            show=False)
            title = f"{adata.var.loc[var_id,'gene_name']} : {var_id}" if var_id in adata.var_names else var_id
            images_list.append((title, fig_to_base64(fig)))

# ----------------------------
# 3D UMAP plotting
# ----------------------------
def get_color_map(key):
    color_key = f"{key}_colors"
    if color_key in adata.uns and key in adata.obs:
        series = adata.obs[key]
        categories = series.cat.categories if hasattr(series, "cat") else series.unique()
        return {cat: adata.uns[color_key][i] for i, cat in enumerate(categories)}
    return None

def make_3d_umap(df, key):
    color_map = get_color_map(key)
    if color_map:
        fig = px.scatter_3d(df, x="UMAP1", y="UMAP2", z="UMAP3", color=key,
                            color_discrete_map=color_map, opacity=1.0, size=[8]*len(df), size_max=8,
                            title=f"3D UMAP colored by {key}")
    else:
        vals = df[key].values
        cubehelix = [rgb2hex(c) for c in sns.cubehelix_palette(dark=0, light=.9, as_cmap=False)]
        colorscale = [[i/(len(cubehelix)-1), cubehelix[i]] for i in range(len(cubehelix))]
        fig = px.scatter_3d(df, x="UMAP1", y="UMAP2", z="UMAP3", color=key,
                            color_continuous_scale=colorscale,
                            range_color=[vals.min(), vals.max()],
                            opacity=1.0, size=[8]*len(df), size_max=8,
                            title=f"3D UMAP colored by {key}")
    fig.update_layout(margin=dict(l=0,r=0,t=40,b=0))
    return fig

def try_add_plot(df, key, plots):
    expr = get_expr(key)
    if expr is None:
        return
    if isinstance(expr, dict):  # multiple gene_ids
        for var_id, val in expr.items():
            col_name = f"{key} : {var_id}"
            df[col_name] = val
            plots.append((col_name, make_3d_umap(df, col_name)))
    else:
        df[key] = expr
        plots.append((key, make_3d_umap(df, key)))

# ----------------------------
# Normalization figures
# ----------------------------
fig_total = px.histogram(adata.obs, x="total_counts", nbins=100, title="Total counts")
fig_total.update_traces(marker=dict(color="rgba(70,130,180,0.5)", line=dict(color="rgba(40,40,40,0.9)", width=1)))

norm_layer = adata.layers.get("norm_counts")
norm_sums = np.array(norm_layer.sum(axis=1)).ravel() if issparse(norm_layer) else norm_layer.sum(axis=1)
counts, edges = np.histogram(norm_sums, bins=100)
centers = (edges[:-1]+edges[1:])/2
fig_layer = px.bar(x=centers, y=counts)
fig_layer.update_traces(marker=dict(color="rgba(70,130,180,0.5)", line=dict(color="rgba(40,40,40,0.9)", width=1)))
fig_layer.update_layout(title=f"Normalized counts ({args.n})", xaxis_title="Sum per cell", yaxis_title="Number of cells", bargap=0, plot_bgcolor="white")

fig_norm = make_subplots(rows=1, cols=2, subplot_titles=("Total counts", f"Normalized counts ({args.n})"))
for trace in fig_total.data: fig_norm.add_trace(trace, row=1, col=1)
for trace in fig_layer.data: fig_norm.add_trace(trace, row=1, col=2)
fig_norm.update_layout(width=1500, height=700)
norm_html = f'<div style="text-align:center;">{pio.to_html(fig_norm, include_plotlyjs=True, full_html=False)}</div>'

# ----------------------------
# Feature selection
# ----------------------------
color_col = "highly_variable" if args.f=="highly_variable_genes" else "highly_deviant"
fig_fs = px.scatter(adata.var, x="means", y="dispersions", color=color_col,
                    color_discrete_sequence=['darkgrey','darkorange'],
                    hover_data=["gene_name","gene_id","means","dispersions", color_col],
                    width=1500, height=700)
fig_fs.update_layout(margin=dict(l=0,r=0,t=50,b=0))
feature_html = f'<div style="text-align:center;">{pio.to_html(fig_fs, include_plotlyjs=False, full_html=False)}</div>'

# ----------------------------
# PCA & UMAP 2D QC + additional obs keys
# ----------------------------
def build_2d_images(basis, qc_keys, extra_keys):
    qc_keys_exist = [k for k in qc_keys if k in adata.obs.columns]
    n_qc = len(qc_keys_exist)
    n_cols = 3
    n_rows = int(np.ceil(n_qc / n_cols)) if n_qc > 0 else 1
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*3, n_rows*3))
    if isinstance(axes, np.ndarray):
        axes_flat = axes.flatten()
    else:
        axes_flat = [axes]
    for i, key in enumerate(qc_keys_exist):
        sc.pl.embedding(adata, basis=basis, color=key, ax=axes_flat[i], show=False)
    for ax in axes_flat[n_qc:]:
        fig.delaxes(ax)
    qc_img = fig_to_base64(fig)
    images = [("QC Metrics", qc_img)]
    for key in extra_keys:
        make_2d_plot(key, basis, images)
    return images

qc_keys = ["log1p_n_genes_by_counts","log1p_total_counts","is_doublet_all",
           "pct_counts_mt","pct_counts_ribo","pct_counts_hb",
           "S_score","G2M_score","phase"]
extra_keys = [k for k in [args.b, args.l]+(args.k or []) if k]
pca2d_images = build_2d_images(args.e, qc_keys, extra_keys)
umap2d_images = build_2d_images("X_umap_2D", qc_keys, extra_keys)
if 'X_diffmap' in adata.obsm.keys():
    diffmap2d_images = build_2d_images("X_diffmap", qc_keys, extra_keys)


def build_2d_html(images, prefix):
    subtab_buttons_html = "".join([f'<button class="subtablinks" onclick="activateSubTab(this); show{prefix}(\'{key}\')">{key}</button>' for key,_ in images])
    subtab_contents_html = "".join([f'<div id="{prefix.lower()}-{key}" class="subtabcontent" style="display:none;"><img src="data:image/png;base64,{img}" style="max-width:100%;"></div>' for key,img in images])
    return f"""
    <div style="display:flex; height:100%;">
        <div class="subtab" style="width:200px; border-right:1px solid #ccc; padding-top:10px;">
            {subtab_buttons_html}
        </div>
        <div id="{prefix.lower()}-container" style="flex:1; height:100%; padding:10px;">
            {subtab_contents_html}
        </div>
    </div>
    <script>
    function show{prefix}(key){{
        document.querySelectorAll(".subtabcontent").forEach(c=>c.style.display="none");
        document.getElementById("{prefix.lower()}-"+key).style.display="block";
        document.querySelectorAll(".subtablinks").forEach(b=>b.classList.remove("active"));
        event.currentTarget.classList.add("active");
    }}
    window.onload = function(){{
        const firstKey = "{images[0][0]}";
        document.getElementById("{prefix.lower()}-"+firstKey).style.display="block";
        document.querySelectorAll(".subtablinks")[0].classList.add("active");
    }};
    </script>
    """

pca2d_html = build_2d_html(pca2d_images, "PCA2D")
umap2d_html = build_2d_html(umap2d_images, "UMAP2D")
if 'X_diffmap' in adata.obsm.keys():
    diffmap2d_html = build_2d_html(diffmap2d_images, "DIFFMAP2D")

# ----------------------------
# UMAP 3D interactive
# ----------------------------
df_umap3d = pd.DataFrame(adata.obsm["X_umap_3D"], columns=["UMAP1","UMAP2","UMAP3"])
umap3d_plots = []
for key in extra_keys:
    try_add_plot(df_umap3d, key, umap3d_plots)
umap3d_json = {key: fig.to_json() for key, fig in umap3d_plots}

umap3d_html = f"""
<div style="display:flex; height:100%;">
    <div class="subtab" style="width:200px; border-right:1px solid #ccc; padding-top:10px;">
        {"".join([f'<button class="subtablinks" onclick="activateSubTab(this); showUMAP3D(\'{key}\')">{key}</button>' for key, _ in umap3d_plots])}
    </div>
    <div id="umap3d-container" style="flex:1; height:100%;">
        <div id="umap3d-plot" style="width:100%; height:100%;"></div>
    </div>
</div>
<script src="https://cdn.plot.ly/plotly-2.30.0.min.js"></script>
<script>
const UMAP_PLOTS = {json.dumps(umap3d_json)};
let currentCamera = null;
function showUMAP3D(key){{
    const fig = JSON.parse(UMAP_PLOTS[key]);
    if(currentCamera && fig.layout.scene){{ fig.layout.scene.camera = currentCamera; }}
    Plotly.react("umap3d-plot", fig.data, {{...fig.layout, autosize:true, margin:{{l:0,r:0,t:40,b:0}}}}, {{responsive:true}});
    const plotDiv = document.getElementById("umap3d-plot");
    plotDiv.on('plotly_relayout', (eventdata) => {{ if(eventdata['scene.camera']) currentCamera = eventdata['scene.camera']; }});
}}
</script>
"""

# ----------------------------
# Build final HTML report
# ----------------------------
tab_buttons = []
tab_contents = []

def add_tab(title, content, active=False):
    tab_buttons.append(f'<button class="tablinks {"active" if active else ""}" onclick="openTab(event,\'{title}\')">{title}</button>')
    tab_contents.append(f'<div id="{title}" class="tabcontent" style="display:{"block" if active else "none"};">{content}</div>')

for title, content, active in [("Normalization", norm_html, True),
                               ("Feature Selection", feature_html, False),
                               ("PCA", pca2d_html, False),
                               ("UMAP 2D", umap2d_html, False),
                               ("UMAP 3D", umap3d_html, False)]:
    add_tab(title, content, active)

if "X_diffmap" in adata.obsm.keys():
    add_tab("DIFFMAP 2D", diffmap2d_html, False)

html = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>{args.integration} - {args.featspace}</title>
<style>
body {{ font-family: Arial; margin:0; display:flex; height:100vh; overflow:hidden; }}
.tab {{ display:flex; flex-direction:column; width:200px; background-color:#f1f1f1; padding-top:20px; border-right:1px solid #ccc; }}
.tab-title {{ font-weight: bold; font-size: 16px; text-align: center; padding: 10px 0; border-bottom: 1px solid #ccc; color: #333; }}
.tab button {{ background:none; border:none; padding:12px 20px; text-align:left; cursor:pointer; font-size:16px; transition:0.3s; }}
.tab button:hover {{ background-color:#ddd; }}
.tab button.active {{ background-color:#007acc; color:white; font-weight:bold; border-left:4px solid #005f99; }}
.tabcontent {{ flex-grow:1; display:none; padding:20px; overflow:auto; }}
.subtab button {{ background:none; border:none; padding:8px 16px; cursor:pointer; text-align:left; width:100%; font-size:14px; transition:0.2s; }}
.subtab button.active {{ background-color:#007acc; color:white; font-weight:bold; border-left:4px solid #005f99; }}
.subtabcontent {{ display:none; }}
</style>
<script>
function openTab(evt, tabName){{
    document.querySelectorAll(".tabcontent").forEach(c=>c.style.display="none");
    document.querySelectorAll(".tablinks").forEach(b=>b.classList.remove("active"));
    document.getElementById(tabName).style.display="block";
    evt.currentTarget.classList.add("active");
}}
function activateSubTab(btn){{
    document.querySelectorAll(".subtablinks").forEach(b=>b.classList.remove("active"));
    btn.classList.add("active");
}}
</script>
</head>
<body>
<div class="tab">
    <div class="tab-title">
        Integration: {args.integration} \n Feature space: {args.featspace}
    </div>
    {"".join(tab_buttons)}
</div>
{"".join(tab_contents)}
</body>
</html>
"""

with open(args.o, "w") as f:
    f.write(html)
