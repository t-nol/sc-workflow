from pyclustree import clustree


def create_cluster_tree(adata, levels, title):
    clustree_fig = clustree(
        adata,
        levels,
        title=title,
        edge_weight_threshold=0.00,
        show_fraction=True,
    )
    clustree_fig.set_size_inches(15, 2*len(levels))
    clustree_fig.set_dpi(100)
    return clustree_fig
