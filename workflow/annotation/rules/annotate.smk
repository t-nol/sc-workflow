rule annotate_with_celltypist:
    input:
        query = rules.convert_zarr_to_h5ad.output,
        models = expand(
            Path(mcfg.temp_dir) / mcfg.config['atlas']['name'] / 'model.{cell_type_column}.pkl',
            cell_type_column=mcfg.config['atlas']['cell_type_columns']
        )
    output:
        h5ad = Path(mcfg.out_dir) / "{filename}.celltypist_annotated.h5ad",
        **(
            {"dotplot": Path(mcfg.reports_dir) / "{filename}.celltypist_dotplots.pdf"}
            if mcfg.config["celltypist"]["dotplot"].get("enable", False)
            else {}
        )
    params:
        norm_layer = mcfg.config["query"]["norm_layer"],
        gene_name_column = mcfg.config["query"]["gene_name_column"],
        majority_voting = "--majority-voting" if mcfg.config["celltypist"]["annotate"]["majority_voting"] else "",
        cluster_key = mcfg.config["celltypist"]["annotate"]["cluster_key"] or "None",
        dotplot = "--dotplot  --dotplot-output  {output.dotplot}" if mcfg.config["celltypist"]["dotplot"].get("enable", False) else "",
        use_as_reference = mcfg.config["celltypist"]["dotplot"]["use_as_reference"],
        use_as_prediction = "majority_voting" if mcfg.config["celltypist"]["annotate"]["majority_voting"] else "predicted_labels",
        filter_prediction = mcfg.config["celltypist"]["dotplot"].get("filter_prediction", 0.0),
        swap_axes = "--swap-axes" if mcfg.config["celltypist"]["dotplot"].get("swap_axes", False) else "",
        title = mcfg.config["celltypist"]["dotplot"].get("title", "CellTypist label transfer")
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.annotate_celltypist.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.annotate_celltypist.txt"
    shell:
        """
        python3 workflow/annotation/script/annotate.py \
            -i {input.query} \
            -m {input.models} \
            -o {output.h5ad} \
            --norm-layer {params.norm_layer} \
            --gene-name-column {params.gene_name_column} \
            {params.majority_voting} \
            --cluster_key {params.cluster_key} \
            {params.dotplot} \
            --use-as-reference {params.use_as_reference} \
            --use-as-prediction {params.use_as_prediction} \
            --filter-prediction {params.filter_prediction} \
            {params.swap_axes} \
            --title "{params.title}" \
            > {log} 2>&1
        """
