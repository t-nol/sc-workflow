rule plot_celltypist_umap:
    input:
        h5ad = rules.annotate_with_celltypist.output.h5ad,
    output:
        plot = Path(mcfg.reports_dir) / "{filename}.celltypist_umap.pdf",
    params:
        use_existing_umap = "--use-existing-umap" if mcfg.config["query"]["use_existing_umap"] else "",
        umap_key = mcfg.config["query"]["umap_key"],
    conda:
        Path(mcfg.env_paths[mcfg.module_name]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.plot_celltypist_umap.txt",
    shell:
        """
        python3 workflow/annotation/script/plot.py \
            -i {input.h5ad} \
            -p {output.plot} \
            {params.use_existing_umap} \
            --umap-key {params.umap_key}
        """
