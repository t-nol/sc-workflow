rule metrics:
    input:
        lambda w, mcfg=mcfg: (Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}")
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.metrics.zarr"))
    params:
        n_counts_max = mcfg.config["filtering"]["n_counts_max"],
        n_counts_min = mcfg.config["filtering"]["n_counts_min"],
        n_genes_min = mcfg.config["filtering"]["n_genes_min"],
        n_cells_min = mcfg.config["filtering"]["n_cells_min"],
        percent_mito_max = mcfg.config["filtering"]["percent_mito_max"],
        use_mad = mcfg.config["filtering"]["use_mad"],
        batch_key = mcfg.config["batch_key"],
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.metrics.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.metrics.txt"
    shell:
        """
        python3 workflow/qc/scripts/metrics.py \
            -i {input} \
            -o {output} \
            -b {params.batch_key} \
            --n_counts_max {params.n_counts_max} \
            --n_counts_min {params.n_counts_min} \
            --n_genes_min {params.n_genes_min} \
            --n_cells_min {params.n_cells_min} \
            --percent_mito_max {params.percent_mito_max} \
            --use_mad {params.use_mad} \
            > {log} 2>&1 
        """