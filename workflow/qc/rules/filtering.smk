rule filtering:
    input:
        ad = lambda w: Path(mcfg.temp_dir) / f"{w.filename}.doublets.zarr"
    output:
        directory(Path(mcfg.out_dir) / "{filename}.zarr")
    params:
        db_filter = mcfg.config["doublets"]["filter_out"]
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.filtering.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.filtering.txt"
    shell:
        """
        python3 workflow/qc/scripts/filtering.py \
            -i {input.ad} \
            -o {output} \
            -f {params.db_filter} \
            > {log} 2>&1 
        """