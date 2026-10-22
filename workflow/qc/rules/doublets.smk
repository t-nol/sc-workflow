rule doublets_collection:
    input:
        ad = lambda w: Path(mcfg.temp_dir) / f"{w.filename}.metrics.zarr",
        csvs = lambda w: [
            Path(mcfg.temp_dir) / f"{w.filename}.doublets" / f"{tool}.csv"
            for tool in mcfg.get_doublet_tools()
        ]
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.doublets.zarr")),
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.doublets_collection.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.doublets_collection.txt"
    shell:
        """
        python3 workflow/qc/scripts/doublets_collection.py \
            -i {input.ad} \
            -o {output} \
            -c {input.csvs} \
            > {log} 2>&1 
        """

rule scrublet:
    input:
        lambda w: Path(mcfg.temp_dir) / f"{w.filename}.metrics.zarr"
    output:
        Path(mcfg.temp_dir) / "{filename}.doublets" / "scrublet.csv"
    params:
        batch_key = mcfg.config["batch_key"]
    conda:
        Path(mcfg.env_paths["scrublet"])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.doublets_scrublet.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.doublets_scrublet.txt"
    shell:
        """
        python3 workflow/qc/scripts/doublets_scrublet.py \
            -i {input} \
            -o {output} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """
