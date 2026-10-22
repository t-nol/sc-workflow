rule report:
    input:
        Path(mcfg.temp_dir) / "{filename}.doublets.zarr",
    output:
        directory(Path(mcfg.reports_dir) / "{filename}.qc"),
    params:
        batch_key = mcfg.config["batch_key"]
    conda:
        Path(mcfg.env_paths["reports"])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.report.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.report.txt"
    shell:
        """
        python3 workflow/qc/scripts/report.py \
            -u {input} \
            -o {output} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """