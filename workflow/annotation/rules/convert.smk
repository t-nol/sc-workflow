rule convert_zarr_to_h5ad:
    input:
        lambda w, mcfg=mcfg: (Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}"),
    output:
        temp(Path(mcfg.temp_dir) / "{filename}.h5ad")
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.convert_zarr_to_h5ad.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.convert_zarr_to_h5ad.txt"
    shell:
        """
        python3 workflow/annotation/script/convert_zarr_to_h5ad.py \
            -i {input} \
            -o {output} \
            > {log} 2>&1
        """

