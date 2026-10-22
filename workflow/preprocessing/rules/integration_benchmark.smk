rule scib_metrics:
    input:
        unintegrated = Path(mcfg.out_dir) / "{filename}.{feature_space}_noint.zarr",
        integrated = Path(mcfg.out_dir) / "{filename}.{feature_space}_{integration}.zarr",
    output:
        temp(Path(mcfg.out_dir) / "{filename}.{feature_space}_{integration}.scib_metrics.csv")
    params:
        feature_space = lambda w: w.feature_space,
        batch_key = mcfg.config["batch_key"],
        label_key = mcfg.config["integration"]["label_key"],
        embed_key = lambda w, mcfg=mcfg: mcfg.get_embedding_key(w.integration),
        integration_method = lambda w: w.integration,
        integration_product = lambda w, mcfg=mcfg: mcfg.get_integration_product(w.integration),
    conda:
        Path(mcfg.env_paths["integration"]) / "scib.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.scib_metrics.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.scib_metrics.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/scib_metrics.py \
            -u {input.unintegrated} \
            -i {input.integrated} \
            -o {output} \
            -f {params.feature_space} \
            -b {params.batch_key} \
            -l {params.label_key} \
            -e {params.embed_key} \
            -m {params.integration_method} \
            -p {params.integration_product} \
            > {log} 2>&1 
        """

rule scib_metrics_collection:
    input:
        lambda w: expand(
            Path(mcfg.out_dir) / f"{w.filename}" ".{feature_space}_{integration}.scib_metrics.csv",
            integration=[m for m in mcfg.config["integration"]["methods"]],
            feature_space=mcfg.config["integration"]["feature_space"]
        ),
    output:
        Path(mcfg.out_dir) / "{filename}.scib_metrics.csv",
    conda:
        Path(mcfg.env_paths["reports"])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.scib_metrics_collection.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.scib_metrics_collection.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/scib_metrics_collection.py \
            -i {input} \
            -o {output} \
            > {log} 2>&1 
        """