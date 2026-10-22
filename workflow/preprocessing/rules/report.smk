rule report_per_integration:
    input:
        Path(mcfg.out_dir) / "{filename}.{feature_space}_{integration}.zarr",
    output:
        Path(mcfg.reports_dir) / "{filename}.preprocessing.{feature_space}_{integration}.html"
    params:
        norm_method = lambda w, mcfg=mcfg: mcfg.config['normalize']['method'],
        feature_method = mcfg.config['feature_selection']['method'],
        pca_n_comps = mcfg.config['pca']['n_comps'],
        batch_key = mcfg.config['batch_key'],
        label_key = mcfg.config['integration']['label_key'] if mcfg.config['integration']['label_key'] else None,
        obs_keys = lambda w: (
            " ".join(mcfg.config['projection']['extra_keys'])
            if mcfg.config['projection']['extra_keys']
            else "None"
        ),
        embed_key = lambda w, mcfg=mcfg: mcfg.get_embedding_key(w.integration),
    conda:
        Path(mcfg.env_paths['reports']),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.report.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.report.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/report_per_integration.py \
            -i {input} \
            -o {output} \
            -n {params.norm_method} \
            -f {params.feature_method} \
            -p {params.pca_n_comps} \
            --integration {wildcards.integration} \
            --featspace {wildcards.feature_space} \
            -b {params.batch_key} \
            -l {params.label_key} \
            -k {params.obs_keys} \
            -e {params.embed_key} \
            > {log} 2>&1
        """


rule report_integration_benchmark:
    input:
        Path(mcfg.out_dir) / "{filename}.scib_metrics.csv",
    output:
        Path(mcfg.reports_dir) / "{filename}.preprocessing.scib_integration_benchmark.html",
    conda:
        Path(mcfg.env_paths["reports"])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.report_benchmark.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.report_benchmark.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/report_integration_benchmark.py \
            -i {input} \
            -o {output} \
            > {log} 2>&1 
        """