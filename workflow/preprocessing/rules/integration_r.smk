#----------------------------------------------------------
# converts .h5ad inputs to .zarr 
# for full seurat preprocessing/integration, for easier loading
#----------------------------------------------------------

rule convert_h5ad_to_zarr:
    input:
        lambda w, mcfg=mcfg: (Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}")
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.zarr"))
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.convert_h5ad_2_zarr.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.convert_h5ad_2_zarr.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/convert_h5ad_to_zarr.py \
            -i {input} \
            -o {output} \
            > {log} 2>&1 
        """


#----------------------------------------------------------
# Post R steps -> take R outputs and add to adata
# one rule for feature_space, embedding space or graph space
#----------------------------------------------------------

rule save_r_output:
    input:
        uncorrected = lambda w: mcfg.get_prepro_input(w.integration, w),
        corrected = lambda w: mcfg.get_post_r_input("embedding_space", w),
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.{integration}.zarr"))
    params:
        space = "embedding_space"
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.post_r_int.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.post_r_int.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/post_r_int.py \
            -i {input.corrected} \
            -u {input.uncorrected} \
            -o {output} \
            -s {params.space} \
            -m {wildcards.integration} \
            > {log} 2>&1 
        """


#----------------------------------------------------------
# Expression correctors (fix the expr matrix)
#----------------------------------------------------------
## limma
## scMerge


#----------------------------------------------------------
# Latent space generators (replaces PCA space)
#----------------------------------------------------------


#----------------------------------------------------------
# PCA Modifiers
#----------------------------------------------------------
## FastMNN
## Liger (iNMF)
## scAlign

rule css:
    input:
        lambda w: mcfg.get_prepro_input("css", w)
    output:
        temp(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.css.npz"),
    params:
        batch_key = mcfg.config['batch_key'],
        corr_method = mcfg.config['integration']['methods'].get('css', {}).get('correlation', []),
    conda:
        Path(mcfg.env_paths["integration"]) / "css.yaml"
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_css.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_css.txt"
    shell:
        """
        Rscript workflow/preprocessing/scripts/integration/int_css.R \
            -i {input} \
            -o {output} \
            -f {wildcards.feature_space} \
            -b {params.batch_key} \
            -c {params.corr_method} \
            > {log} 2>&1 
        """

rule cca:
    input:
        lambda w, mcfg=mcfg: (Path(mcfg.temp_dir) / f"{w.filename}.zarr")
    output:
        temp(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.cca.npz")
    params:
        batch_key = mcfg.config['batch_key'],
        n_top_genes = mcfg.config['feature_selection']['n_top_genes']
    conda:
        Path(mcfg.env_paths["integration"]) / "seurat.yaml"
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_cca.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_cca.txt"
    shell:
        """
        Rscript workflow/preprocessing/scripts/integration/int_cca.R \
            -i {input} \
            -o {output} \
            -n {params.n_top_genes} \
            -f {wildcards.feature_space} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """

rule rpca:
    input:
        lambda w, mcfg=mcfg: (Path(mcfg.temp_dir) / f"{w.filename}.zarr")
    output:
        temp(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.rpca.npz")
    params:
        batch_key = mcfg.config['batch_key'],
        n_top_genes = mcfg.config['feature_selection']['n_top_genes']
    conda:
        Path(mcfg.env_paths["integration"]) / "seurat.yaml"
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_rpca.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_rpca.txt"
    shell:
        """
        Rscript workflow/preprocessing/scripts/integration/int_rpca.R \
            -i {input} \
            -o {output} \
            -n {params.n_top_genes} \
            -f {wildcards.feature_space} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# Neighbor Graph Modifiers
#----------------------------------------------------------
## CONOS
