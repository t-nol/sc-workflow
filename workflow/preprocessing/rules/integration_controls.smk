#----------------------------------------------------------
# Expression modifiers (changes the expr matrix)
#----------------------------------------------------------

rule ctrl_shuffle:
    input:
        lambda w: mcfg.get_prepro_input("shuffle", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.feature_space.shuffle.zarr"))
    params:
        layer_key = lambda w, mcfg=mcfg: mcfg.get_layer_key(w.integration),
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.shuffle.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.shuffle.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration_controls/feature_shuffle.py \
            -i {input} \
            -o {output} \
            -l {params.layer_key} \
            > {log} 2>&1 
        """

rule ctrl_shuffle_per_batch:
    input:
        lambda w: mcfg.get_prepro_input("shuffle-per-batch", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.feature_space.shuffle-per-batch.zarr"))
    params:
        layer_key = lambda w, mcfg=mcfg: mcfg.get_layer_key(w.integration),
        shuffle_by = mcfg.config["batch_key"]
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.shuffle-per-batch.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.shuffle-per-batch.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration_controls/feature_shuffle.py \
            -i {input} \
            -o {output} \
            -l {params.layer_key} \
            -s {params.shuffle_by} \
            > {log} 2>&1 
        """

rule ctrl_shuffle_per_label:
    input:
        lambda w: mcfg.get_prepro_input("shuffle-per-label", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.feature_space.shuffle-per-label.zarr"))
    params:
        layer_key = lambda w, mcfg=mcfg: mcfg.get_layer_key(w.integration),
        shuffle_by = mcfg.config["integration"]["label_key"]
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.shuffle-per-label.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.shuffle-per-label.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration_controls/feature_shuffle.py \
            -i {input} \
            -o {output} \
            -l {params.layer_key} \
            -s {params.shuffle_by} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# Latent space generators (replaces PCA space)
#----------------------------------------------------------

rule ctrl_batchwise_pca:
    input:
        lambda w: mcfg.get_prepro_input("batchwise-pca", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.batchwise-pca.zarr"))
    params:
        n_comps = mcfg.config['pca']['n_comps'],
        rcc = mcfg.config['pca']['regress_cell_cycle'],
        layer_key = lambda w, mcfg=mcfg: mcfg.get_layer_key(w.integration),
        batch_key = mcfg.config["batch_key"]
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.batchwise-pca.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.batchwise-pca.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/preprocessing/compute_pca.py \
            -i {input} \
            -o {output} \
            -l {params.layer_key} \
            -n {params.n_comps} \
            -f {wildcards.feature_space} \
            -r {params.rcc} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """

rule ctrl_perfect_embed_label:
    input:
        lambda w: mcfg.get_prepro_input("perfect-embed-label", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.perfect-embed-label.zarr"))
    params:
        layer_key = lambda w, mcfg=mcfg: mcfg.get_layer_key(w.integration), 
        label_key = mcfg.config["integration"]['label_key']
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_ctrl_embed_label.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_ctrl_embed_label.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration_controls/perfect_label_embed.py \
            -i {input} \
            -o {output} \
            -f {wildcards.feature_space} \
            -l {params.layer_key} \
            -k {params.label_key} \
            > {log} 2>&1 
        """

rule ctrl_perfect_embed_label_jitter:
    input:
        lambda w: mcfg.get_prepro_input("perfect-embed-label-jitter", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.perfect-embed-label-jitter.zarr"))
    params:
        layer_key = lambda w, mcfg=mcfg: mcfg.get_layer_key(w.integration), 
        label_key = mcfg.config["integration"]['label_key'],
        jitter = 0.01 #1% of feature scale [0,1]
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_ctrl_embed_label_jitter.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_ctrl_embed_label_jitter.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration_controls/perfect_label_embed.py \
            -i {input} \
            -o {output} \
            -f {wildcards.feature_space} \
            -l {params.layer_key} \
            -k {params.label_key} \
            -j {params.jitter} \
            > {log} 2>&1 
        """
