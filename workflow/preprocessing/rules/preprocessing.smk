#----------------------------------------------------------
# Expression Normalization
#----------------------------------------------------------

rule normalize:
    input:
        lambda w, mcfg=mcfg: (Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}")
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.norm.zarr"))
    params:
        method = mcfg.config['normalize']['method']
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.norm.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.norm.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/preprocessing/normalization.py \
            -i {input} \
            -o {output} \
            -m {params.method} \
            > {log} 2>&1 
        """


#----------------------------------------------------------
# Cell Phase Scoring
#----------------------------------------------------------

rule cell_phase_scoring:
    input:
        lambda w: Path(mcfg.temp_dir) / f"{w.filename}.norm.zarr"
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.phase.zarr"))
    params:
        method = mcfg.config["cell_phase_scoring"]["method"],
        gene_lists = mcfg.get_gene_lists_path(mcfg.config["cell_phase_scoring"]["method"])
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.cell_phase_scoring.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.cell_phase_scoring.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/preprocessing/cell_phase_scoring.py \
            -i {input} \
            -o {output} \
            -m {params.method} \
            -l {params.gene_lists} \
            > {log} 2>&1 
        """


#----------------------------------------------------------
# Feature Selection
#----------------------------------------------------------

rule feature_selection:
    input:
        lambda w: Path(mcfg.temp_dir) / f"{w.filename}.phase.zarr"
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.feature_selection.zarr"))
    params:
        method = mcfg.config['feature_selection']['method'],
        n_top_genes = mcfg.config['feature_selection']['n_top_genes'],
        batch_key = lambda wc: (
            f"-b {mcfg.config['batch_key']}"
            if mcfg.config['feature_selection']['batch_aware']
            else ""
        )
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.feature_selection.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.feature_selection.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/preprocessing/feature_selection.py \
            -i {input} \
            -o {output} \
            -m {params.method} \
            -n {params.n_top_genes} \
            {params.batch_key} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# PCA
#----------------------------------------------------------

rule pca:
    input:
        lambda w: mcfg.get_prepro_input("pca", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.pca.zarr"))
    params:
        n_comps = mcfg.config['pca']['n_comps'],
        rcc = mcfg.config['pca']['regress_cell_cycle'],
        layer_key = lambda w, mcfg=mcfg: mcfg.get_layer_key(w.integration) 
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.pca.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.pca.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/preprocessing/compute_pca.py \
            -i {input} \
            -o {output} \
            -l {params.layer_key} \
            -n {params.n_comps} \
            -f {wildcards.feature_space} \
            -r {params.rcc} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# Neighbors
#----------------------------------------------------------

rule neighbors:
    input:
        lambda w: mcfg.get_prepro_input("neighbors", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.graph_space.neighbors.zarr"))
    params:
        knn_neighbors = mcfg.config['knn']['n_neighbors'],
        embed_key = lambda w, mcfg=mcfg: mcfg.get_embedding_key(w.integration) 
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.neighbors.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.neighbors.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/preprocessing/compute_neighbors.py \
            -i {input} \
            -o {output} \
            -n {params.knn_neighbors} \
            -e {params.embed_key} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# Dim Reduction
#----------------------------------------------------------

rule projection:
    input:
        lambda w: mcfg.get_prepro_input("projection", w)
    output:
        directory(Path(mcfg.out_dir) / "{filename}.{feature_space}_{integration}.zarr")
    params:
        methods = mcfg.config['projection']['methods'],
        do_3d = mcfg.config['projection']['do_3d']
    conda:
        Path(mcfg.env_paths[mcfg.module_name])
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.projection.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.projection.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/preprocessing/compute_projection.py \
            -i {input} \
            -o {output} \
            -m {params.methods} \
            -d {params.do_3d} \
            > {log} 2>&1 
        """