#----------------------------------------------------------
# Expression correctors (fix the expr matrix)
#----------------------------------------------------------

rule mnn:
    input:
        lambda w: mcfg.get_prepro_input("mnn", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.feature_space.mnn.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
    conda:
        Path(mcfg.env_paths["integration"]) / "mnn.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_mnn.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_mnn.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_mnn.py \
            -i {input} \
            -o {output} \
            -f {wildcards.feature_space} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """

rule combat:
    input:
        lambda w: mcfg.get_prepro_input("combat", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.feature_space.combat.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
        label_key = lambda wc: (
            f"-l {mcfg.config["integration"]['label_key']}"
            if mcfg.config['integration']['label_key'] is not None
            else ""
        ),
        covariates = lambda wc: (
            f"-c {mcfg.config["integration"]['methods'].get('combat', {}).get('covariates', [])}"
            if mcfg.config["integration"]['methods'].get('combat', {}).get('covariates', []) is not None
            else ""
        ),
    conda:
        Path(mcfg.env_paths["integration"]) / "combat.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_combat.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_combat.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_combat.py \
            -i {input} \
            -o {output} \
            -f {wildcards.feature_space} \
            -b {params.batch_key} \
            {params.covariates} \
            {params.label_key} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# Latent space generators (replaces PCA space)
#----------------------------------------------------------

rule scvi:
    input:
        lambda w: mcfg.get_prepro_input("scvi", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.scvi.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
    conda:
        Path(mcfg.env_paths["integration"]) / "scvi.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_scvi.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_scvi.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_scvi.py \
            -i {input} \
            -o {output} \
            -f {wildcards.feature_space} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """

rule scanvi:
    input:
        lambda w: mcfg.get_prepro_input("scanvi", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.scanvi.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
        label_key = mcfg.config["integration"]["label_key"],
    conda:
        Path(mcfg.env_paths["integration"]) / "scvi.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_scanvi.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_scanvi.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_scvi.py \
            -i {input} \
            -o {output} \
            -f {wildcards.feature_space} \
            -b {params.batch_key} \
            -l {params.label_key} \
            > {log} 2>&1 
        """

rule concord:
    input:
        lambda w: mcfg.get_prepro_input("concord", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.concord.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
    conda:
        Path(mcfg.env_paths["integration"]) / "concord.yaml"
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_concord.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_concord.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_concord.py \
            -i {input} \
            -o {output} \
            -b {params.batch_key} \
            -f {wildcards.feature_space} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# PCA Modifiers
#----------------------------------------------------------

rule harmony:
    input:
        lambda w: mcfg.get_prepro_input("harmony", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.harmony.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
    conda:
        Path(mcfg.env_paths["integration"]) / "harmony.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_harmony.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_harmony.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_harmony.py \
            -i {input} \
            -o {output} \
            -b {params.batch_key} \
            > {log} 2>&1 
        """

rule scanorama:
    input:
        lambda w: mcfg.get_prepro_input("scanorama", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.embedding_space.scanorama.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
        n = mcfg.config["integration"]["methods"].get("scanorama", {}).get("n_neighbors", None)
    conda:
        Path(mcfg.env_paths["integration"]) / "scanorama.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_scanorama.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_scanorama.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_scanorama.py \
            -i {input} \
            -o {output} \
            -b {params.batch_key} \
            -n {params.n} \
            > {log} 2>&1 
        """

#----------------------------------------------------------
# Neighbor Graph Modifiers
#----------------------------------------------------------

rule bbknn:
    input:
        lambda w: mcfg.get_prepro_input("bbknn", w)
    output:
        temp(directory(Path(mcfg.temp_dir) / "{filename}.{feature_space}_{integration}.graph_space.bbknn.zarr"))
    params:
        batch_key = mcfg.config["batch_key"],
        label_key = lambda wc: (
            f"-l {mcfg.config["integration"]['label_key']}"
            if mcfg.config['integration']['label_key'] is not None
            else ""
        ),
        n_comps = mcfg.config["pca"]["n_comps"],
    conda:
        Path(mcfg.env_paths["integration"]) / "bbknn.yaml",
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{feature_space}_{integration}.int_bbknn.txt"
    log:
        Path(mcfg.log_dir) / "{filename}.{feature_space}_{integration}.int_bbknn.txt"
    shell:
        """
        python3 workflow/preprocessing/scripts/integration/int_bbknn.py \
            -i {input} \
            -o {output} \
            -n {params.n_comps} \
            -b {params.batch_key} \
            {params.label_key} \
            > {log} 2>&1 
        """
