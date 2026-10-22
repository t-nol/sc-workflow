rule train_celltypist_model:
    input:
        atlas = mcfg.config['atlas']['path'],
    output:
        model = Path(mcfg.temp_dir) / mcfg.config['atlas']['name'] / 'model.{cell_type_column}.pkl'
    params:
        gene_id_column = mcfg.config['atlas']['gene_id_column'],
        top_genes = mcfg.config['celltypist']['train']['top_genes'],
        feature_selection_flag = (
            "--feature-selection"
            if mcfg.config['celltypist']['train']['feature_selection']
            else ""
        ),
        check_expression_flag = (
            "--check-expression"
            if mcfg.config['celltypist']['train']['check_expression']
            else ""
        ),
    conda:
        Path(mcfg.env_paths[mcfg.module_name]),
    threads: 10
    benchmark:
        Path(mcfg.bench_dir) / "train_celltypist_model.{cell_type_column}.txt",
    log:
        Path(mcfg.log_dir) / "train_celltypist_model.{cell_type_column}.txt",
    shell:
        """
        python3 workflow/annotation/script/train_model.py \
           -i {input.atlas} \
           -o {output.model} \
           --gene-id-column {params.gene_id_column} \
           --cell-type-column {wildcards.cell_type_column} \
           --top-genes {params.top_genes} \
           {params.feature_selection_flag} \
           {params.check_expression_flag} \
           --n-jobs {threads} \
            > {log} 2>&1
        """
