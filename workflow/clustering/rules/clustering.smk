rule main:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
    output:
        csv = Path(mcfg.out_dir) / "{filename}" / "main" / "{alg}_{res}_clusters.csv"
    params:
        alg = ALGO,
    conda:
        Path(mcfg.env_paths[mcfg.module_name]),
    benchmark:
        Path(mcfg.bench_dir) / '{filename}.main_clustering_{alg}_{res}_clusters.txt',
    shell:
        """
        python3 workflow/clustering/scripts/main_clustering.py \
            -i {input.zarr} \
            -o {output.csv} \
            -a {params.alg} \
            -r {wildcards.res}
        """

rule sub:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        main = Path(mcfg.out_dir) / "{filename}" / "main" / "{alg}_{res}_clusters.csv"
    output:
        csv = Path(mcfg.out_dir) / "{filename}" / "sub" / "{alg}_{res}_sub_{sub}_clusters.csv"
    params:
        alg = ALGO,
    conda:
        Path(mcfg.env_paths[mcfg.module_name]),
    benchmark:
        Path(mcfg.bench_dir) / '{filename}.sub_clustering_{alg}_{res}_sub_{sub}_clusters.txt',
    shell:
        """
        python3 workflow/clustering/scripts/sub_clustering.py \
            -i {input.zarr} \
            -a {params.alg} \
            -m {input.main} \
            -s {wildcards.sub} \
            -o {output.csv} 
        """

rule clustering:
    input:
        # main clustering markers
        markers_main = expand(
            Path(mcfg.out_dir) / "{filename}/main/{alg}_{res}_clusters.csv",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS
        ),
        # sub clustering markers (only if sub resolutions exist)
        markers_sub = expand(
            Path(mcfg.out_dir) / "{filename}/sub/{alg}_{res}_sub_{sub}_clusters.csv",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if mcfg.config['sub_res'] != None else [],
        
