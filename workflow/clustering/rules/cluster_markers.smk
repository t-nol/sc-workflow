rule get_cluster_markers:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        clusters = Path(mcfg.out_dir) / "{filename}" / "{level}" / "{clust_name}_clusters.csv"
    output:
        csv = Path(mcfg.out_dir) / "{filename}" / "markers" / "{level}" / "{clust_name}_markers.csv",
        dendro = Path(mcfg.out_dir) / "{filename}" / "markers"  / "{level}" / "{clust_name}_dendrogram.pkl"
    params:
        clust_key = "{clust_name}",
    conda:
        Path(mcfg.env_paths[mcfg.module_name]),
    benchmark:
        Path(mcfg.bench_dir) / '{filename}.{level}_cluster_markers_{clust_name}.txt',
    shell:
        """
        python3 workflow/clustering/scripts/cluster_markers.py \
            -i {input.zarr} \
            -c {input.clusters} \
            -o {output} \
            -k {params.clust_key}
        """

rule get_subcluster_markers:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        clusters = Path(mcfg.out_dir) / "{filename}" / "{level}" / "{clust_name}_clusters.csv"
    output:
        csv = Path(mcfg.out_dir) / "{filename}" / "markers" / "{level}" / "{clust_name}_submarkers.csv",
        dendro = Path(mcfg.out_dir) / "{filename}" / "markers"  / "{level}" / "{clust_name}_subdendrogram.pkl"
    params:
        clust_key = "{clust_name}",
    conda:
        Path(mcfg.env_paths[mcfg.module_name]),
    benchmark:
        Path(mcfg.bench_dir) / '{filename}.{level}_subcluster_markers_{clust_name}.txt',
    shell:
        """
        python3 workflow/clustering/scripts/cluster_markers.py \
            -i {input.zarr} \
            -c {input.clusters} \
            -o {output} \
            -k {params.clust_key} \
            -s True
        """

rule markers_all:
    input:
        markers_main = expand(
            Path(mcfg.out_dir) / "{filename}/markers/main/{alg}_{res}_markers.csv",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS
        ) if mcfg.config['markers'] else "",
        
        markers_sub = expand(
            Path(mcfg.out_dir) / "{filename}/markers/sub/{alg}_{res}_sub_{sub}_markers.csv",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if (
            mcfg.config['markers']
            and mcfg.config['sub_res'] != [0]
         ) else [],
        
        submarkers_sub = expand(
            Path(mcfg.out_dir) / "{filename}/markers/sub/{alg}_{res}_sub_{sub}_submarkers.csv",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if (
            mcfg.config['markers']
            and mcfg.config['sub_res'] != [0]
         ) else [],
         