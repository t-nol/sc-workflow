rule gprofiler:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        markers = Path(mcfg.out_dir) / "{filename}" / "markers" / "{level}" / "{clust_name}_markers.csv",
    output:
        enrich = Path(mcfg.out_dir) / "{filename}" / "markers" / "{level}" / "{clust_name}_enrichment.csv",
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{level}_gprofiler_{clust_name}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/run_gprofiler.py \
            -i {input.zarr} \
            -m {input.markers} \
            -o {output.enrich}
        """


rule gprofiler_subclustering:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        markers = Path(mcfg.out_dir) / "{filename}" / "markers" / "{level}" / "{clust_name}_submarkers.csv",
    output:
        enrich = Path(mcfg.out_dir) / "{filename}" / "markers" / "{level}" / "{clust_name}_subenrichment.csv",
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.{level}_gprofiler_{clust_name}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/run_gprofiler.py \
            -i {input.zarr} \
            -m {input.markers} \
            -o {output.enrich} 
        """
