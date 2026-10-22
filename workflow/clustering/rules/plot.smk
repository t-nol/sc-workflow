# ----------------------
# Plot main summary
# ----------------------
rule plot_main_summary:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        main=lambda wc: expand(
            Path(mcfg.out_dir) / f"{wc.filename}/main/{{alg}}_{{res}}_clusters.csv",
            alg=ALGO,
            res=MAIN_RESOLUTIONS
        )
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/main/summary.pdf",
    params:
        alg = ALGO,
        main_res = MAIN_RESOLUTIONS,
    conda:
        Path(mcfg.env_paths["reports"])
    shell:
        """
        python3 workflow/clustering/scripts/plot_summary.py \
            -i {input.zarr} \
            -m {input.main} \
            -a {params.alg} \
            --main-res {params.main_res} \
            -o {output.pdf}
        """

# ----------------------
# Plot sub summary
# ----------------------
rule plot_sub_summary:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        main = Path(mcfg.out_dir) / "{filename}/main/{alg}_{res}_clusters.csv",
        subs = lambda wc: expand(
            Path(mcfg.out_dir) / f"{wc.filename}/sub/{wc.alg}_{wc.res}_sub_{{sub}}_clusters.csv",
            sub=SUB_RESOLUTIONS
        )
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/sub/summary_{alg}_{res}.pdf",
    params:
        alg = ALGO,
        main_res = "{res}",
        sub_res = SUB_RESOLUTIONS,
    conda:
        Path(mcfg.env_paths["reports"])
    shell:
        """
        python3 workflow/clustering/scripts/plot_summary.py \
            -i {input.zarr} \
            -m {input.main} \
            -s {input.subs} \
            -a {params.alg} \
            --main-res {params.main_res} \
            --sub-res {params.sub_res} \
            -o {output.pdf}
        """

# ----------------------
# Plot DEG
# ----------------------
rule plot_deg:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        clusters = Path(mcfg.out_dir) / "{filename}/{level}/{clust_name}_clusters.csv",
        markers_csv = Path(mcfg.out_dir) / "{filename}/markers/{level}/{clust_name}_markers.csv",
        dendro = Path(mcfg.out_dir) / "{filename}/markers/{level}/{clust_name}_dendrogram.pkl",
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/{level}/{clust_name}_markers.pdf",
    params:
        clust_key = "{clust_name}",
        top_n = mcfg.config['top_n'],
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.plot_markers_{level}_{clust_name}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/plot_markers.py \
            -i {input.zarr} \
            -c {input.clusters} \
            -m {input.markers_csv} \
            -d {input.dendro} \
            -k {params.clust_key} \
            -n {params.top_n} \
            -o {output.pdf}
        """

rule plot_subdeg:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        clusters = Path(mcfg.out_dir) / "{filename}/{level}/{clust_name}_clusters.csv",
        markers_csv = Path(mcfg.out_dir) / "{filename}/markers/{level}/{clust_name}_submarkers.csv",
        dendro = Path(mcfg.out_dir) / "{filename}/markers/{level}/{clust_name}_dendrogram.pkl",
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/{level}/{clust_name}_submarkers.pdf",
    params:
        clust_key = "{clust_name}",
        top_n = mcfg.config['top_n'],
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.plot_markers_{level}_{clust_name}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/plot_markers.py \
            -i {input.zarr} \
            -c {input.clusters} \
            -m {input.markers_csv} \
            -d {input.dendro} \
            -k {params.clust_key} \
            -n {params.top_n} \
            -o {output.pdf} \
            -s True
        """


# ----------------------
# Plot DEG scatter
# ----------------------
rule plot_deg_scatter:
    input:
        markers_csv = Path(mcfg.out_dir) / "{filename}/markers/{level}/{clust_name}_markers.csv",
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/{level}/{clust_name}_markers_scatter.pdf",
    params:
        clust_key = "{clust_name}",
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.plot_marker_scatter_{level}_{clust_name}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/plot_markers_scatter.py \
            -m {input.markers_csv} \
            -k {params.clust_key} \
            -o {output.pdf}
        """

# ----------------------
# Plot g:Profiler enrichment
# ----------------------
rule plot_gprofiler:
    input:
        enrich = Path(mcfg.out_dir) / "{filename}/markers/{level}/{clust_name}_enrichment.csv",
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/{level}/{clust_name}_enrichment.pdf",
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.plot_gprofiler_{level}_{clust_name}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/plot_gprofiler.py \
            -i {input.enrich} \
            -o {output.pdf}
        """

rule plot_gprofiler_subclustering:
    input:
        enrich = Path(mcfg.out_dir) / "{filename}/markers/{level}/{clust_name}_subenrichment.csv",
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/{level}/{clust_name}_subenrichment.pdf",
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.plot_gprofiler_subclustering_{level}_{clust_name}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/plot_gprofiler.py \
            -i {input.enrich} \
            -o {output.pdf}
        """


# ----------------------
# Plot user-defined markers
# ----------------------
rule plot_user_markers:
    input:
        zarr = lambda w, mcfg=mcfg: Path(mcfg.input_dir) / f"{w.filename}.{mcfg.input_type}",
        clusters = Path(mcfg.out_dir) / "{filename}/{level}/{alg}_{res}_clusters.csv",
        dendro = Path(mcfg.out_dir) / "{filename}/markers/{level}/{alg}_{res}_dendrogram.pkl",
        markers_yaml = lambda w: MARKER_STEMS[w.stem],
    output:
        pdf = Path(mcfg.reports_dir) / "{filename}.clustering/{level}/{stem}/{alg}_{res}.pdf",
    params:
        clust_key = "{alg}_{res}",
        top_n = mcfg.config['top_n'],
    conda:
        Path(mcfg.env_paths["reports"]),
    benchmark:
        Path(mcfg.bench_dir) / "{filename}.plot_{stem}_{level}_{alg}_{res}.txt",
    shell:
        """
        python3 workflow/clustering/scripts/plot_user_markers.py \
            -i {input.zarr} \
            -c {input.clusters} \
            -m {input.markers_yaml} \
            -d {input.dendro} \
            -k {params.clust_key} \
            -n {params.top_n} \
            -o {output.pdf}
        """

# ----------------------
# Aggregate all plots
# ----------------------
rule plot_all:
    input:
        main = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/main/summary.pdf",
            filename=mcfg.input_files
        ),
        sub = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/sub/summary_{alg}_{res}.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS
        ) if mcfg.config['sub_res'] is not None else [],
        
        degs_main = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/main/{alg}_{res}_markers.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS
        ) if mcfg.config['markers'] else [],
        degs_sub = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/sub/{alg}_{res}_sub_{sub}_markers.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if mcfg.config['markers'] and mcfg.config['sub_res'] is not None else [],

        subdegs_sub = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/sub/{alg}_{res}_sub_{sub}_submarkers.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if mcfg.config['markers'] and mcfg.config['sub_res'] is not None else [],
        
        scatter_main = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/main/{alg}_{res}_markers_scatter.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS
        ) if mcfg.config['markers'] else [],
        scatter_sub = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/sub/{alg}_{res}_sub_{sub}_markers_scatter.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if mcfg.config['markers'] and mcfg.config['sub_res'] is not None else [],
        
        gprofile_main = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/main/{alg}_{res}_enrichment.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS
        ) if mcfg.config['enrichment'] == 'gprofiler' else [],
        gprofile_sub = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/sub/{alg}_{res}_sub_{sub}_enrichment.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if mcfg.config['enrichment'] == 'gprofiler' and mcfg.config['sub_res'] is not None else [],

        gprofile_sub_subclustering = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/sub/{alg}_{res}_sub_{sub}_subenrichment.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS
        ) if mcfg.config['enrichment'] == 'gprofiler' and mcfg.config['sub_res'] is not None else [],
        
        user_markers_main = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/main/{stem}/{alg}_{res}.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            stem=MARKER_STEMS,
        ) if mcfg.config['user_markers'] is not None else [],
        user_markers_sub = expand(
            Path(mcfg.reports_dir) / "{filename}.clustering/sub/{stem}/{alg}_{res}_sub_{sub}.pdf",
            filename=mcfg.input_files,
            alg=ALGO,
            res=MAIN_RESOLUTIONS,
            sub=SUB_RESOLUTIONS,
            stem=MARKER_STEMS,
        ) if mcfg.config['user_markers'] is not None and mcfg.config['sub_res'] is not None else []
