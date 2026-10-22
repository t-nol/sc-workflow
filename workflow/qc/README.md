# QC MODULE

## Goal

Quality control (QC) on a single-cell AnnData object, optionally stratified by batch (as defined by `batch_key` in config). 

## Steps

- **metrics computing** : Computes standard QC metrics (library size, gene counts, mitochondrial/ribosomal content) and filters genes expressed in too few cells. Cells are then flagged as low quality either using robust MAD-based outlier detection or fixed user-defined thresholds in config. QC is applied per batch and the results are merged, retaining cells and genes that pass QC in at least one batch in `qc_pass_global` columns and batch specific information in `qc_pass_{batch}` columns.

- **doublet prediction** : Computes `{tool}_doublet_score` and `{tool}_doublet_prediction` for each available doublet detection tool specified in config. Two new observations are created : `is_doublet_any` tracks wether a cell is defined as a doublet by any (at least one) of the tools used and `is_doublet_all` tracks wether a cell is defined as a doublet by all tools used.

- **filtering** : Filter cells and features in full object based on global qc results (obs and var `qc_pass_global` columns) as well as based on doublet predictions results and `doublets/filter_out` params (null by default, can be any or all).

- **report** : Generates an interactive HTML report with tabs per batch showing cells and features passing each QC step.

## Config parameters

Default configuration values are defined in the `config.defaults.yaml` file. To override any default setting, specify the desired values under the qc: section of your config.yaml file.

- **filtering** :
    - *n_counts_max*: Maximum total count per cell for a cell to be retained.
    - *n_counts_min*: Minimum total count per cell for a cell to be retained.
    - *n_genes_min*: Minimum number of expressed features required for a cell to be retained.
    - *n_cells_min*: Minimum number of cells in which a feature must be expressed for that feature to be retained.
    - *percent_mito_max*: Maximum allowed percentage of mitochondrial counts per cell. Note that this parameter will still be used if `use_mad: True` to make sure to filter bad quality cells in datasets that would have a high percentage of mitochondrial counts.
    - *use_mad*: *Optional* Whether to use the median absolute deviation (MAD) to identify outliers in standard quality control metrics. If set to `True`, MAD-based thresholds are used and all other filtering parameters are ignored

- **doublets** :
    - *tool*: List of tools to use for doublet detection. Possible values include : scrublet.
    - *filter_out*: How to filter out identified doublets. Use `null` to keep all doublets, `any` to remove doublets identified by at least one of the listed tools or `all` to only remove cells marked as doublets by all tools used.
    null #null to keep all cells; 'all' or 'any' to either remove cells markeds as doublet by all tools used, or by at least one (any)

- **batch_key** : *Optional* AnnData.obs key used to split the dataset by batch before quality control. Lease as `null` if not used.
