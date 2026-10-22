#!/usr/bin/env bash
# css_post_deployment.sh
# Install missing R packages for CSS integration inside the conda environment

# Ensure R installs packages in the current conda environment
R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"

Rscript -e '
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

devtools::install_github("quadbiolab/simspec")
'

