#!/usr/bin/env bash
set -e -x

pipeline="/home/theo-n/Documents/codes/sc-process"
configfile="/home/theo-n/Documents/codes/sc-process/configs/config_example.yaml"

snakemake \
  --configfile $configfile \
  --snakefile $pipeline/workflow/Snakefile \
  --use-conda \
  --printshellcmds \
    $@



# think to add profiles for glicid/slurm 
#  --profile .profiles/local \

#bash run_example.sh --cores 1 --dag | dot -Tsvg > dag_example.svg
#bash run_example.sh --cores 1 --rulegraph | dot -Tsvg > rulegraph_example.svg
