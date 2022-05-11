#!/bin/bash

ncores=1

DGL=/mnt/archgen/DAG_GL
SD=/mnt/archgen/DAG_GL/scripts

mkdir -p \
renamedreads \
infos \
smdir \
outdir

snakemake \
-s ${SD}/ATLAS_GLIMPSE.smk \
--cores $ncores \
--dag \
--configfile ${SD}/atlas_config.yaml $1 $2\
| dot -Tpdf > dag.pdf
