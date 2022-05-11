#!/bin/bash

for i in $(cut -f 1 sample_list.txt | grep -v Sample); do mkdir -p MLVCFS/${i}; done



ncores=$1
DGL=/mnt/archgen/DAG_GL
SD=/mnt/archgen/DAG_GL/scripts
PATH=/home/luca_traverso/.conda/envs/snakemake_lt/bin:/usr/bin/


mkdir -p \
renamedreads \
infos \
smdir \
outdir

snakemake \
-s ${SD}/ATLAS_GLIMPSE.smk \
--cores $ncores \
--configfile ${SD}/atlas_config.yaml  \
$2 $3 \
$4 $5 
