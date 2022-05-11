#!/bin/bash

for i in $(cut -f 1 sample_list.txt | grep -v Sample); do mkdir -p MLVCFS/${i}; done

DGL=/mnt/archgen/DAG_GL
SD=/mnt/archgen/DAG_GL/scripts
#export PATH=/home/luca_traverso/.conda/envs/snakemake_lt/bin:/usr/bin/:/opt/sge/bin/:/bin:/usr/sbin:/usr/bin/X11:/sbin:/opt/sge/bin/lx-amd64

cd ${DGL}

snakemake -s ${SD}/ATLAS_GLIMPSE.smk \
--wait-for-files \
--use-conda \
--keep-going \
-j 20 \
--verbose \
--configfile ${SD}/atlas_config.yaml \
--cluster-config ${SD}/cluster2.config \
--latency-wait 20 \
--cluster "qsub -V \
-N {cluster.name} \
-pe smp {cluster.ncpus} \
-l virtual_free={cluster.mem},h_vmem={cluster.mem} \
-q archgen.q \
-cwd " \
$1 $2 $3 $4
#,h=bionode55.eva.mpg.de
 
bash scripts/clean.sh
