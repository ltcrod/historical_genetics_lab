#!/bin/bash

tag=$(date | cut -f2-4 -d" " | sed 's/\ /_/g')
EO="/mnt/archgen/Autorun_eager/eager_outputs"
GLD="/mnt/archgen/DAG_GL"
AR="/mnt/archgen/Autorun/Results"

mkdir -p old_triggers
mv listruns_*.txt old_triggers

#1 update list of runs
l ${AR}/Human_1240k/ ${AR}/Human_Shotgun/ \
  | grep \/ \
  | sort \
  | uniq \
  | grep -v \: \
  | sed 's/\///g' \
  > listruns_${tag}.txt


for run in $(cat listruns_${tag}.txt)
do
#2 retrieve sample lists
  l /mnt/archgen/Autorun/Results/Human_1240k/${run}/ \
  /mnt/archgen/Autorun/Results/Human_Shotgun/${run}/ \
    | grep \/ \
    | sort \
    | uniq \
    | grep -v \: \
    | sed 's/\///g' \
    #| sed 's/\..*//g' \
    > .${run}_${tag}_samples.txt

  for sample in $(cat .${run}_${tag}_samples.txt)
  do
      id=$(echo $sample | sed 's/\..*//g')
      #file containing the infos post-alignment:
      ARRES=${AR}"/"${RUN}"/"${sample}"/Results.txt"
      #newest multiqc for that sample:
      MQC=$(ls -t ls ${EO}/*/*/${sample}/multiqc/multiqc_report.html \
        | head -1)
      #newest ATLAS report:
