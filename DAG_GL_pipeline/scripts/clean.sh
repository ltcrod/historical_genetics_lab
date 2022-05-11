#!/bin/bash

mkdir -p old_logs

for i in makepoolsPMD make_RG_info splitmerge makepoolsPMD PMDestimation \
samtools_merge update_recal MLE
do
  mv ${i}.e* old_logs/
  mv ${i}.o* old_logs/
  mv ${i}.pe* old_logs/
  mv ${i}.po* old_logs/
done
