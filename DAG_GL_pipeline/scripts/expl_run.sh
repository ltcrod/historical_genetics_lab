#!/bin/bash

SAMPLE=GBZ005.bam
RH=${SAMPLE}.runexpl.bam
cp ${SAMPLE} ${RH}

for i in $(samtools view -H $SAMPLE | grep RG | cut -f2 | sed 's/ID\://g')
do
echo $i

echo changing the reads - $i
samtools view $RH \
  | sed "s/${i}/${i}.mark/g" \
  > reread_$SAMPLE
samtools view reread_$SAMPLE | head -2
echo

echo fixing the header - $i
samtools view -H $RH \
  | sed "s/${i}/${i}.mark/g" \
  > ${SAMPLE}_header.sam
grep RG ${SAMPLE}_header.sam

echo samtools reheader - $i
samtools reheader \
  header.sam \
  reread_$SAMPLE \
  > ${RH}
samtools view $RH | head -2
echo
echo -----------

done
