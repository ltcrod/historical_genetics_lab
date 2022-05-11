#!/usr/bin/env bash

ID=$1
#ID=keszthelyKFP_relatable.csv

rels=$(awk '$3!="un"' filtered0/merged_relatable_allLikelihoods.csv | tail -n +2 | awk '{print $3}' | sort | uniq)
rinds=$(awk '$3!="un"' filtered0/merged_relatable_allLikelihoods.csv  | tail -n +2 | awk '{print $2}' | tr '_' '\n' | sort | uniq)

echo -e "id\trelatedness\trelatedindividualID" > ${ID}_relatedness_per_individual.txt
for id in $rinds;do
    for rel in $rels;do
        relinds=$(grep $id <(awk '$3!="un"' filtered0/merged_relatable_allLikelihoods.csv | awk  '{ $5 = sprintf("%.2f", $5) }1' ) \
         | awk '{print $2"\t"$3"\t"$5}' \
         | tr '_' '\t' \
         | grep $rel \
         |  awk '{print $1"\n"$2"("$4")"}' \
         | grep -v $id \
         | paste -s -d';')

        echo -e "$id\t$rel\t$relinds" >> ${ID}_relatedness_per_individual.txt
    done
done

#open R, long to wide conversion
Rscript comp.R $ID
