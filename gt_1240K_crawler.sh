#!/bin/bash

#THE PSEUDOAPLOID GT CRAWLER
#Usage bash gt_1240K_crawler.sh <list_pops.txt> <outprefix>
#list_pops.txt is a list of populations to be included, one per line

list=$1
outprefix=$2
trident="/home/luca_traverso/scripts/./trident-Linux"

#Reich's dataset pseudoaploid 1240K GTs
reich1240k="/mnt/archgen/users/luca_traverso/Keszthely/mking_ds/220531/dataset_kesz_qpadm/v50.0_1240k_public.trueGT"
date=$(date +%Y%m%d)

#HistoGenes pipeline latest pseudohaploid 1240K GTs
lastHG=$(ls /mnt/archgen/HistoGenes/analysis/genotypes/1240K*b{0..100}*.geno -t 2>/dev/null \
| tail -1 \
| sed 's/\.geno//g')

poseidondir="/mnt/archgen/poseidon/published_data/"


if [ -z "$lastHG" ]; then
    echo "No 1240K genotypes found from our pipeline, did the dir structure change? Please check"
    exit 1

    else echo "Last 1240K genotype from our pipeline is $lastHG, apparently"
fi

if [ -z "$reich1240k" ]; then
    echo "No 1240K genotypes for reich's dataset found. Please edit the script line 9"
    exit 1
fi

##################################################################################

#Crawler
#creates lists for each dataset

while read line; 
do
    echo "Processing $line"
    nHG=$(grep -cw $line ${lastHG}.ind)
    nReich=$(grep -cw $line ${reich1240k}.ind)
    nposeidon=$(grep -w $line ${poseidondir}/*/*.fam | wc -l | cut -d ' ' -f1)

    if [ $nHG -eq 0 ]; then
        echo "No 1240K genotypes found for $line in histogenes"
    fi
    if [ $nReich -eq 0 ]; then
        echo "No 1240K genotypes found for $line in Reich's dataset"
    fi
    if [ $nposeidon -eq 0 ]; then
        echo "No 1240K genotypes found for $line in Poseidon"
    fi

    if [ $nHG -eq 0 ] && [ $nReich -eq 0 ] && [ $nposeidon -eq 0 ]; then
                echo "no $line, exiting"
                echo "no $line" > ${outprefix}_${date}.log
                exit 1
    fi



    if [ $nHG -gt 0 ]; then 
        echo "Found $nHG 1240K genotypes for $line in our pipeline"
        echo $line >> .$outprefix.$date.list_pops.HG.tmp

    elif [ $nReich -gt 0 ]; then
        echo "Found $nReich 1240K genotypes for $line in Reich's dataset"
        echo $line >> .$outprefix.$date.list_pops.Reich.tmp

    elif [ $nposeidon -gt 0 ]; then
        echo "Found $nposeidon 1240K genotypes for $line in Poseidon"
        echo "Please then do this manually"
        echo "they can be added to the HistoGenes populations and a new file can be created in the analysis/genotypes directory"
        echo "----"
        echo "This population is present in the following files:"
        grep -w $line ${poseidondir}/*/*.fam | cut -d ':' -f1 | sort -u
        rm .$outprefix.$date.list_pops.*.tmp
        exit 1


    else
        echo "No 1240K genotypes found for $line"
        rm .$outprefix.$date.list_pops.*.tmp
        exit 1
    fi

done < $list

if [ -s .$outprefix.$date.list_pops.HG.tmp ]; then
    echo "extracting pops from $lastHG"
    $trident forge --forgeFile .$outprefix.$date.list_pops.HG.tmp \
    -p ${lastHG}.geno \
    --outFormat EIGENSTRAT \
    --onlyGeno \
    -o . \
    -n ${date}_HG_pops
fi

if [ -s .$outprefix.$date.list_pops.Reich.tmp ]; then
    echo "extracting pops from $reich1240k"
        $trident forge --forgeFile .$outprefix.$date.list_pops.Reich.tmp \
    -p ${reich1240k}.geno \
    --outFormat EIGENSTRAT \
    --onlyGeno \
    -o . \
    -n ${date}_REICH_pops
fi

nDS=$(ls .$outprefix.$date.list_pops.*.tmp | wc -l | cut -d ' ' -f1)
echo "Found $nDS datasets"


if [ $nDS -eq 1 ]
then
    echo "pop file = $outprefix.$date.list_pops.{}"
    for i in geno snp ind 
    do
        mv ${date}_*_pops.${i} $outprefix.$date.${i}
    done

elif [ $nDS -gt 1 ]
then
    echo "merging from the two datasets"

    $trident forge -f "" \
    -p ${date}_REICH_pops.geno \
    -p ${date}_HG_pops.geno \
    --outFormat EIGENSTRAT \
    --onlyGeno \
    -o . \
    -n ${date}_1240K_pops.merged 
else
    echo "No populations found"
fi



# TODO trident poseidon if exists


# rm ${date}_HG_pops* 2> /dev/null
# rm .$outprefix.$date.list_pops.*.tmp 2> /dev/null
# rm ${date}_REICH_pops.* 2> /dev/null
# rm ${date}_HG_pops.* 2> /dev/null



echo --------
echo "Done"
echo --------