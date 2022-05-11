#!/bin/bash

vcf=$1 #vcf.gz
href=$2 #pass from config
outname=$3
sample=$4
wd=$(pwd)
echo $wd
#href="/mnt/archgen/Reference_Genomes/Human/hs37d5/hs37d5.fa"



mkdir -p ./GLIMPSE/; cd ./GLIMPSE/
rm -r ${sample}/ 2>/dev/null
mkdir -p ${sample}
cd ${sample}

for i in phase ligate sample ANNOT
do
    mkdir -p GLIMPSE_${i}
done

#!!!!!!
#add paths as variables through pwd+ 
#substitute $vcf with a new variable from the input



for J in $(seq 1 22)
do

    #split into chromosomes 
    tabix -h \
        ${wd}/${vcf} \
        ${J} \
        > ${sample}.chr${J}_MaximumLikelihood.vcf
    
    bgzip -c \
        ${sample}.chr${J}_MaximumLikelihood.vcf \
        > ${sample}.chr${J}_MaximumLikelihood.vcf.gz
    
    tabix -p vcf \
        ${sample}.chr${J}_MaximumLikelihood.vcf.gz
    
    bcftools index \
        -f ${sample}.chr${J}_MaximumLikelihood.vcf.gz
        


    cd GLIMPSE_phase
    #phase by chunk
    VCF="../"${sample}".chr"${J}"_MaximumLikelihood.vcf.gz"
    REF="/mnt/archgen/users/gnecchi/Jena/References/1KGP3/chr"${J}".1kg.phase3.v5a.vcf.gz"
    MAP="/mnt/archgen/users/childebayeva/bin/shapeit4/maps/chr"${J}".b37.gmap.gz"
    CHK="/mnt/archgen/users/gnecchi/Jena/References/1KGP3/Chunks/chunks.chr"${J}".txt"

    #!!!! make list
    for chunk in $(cut -f1 ${CHK})
    do
        OUT=${sample}".chr"${J}"."${chunk}".bcf"
        IRG=$(awk -v CHK=$chunk '$1 == CHK { print $3 }' $CHK)
        ORG=$(awk -v CHK=$chunk '$1 == CHK { print $4 }' $CHK)

        /projects1/clusterhomes/gnecchiruscone/GLIMPSE/static_bins/GLIMPSE_phase_static \
            --input ${VCF} \
            --reference ${REF} \
            --map ${MAP} \
            --input-region ${IRG} \
            --output-region ${ORG} \
            --output ${OUT}

        bcftools index -f ${OUT}
    done

    cd ..
    LST="GLIMPSE_ligate/list.chr"${J}".txt"
    ls GLIMPSE_phase/${sample}".chr"${J}.*.bcf > ${LST}
    OUT="GLIMPSE_ligate/"${sample}".chr"${J}".merged.bcf"

    /projects1/clusterhomes/gnecchiruscone/GLIMPSE/static_bins/GLIMPSE_ligate_static \
        --input ${LST} \
        --output ${OUT}

    bcftools index \
        -f ${OUT}



    VCF="GLIMPSE_ligate/"${sample}".chr"${J}".merged.bcf"
    OUT="GLIMPSE_sample/"${sample}".chr"${J}".phased.bcf"
    /projects1/clusterhomes/gnecchiruscone/GLIMPSE/static_bins/GLIMPSE_sample_static \
        --input ${VCF} \
        --solve \
        --output ${OUT}

    bcftools index \
        -f ${OUT}



    #sample
    VCF="GLIMPSE_ligate/"${sample}".chr"${J}".merged.bcf"
    OUT="GLIMPSE_sample/"${sample}".chr"${J}".sample.bcf"
    
    /projects1/clusterhomes/gnecchiruscone/GLIMPSE/static_bins/GLIMPSE_sample_static \
    --input ${VCF} \
    --sample \
    --output ${OUT}

    bcftools index -f ${OUT}

    bcftools annotate \
    -a GLIMPSE_ligate/${sample}.chr${J}.merged.bcf \
    -c FORMAT/GP \
    GLIMPSE_sample/${sample}.chr${J}.phased.bcf \
    -Oz \
    -o GLIMPSE_ANNOT/${sample}.chr${J}.MERGED.bcf
    
    bcftools index GLIMPSE_ANNOT/${sample}.chr${J}.MERGED.bcf



done

bcftools concat \
GLIMPSE_ANNOT/${sample}.chr*.MERGED.bcf \
-Oz \
-o ${wd}/${outname}
