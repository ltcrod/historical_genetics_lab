#!/bin/bash

tag=$(date | cut -f1-4 -d' ' | sed 's/\ /_/g')


for i in $(ls MLVCFS/ \
 | grep _MaximumLikelihood.vcf.gz \
 | grep -v .tbi \
 | sed 's/_MaximumLikelihood\.vcf\.gz//g')
do
    #when directories ready & true runs, do "for dir in TF SG"

    echo    $i   
    echo creating the output summary directory outdir/${i}
    mkdir -p /mnt/archgen/Autorun_eager/eager_outputs_old/TF/${i}/genotypes/

    ln -S -f /mnt/archgen/DAG_GL/MLVCFS/${i}/${i}_MaximumLikelihood.vcf.gz \
     /mnt/archgen/Autorun_eager/eager_outputs_old/TF/${i}/genotypes/

    ln -S -f /mnt/archgen/DAG_GL/MLVCFS/${i}/${i}_i*puted.vcf.gz \
     /mnt/archgen/Autorun_eager/eager_outputs_old/TF/${i}/genotypes/

    #ln the infos in the outdir directory
    echo logging in outdir/${i}
    mkdir -p outdir/${i}
    ln -s -f MLVCFS/${i}/${i}*vcf* outdir/${i}
    ln -s -f infos/*${i}* outdir/${i}


    echo removing GLIMPSE intermediates for $i
    echo directory last cleaned ${tag} \
     > /mnt/archgen/DAG_GL/GLIMPSE/out/glimpse_process/${i}/lastclean

    rm /mnt/archgen/DAG_GL/GLIMPSE/out/glimpse_process/${i}/*vcf*
    rm /mnt/archgen/DAG_GL/GLIMPSE/out/glimpse_process/${i}/GLIMPSE_ligate/*bcf* 
    rm /mnt/archgen/DAG_GL/GLIMPSE/out/glimpse_process/${i}/GLIMPSE_ligate/*chr*.txt 
    rm /mnt/archgen/DAG_GL/GLIMPSE/out/glimpse_process/${i}/GLIMPSE_phase/*bcf* 
    rm /mnt/archgen/DAG_GL/GLIMPSE/out/glimpse_process/${i}/GLIMPSE_sample/*bcf* 

    echo ----------------------------------------

done

echo ----
echo ----
echo making a log directory for the last run
runname=$(ls -t PMD*.txt | head -1 | sed 's/PMD_//g' | sed 's/\.txt//g')

mkdir -p outdir/$runname 
mv *${runname}* outdir/${runname}/
