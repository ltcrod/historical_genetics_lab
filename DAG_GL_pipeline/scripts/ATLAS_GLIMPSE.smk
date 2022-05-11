shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")

import os
import shutil
import pandas as pd
import sys

###########################################################

configfile: "/mnt/archgen/DAG_GL/scripts/atlas_config.yaml"
runname=config['runname']
renamedreads="renamedreads/"
infofile=runname+".results.txt"
infos="infos/"
smdir="s_m/"
outdir="MLVCFS/"

df = pd.read_csv(config['sample_file'], sep='\t', index_col='Sample', comment='#')
SAMPLES= df.index
BAM_PATH= df["Path"]

###########################################################
listsamp=[]
for bam in SAMPLES:
    listsamp.append(bam)

###########################################################

FXDBAMS=[]
for bam in SAMPLES:
    FXDBAMS.append(renamedreads+bam+"_fRG.bam")

RGINFO=[]
for bam in SAMPLES:
    RGINFO.append(infos+bam+".RGINFO.txt")

SM=[]
for bam in SAMPLES:
    SM.append(smdir+bam+"_mergedReads.bam")

PMDEXP=[]
for bam in SAMPLES:
    PMDEXP.append(infos+bam+"_PMD_exponential.txt")

outvcf=["yadayada.vcf"]

outrecal=["upd.txt"]


# COMPROUT=[]
# for bam in SAMPLES:
#     COMPROUT.append(outdir+bam+"_MaximumLikelihood.vcf.gz")

IMPUTED=[]
for bam in SAMPLES:
    IMPUTED.append(outdir+bam+"/"+bam+"_imputed.vcf.gz")

log=[]
for bam in SAMPLES:
    log.append(outdir+bam+"/"+bam+".imput.log")


TARGETS=[]
# TARGETS.append(RGINFO)
# TARGETS.append(FXDBAMS)
#TARGETS.append(COMPROUT)
TARGETS.append(log)
TARGETS.append(IMPUTED)

###########################################################


localrules: all
rule all:
    input:
        TARGETS
###########################################################
rule assign2run:
    input:
        bam="samples/"+"{bam}.bam"
    output:
        bamass="samples/"+"{bam}.assigned.bam"
    params:
        runname=runname
    shell:
        """
X=$(samtools view -H {input.bam} | grep @RG | wc -l)
if (( $X >= 2 ))
then
    let X-=1
    samtools split {input.bam}
    for i in $(ls {wildcards.bam}*_{{0..20}}.bam 2>/dev/null| sed 's/\.bam//g')
    do
        ID=$(samtools view -H ${{i}}.bam | grep @RG | cut -f2 | sed 's/.*\://g')
        LB=$(samtools view -H ${{i}}.bam | grep @RG | cut -f3 | sed 's/.*\://g')
        PL=$(samtools view -H ${{i}}.bam | grep @RG | cut -f4 | sed 's/.*\://g')
        SM=$(samtools view -H ${{i}}.bam | grep @RG | cut -f5 | sed 's/.*\://g')
        PU=$(samtools view -H ${{i}}.bam | grep @RG | cut -f6 | sed 's/.*\://g')
        PI=$(samtools view -H ${{i}}.bam | grep @RG | cut -f7 | sed 's/.*\://g')
        PM=$(samtools view -H ${{i}}.bam | grep @RG | cut -f8 | sed 's/.*\://g')
        DS=$(samtools view -H ${{i}}.bam | grep @RG | cut -f9 | sed 's/.*\://g')

        /projects1/tools_new/picardtools/1.140/bin/picard AddOrReplaceReadGroups \
        I=${{i}}.bam \
        O=${{i}}.2merge.bam \
        RGID=${{ID}}.{params.runname} \
        RGLB=$LB \
        RGPL=$PL \
        RGSM=$SM \
        RGPU=$PU \
        RGPI=$PI \
        RGPM=$PM \
        RGDS=$DS

        samtools index ${{i}}.2merge.bam
        echo ${{i}}.2merge.bam >> .list{wildcards.bam}merge

        rm ${{i}}.bam ${{i}}.bai
    done

    samtools merge {output.bamass} -b .list{wildcards.bam}merge
    rm .list{wildcards.bam}merge

    rm {wildcards.bam}*_*.2merge.bam {wildcards.bam}*_*.2merge.bai

elif (( $X == 1 ))
then
    ID=$(samtools view -H {input.bam} | grep @RG | cut -f2 | sed 's/.*\://g')
    LB=$(samtools view -H {input.bam} | grep @RG | cut -f3 | sed 's/.*\://g')
    PL=$(samtools view -H {input.bam} | grep @RG | cut -f4 | sed 's/.*\://g')
    SM=$(samtools view -H {input.bam} | grep @RG | cut -f5 | sed 's/.*\://g')
    PU=$(samtools view -H {input.bam} | grep @RG | cut -f6 | sed 's/.*\://g')
    PI=$(samtools view -H {input.bam} | grep @RG | cut -f7 | sed 's/.*\://g')
    PM=$(samtools view -H {input.bam} | grep @RG | cut -f8 | sed 's/.*\://g')
    DS=$(samtools view -H {input.bam} | grep @RG | cut -f9 | sed 's/.*\://g')

    /projects1/tools_new/picardtools/1.140/bin/picard AddOrReplaceReadGroups \
    I={input.bam} \
    O={output.bamass} \
    RGID=${{ID}}.{params.runname} \
    RGLB=$LB \
    RGPL=$PL \
    RGSM=$SM \
    RGPU=$PU \
    RGPI=$PI \
    RGPM=$PM \
    RGDS=$DS

fi

TLC=$(echo {wildcards.bam} | sed 's/[0-9]*//g')

for i in TF SG
do

mkdir -p /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null

mkdir -p /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null

done
        """

rule make_RG_info:
    input:
        bam="samples/"+"{bam}.assigned.bam"
    output:
        readgroupinfo=infos+"{bam}.RGINFO.txt",
        RGs=infos+"{bam}.RG.txt"
    params:
        runname=runname
    shell:
        """
PI=200



samtools view  -H \
{input.bam} \
| grep @RG \
| cut -f 2 \
| sed 's/ID\://g' > {output.RGs}

for ID in $(cat {output.RGs} | sed 's/\.{params.runname}//g')
do
    if [[ $ID == *UM ]]
    then
        echo -e ${{ID}}.{params.runname}'\t'paired'\t' >> {output.readgroupinfo};
    else
        echo -e ${{ID}}.{params.runname}'\t'single'\t'${{PI}} >> {output.readgroupinfo}
    fi
done
        """



rule splitmerge:
    input:
        bam="samples/"+"{bam}.assigned.bam",
        readgroupinfo=infos+"{bam}.RGINFO.txt"
    output:
        bamSeM=smdir+"{bam}_mergedReads.bam"
    shell:
        """
samtools index {input.bam}

{config[atlas]} task=splitMerge \
bam={input.bam} \
out={output.bamSeM} \
readGroupSettings={input.readgroupinfo} \
allowForLarger

mv {output.bamSeM}_mergedReads.bam {output.bamSeM}
mv {output.bamSeM}_mergedReads.bam.bai {output.bamSeM}.bai

        """

rule makepoolsPMD:
    input:
        bamSeM=smdir+"{bam}_mergedReads.bam"
    output:
        poolRG=infos+"{bam}_PMD_pools.txt"
    shell:
        """
samtools view \
-H \
{input.bamSeM} \
| grep @RG \
| cut -f 2 \
| tr '\n' '\t' \
| sed 's/ID\://g' \
> {output.poolRG}

ln -S {output.poolRG} /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null

ln -S {output.poolRG} /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null
        """

rule PMDestimation:
    input:
        bamSeM=smdir+"{bam}_mergedReads.bam",
        poolRG=infos+"{bam}_PMD_pools.txt"
    params:
        prefix=infos+'{bam}'
    output:
        exponentialPMD=infos+"{bam}_PMD_input_Empiric.txt"
    message:
        "PMD-estimation on Bamfile with split readgroups on {wildcards.bam}"
    shell:
        """
{config[atlas]} \
    task=PMD \
    bam={input.bamSeM} \
    poolReadGroups={input.poolRG} \
    out={params.prefix} \
    fasta={config[ref]} \
    length=70 \
    minQual=30

echo {output.exponentialPMD}
        """




#creating the lists
formerge=[]
for bam in SAMPLES:
    formerge.append(smdir+bam+"_mergedReads.bam")

mergedpmd=[]
for bam in SAMPLES:
    mergedpmd.append(infos+bam+"_PMD_input_Empiric.txt")



#merging ALL the vcfs
rule samtools_merge:
    input:
        [formerge]
    params:
        "-R 1"
    output:
        mergbam="merged_"+runname+".bam"
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@
    wrapper:
        "v1.0.0/bio/samtools/merge"




rule update_recal:
    input:
        pmds=[mergedpmd],
        mergbam="merged_"+runname+".bam"
    output:
        updrecal=infos+"enigma_subset.txt"
    params:
        runname=runname,
        bedf=config['recalInput']
    message:
        """
        """
    shell:
        """
        echo 1
tag=$(date | cut -f1-4 -d' ' | sed 's/\ /_/g')

if [ -f RGs_{params.runname}_pool.txt ]
then
    rm RGs_{params.runname}_pool.txt
fi

cat {input.pmds} > PMD_{params.runname}.txt

samtools index {input.mergbam}

awk '{{print $1}}' PMD_{params.runname}.txt  \
| sort | uniq > all_rg_lines.tmp
echo all_rgs:
cat all_rg_lines.tmp
echo ---

rm blacklist.tmp || true
sleep 1

grep \
-w \
-f all_rg_lines.tmp \
recalibration_enigma.txt > blacklist.tmp || true


echo blacklist:
cat blacklist.tmp
echo ---

grep -v -f blacklist.tmp all_rg_lines.tmp > rg_to_use.tmp
echo rgs to use:
cat rg_to_use.tmp
echo ---

cat rg_to_use.tmp \
| tr '\\n' ' ' \
> RGs_{params.runname}_pool.txt
echo onlynewrunRGs:
#cat newrun_RGs.tmp
echo ---


echo ----------------------
echo launching atlas recal
echo ----------------------


{config[atlas]} \
    task=recal \
    bam={input.mergbam} \
    pmdFile=PMD_{params.runname}.txt \
    out={params.runname} \
    regions={params.bedf} \
    chr=1 \
    poolReadGroups=RGs_{params.runname}_pool.txt

if [ -f {params.runname}_recalibrationEM.txt ]
then
    if [ -s {params.runname}_recalibrationEM.txt ]
    then
        echo {params.runname}_recalibrationEM.txt correctly created > mail
         /usr/sbin/sendmail luca_traverso@eva.mpg.de < mail
    else
        echo {params.runname}_recalibrationEM.txt empty - exiting > mail
         /usr/sbin/sendmail luca_traverso@eva.mpg.de < mail
         exit 1
    fi
fi

cat recalibration_enigma.txt \
{params.runname}_recalibrationEM.txt \
> updated_enigma.txt



mv recalibration_enigma.txt recalibration_enigma_changed_${{tag}}.txt


cp updated_enigma.txt recalibration_enigma.txt
grep -w -f all_rg_lines.tmp recalibration_enigma.txt > {output.updrecal}

echo {params.runname} >> listruns_analysed.txt
sort listruns_analysed.txt | uniq > .rm
sleep 1s
mv .rm listruns_analysed.txt
        """


rule MLE:
    input:
        empiricPMD=infos+"{bam}_PMD_input_Empiric.txt",
        updrecal=infos+"enigma_subset.txt",
        bamSeM=smdir+"{bam}_mergedReads.bam"
    params:
        outpref=outdir+"{bam}/{bam}",
        ref=config['ref'],
        refal=config['alleles']
    output:
        mlvcf=outdir+"{bam}/{bam}_MaximumLikelihood.vcf"
    shell:
        """
{config[atlas]} \
    task=call \
    method=MLE \
    bam={input.bamSeM} \
    fasta={params.ref} \
    infoFields=DP \
    formatFields=GT,AD,AB,AI,DP,GQ,PL \
    alleles={params.refal} \
    recal={input.updrecal} \
    pmdFile={input.empiricPMD} \
    minQual=30 \
    out={params.outpref}


bgzip -d {output.mlvcf}.gz
        """

rule compress_BGZIP:
    input:
        mlvcf=outdir+"{bam}/{bam}_MaximumLikelihood.vcf"
    output:
        gzmlvcf=outdir+"{bam}/{bam}_MaximumLikelihood.vcf.gz"
    shell:
        """
bgzip -c {input.mlvcf} > {output.gzmlvcf}
tabix -p vcf {output.gzmlvcf}

ln -S {output.gzmlvcf} /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null

ln -S {output.gzmlvcf} /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null
        """

rule GLIMPSE:
    input:
        gzmlvcf=outdir+"{bam}/{bam}_MaximumLikelihood.vcf.gz"
    output:
        imputed=outdir+"{bam}/{bam}_imputed.vcf.gz",
        log=outdir+"{bam}/{bam}.imput.log"
    shell:
        """
bash scripts/GLIMPSE_bash.sh {input.gzmlvcf} \
    {config[ref]} \
    {output.imputed} \
    {wildcards.bam} 2>/dev/null

if [[ -f "{output.imputed}" ]]
then
    echo {wildcards.bam} not empty > {output.log}
else
    touch {output.imputed}
    echo {wildcards.bam} empty > {output.log}
fi

ln -S -f {output.imputed} /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null

ln -S -f {output.imputed} /mnt/archgen/Autorun_eager/eager_outputs\
/${{i}}/${{TLC}}/{wildcards.bam}/GTL_output/ 2>/dev/null
         """
