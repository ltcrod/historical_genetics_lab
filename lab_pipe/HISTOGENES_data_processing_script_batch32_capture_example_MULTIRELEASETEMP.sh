#######################################################################################################
## 1. Process 1240K capture data of HistoGenes batch 29 combine 2 sequencing runs came out same day  ##
#######################################################################################################

cd /mnt/archgen/HistoGenes/eager/

rootpt="/mnt/archgen/HistoGenes/" # root path
fqpt="/mnt/archgen/data_releases/2022/220429_K00233_0265_AHNN3FBBXY_SRdi_JR_2/" ## path to sequencing run in Jena
fqpt1="/mnt/archgen/data_releases/2022/220504_K00233_0269_AHNN25BBXY_SRdi_SN_JR/" ## path to sequencing run in Jena

nbtch=30 ## batch
libtyp="ssLib" ## dsLib or ssLib
seqtyp="capture" ## shotgun, capture or capture2
nmask=2 # numbers of bp to mask, 2 dsLib, 10 ssLib
#listfile=known_anno_b21.txt


idf="HISTOGENES_data_batch"${nbtch}"_list.txt"
inp="HISTOGENES_data_batch"${nbtch}"_inputEAGER2"
outdir="HISTOGENES_Batch"${nbtch}"_"${seqtyp}
col=3

# Retrieve sample IDs and add it to the list file changed to work without a idf file existing prior
ls ${fqpt} | awk '$1 ~ /MI6/ && $1 ~ /TF/ || $1 ~ /TKI/ && $1 ~ /TF/ || $1 ~ /RKF/ && $1 ~ /TF/  || $1 ~ /RFK/ && $1 ~ /TF/ || $1 ~ /AJD/ && $1 ~ /TF/ || $1 ~ /MJA/ && $1 ~ /TF/ || $1 ~ /BGE/ && $1 ~ /TF/ || $1 ~ /MBS/ && $1 ~ /TF/ || $1 ~ /LNO/ && $1 ~ /TF/' | sed s/".TF1.1"/""/g | sed s/".TF1.2"/""/g | sed s/".TF2.1"/""/g | sed s/".TF2.2"/""/g | sed s/".TF2.3"/""/g | cat - | cut -f 1 | tail -n +1 | sort -k1,1 | uniq > bbb1

echo -e 'ID\tShotgun\tCapture' > bbb2
while read iid; do
    cid=($(ls ${fqpt} | awk -v iid="$iid" 'BEGIN {cid="NA"} {if ($1 ~ iid) cid=$1} END {print cid}'))
    echo -e ${iid}"\t"${iid}".SG1.1\t"${cid} >> bbb2
    echo -e ${cid} >> ${rootpt}eager2/${inp}.txt
done < bbb1

ls ${fqpt1} | awk '$1 ~ /MI6/ && $1 ~ /TF/ || $1 ~ /TKI/ && $1 ~ /TF/ || $1 ~ /RKF/ && $1 ~ /TF/  || $1 ~ /RFK/ && $1 ~ /TF/ || $1 ~ /AJD/ && $1 ~ /TF/ || $1 ~ /MJA/ && $1 ~ /TF/ || $1 ~ /BGE/ && $1 ~ /TF/ || $1 ~ /MBS/ && $1 ~ /TF/ || $1 ~ /LNO/ && $1 ~ /TF/' | sed s/".TF1.1"/""/g | sed s/".TF1.2"/""/g | sed s/".TF2.1"/""/g | sed s/".TF2.2"/""/g | sed s/".TF2.3"/""/g | cat - | cut -f 1 | tail -n +1 | sort -k1,1 | uniq > bbb1
while read iid; do
    cid=($(ls ${fqpt1} | awk -v iid="$iid" 'BEGIN {cid="NA"} {if ($1 ~ iid) cid=$1} END {print cid}'))
    echo -e ${iid}"\t"${iid}".SG1.1\t"${cid} >> bbb2
    echo -e ${cid} >> ${rootpt}eager2/${inp}.txt
done < bbb1

cat bbb2 | (sed -u 1q; sort -k 1,1) > ${idf}; rm bbb*

cd /mnt/archgen/HistoGenes/eager2/

sort -k1,1 ${inp}.txt > aaa

mv aaa ${inp}.txt

bash /mnt/archgen/tools/pandora2eager/0.2.2/pandora2eager.sh ${inp}.txt > ${inp}.tsv


###########################################
## screen functions example. don't run!  ##
###########################################

#Launch screen:
screen
#Detach from Linux Screen Session
Ctrl+a+d
#Reattach to a Linux Screen
screen -r
#In case you have multiple screen sessions running on your machine, you will need to append the screen session ID after the r switch.
screen -r session_name
#Starting Named Session. Named sessions are useful when you run multiple screen sessions.
screen -S session_name
#To find the session ID list the current running screen sessions with:
screen -ls
#Scroll the screen
Shift +Fn + UP or DOWN on a Macbook will allow you to scroll.


###################
## submit eager2 ##
###################
screen

rootpt="/mnt/archgen/HistoGenes/"
fqpt="/mnt/archgen/data_releases/2022/220429_K00233_0265_AHNN3FBBXY_SRdi_JR_2/" ## path to sequencing run in Jena
fqpt1="/mnt/archgen/data_releases/2022/220504_K00233_0269_AHNN25BBXY_SRdi_SN_JR/" ## path to sequencing run in Jena
nbtch=30 ## batch number
libtyp="ssLib" ## dsLib or ssLib
seqtyp="capture" ## shotgun, capture or capture2
nmask=2 # numbers of bp to mask, 2 dsLib, 10 ssLib
idf=${rootpt}"eager/HISTOGENES_data_batch"${nbtch}"_list.txt"
inp="HISTOGENES_data_batch"${nbtch}"_inputEAGER2"
outdir="HISTOGENES_Batch"${nbtch}"_"${seqtyp}


nextflow run nf-core/eager -r 2.3.2 \
-profile eva,archgen,big_data,humanTK_SE \
-c /mnt/archgen/HistoGenes/1-Scripts/human_eva_bwa.config \
-with-tower \
-name ${outdir} \
--input ${rootpt}'eager2/'${inp}'.tsv' \
--outdir ${rootpt}'eager2/'${outdir}'/' \
-w ${rootpt}'eager2/'${outdir}'/work/' \
--email gnecchiruscone@shh.mpg.de


######################################
## once eager run its done proceed  ##
######################################


## Generate per sample directory and symlink the results from eager2 into per-sample folder structure of eager
cd /mnt/archgen/HistoGenes/eager/
col=3
while read line; do
    iid=($(echo ${line} | awk '{print $1}'))
    cid=($(echo ${line} | awk -v c=$col '{print $c}'))
    cid2=($(ls ${rootpt}eager2/${outdir}/deduplication/ | awk -v id="$iid" '$1 ~ id'))
    mkdir -p ./${iid}/${seqtyp}/${cid}/1-AdapClip/; mkdir -p ./${iid}/${seqtyp}/${cid}/5-DeDup/
    ln -s ${rootpt}eager2/${outdir}/deduplication/${cid2}/* ./${iid}/${seqtyp}/${cid}/5-DeDup/
    ln -s ${rootpt}eager2/${outdir}/adapterremoval/output/${cid}* ./${iid}/${seqtyp}/${cid}/1-AdapClip/
done < <(awk -v c=$col '$c != "NA"' ${idf} | tail -n +2)


##################################################
## 2. Run mapDamage with 100,000 reads or less  ##
##################################################

cd ${rootpt}analysis/mapDamage/

rm output_\$\{nbtch\}

qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_data_processing_script_mapDamage_210104.sh ${nbtch} ${seqtyp} ${libtyp} ${rootpt}

###############################################################################
## 3 Run endogenous content on 1240K off target reads                        ##
###############################################################################

cd ${rootpt}analysis/offtar_endo/

rm output_\$\{nbtch\}

qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_data_processing_script_endo_off_210521.sh ${nbtch} ${seqtyp} ${rootpt}

###############################
## 4. Mask 2-bp of each end  ##
###############################

cd ${rootpt}eager/

rm output

qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_Mask_210104.sh ${nbtch} ${seqtyp} ${nmask} ${rootpt}

##########################################################################
## 5. Calculate coverage on 1240K sites & mtDNA using 2-bp masked data  ##
##########################################################################

cd ${rootpt}analysis/coverage/


rm output_\$\{nbtch\}

qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_run_coverage_210104.sh ${nbtch} ${seqtyp} ${nmask} ${rootpt}


#######################################################
## 6. X-based contamination estimation using ANGSD   ##
#######################################################

cd ${rootpt}analysis/Xcont/

rm output_\$\{nbtch\} 

qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_data_processing_script_Xcont_210104.sh ${nbtch} ${seqtyp} ${nmask} ${rootpt}


###############################################
## 7. Run Schmutzi using 1240K capture data  ##
###############################################

cd ${rootpt}analysis/schmutzi/

rm output_\$\{nbtch\} 

smtz_lib="single" # double or single for schmutzi

rm output_\$\{nbtch\}_gt

qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_data_processing_script_SE_schmutzi_210104.sh ${nbtch} ${seqtyp} ${smtz_lib} ${rootpt}

##############################################
## 8. Genotype calling using pileupCaller   ##
##############################################

cd ${rootpt}analysis/genotypes/

libtyp="ssLib" ## dsLib or ssLib, for the masked_all + unmasked_Tvs or --singleStrandMode option

rm output_\$\{nbtch\}_gt


qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_data_processing_GenCall_210104.sh ${nbtch} ${seqtyp} ${libtyp} ${rootpt}

######################################################
## 9. Determine Y haplogroup using Yhaplo program   ##
##    This time, use the majority calling mode      ##
######################################################

cd ${rootpt}analysis/Yhap/

rm output_\$\{nbtch\}_gt

libtyp="ssLib" ## dsLib or ssLib, for the masked_all + unmasked_Tvs or --singleStrandMode option

qsub -S /bin/bash -cwd -pe smp 2 -l h_vmem=8000M TEMPLATE_data_processing_YGenCall_210106.sh ${nbtch} ${seqtyp} ${libtyp} ${rootpt}

cat ./*/haplogroups.HISTOGENES.batch${nbtch}.yHaplo.*.txt | sort -k1,1 | awk '{OFS="\t"} {if ($2 == $3) print $1,$4"("$2")"; else print $1,$4"("$2","$3")"}'

#################################################################################################################################
## 10. Merge data with curated 1240K of all batches!                                                                           ##
##    ACHTUNG: Use this script only if you are sure the alleles between the different dataset are coded on the same strand     ##
##    if not or not sure use mergeit (https://github.com/DReichLab/AdmixTools/blob/master/convertf/README)                     ##
#################################################################################################################################
rootpt="/mnt/archgen/HistoGenes/"
cd ${rootpt}analysis/genotypes/

pt1=($(pwd)"/")
fn1=${rootpt}"analysis/genotypes/1240K.HISTOGENES.220420"
fn2="HISTOGENES.batch29.1240K.ssMode"
fn3="HISTOGENES.batch30.1240K.ssMode"
#fn4="HISTOGENES.batch28.1240K.ssMode"

of1="1240K.HISTOGENES.220726"
tn1="temp1_"${of1}

## Merge 1240K data into one
paste ${fn2}.geno ${fn3}.geno -d '' > ${tn1}_1.geno
cat ${fn2}.ind ${fn3}.ind > ${tn1}_1.ind
cp ${fn2}.snp ${tn1}_1.snp

awk '{print $1}' ${fn1}.snp > ${tn1}_2  ## Extract SNP IDs for HumanOrigins markers
paste ${tn1}_1.snp ${tn1}_1.geno | fgrep -wf ${tn1}_2 - | awk '{print $7}' > ${tn1}_2.geno  ## Extract HumanOrigins markers
paste ${tn1}_1.snp ${tn1}_1.geno | fgrep -wf ${tn1}_2 - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' > ${tn1}_2.snp

## Mark overlapping IIDs to update data
cn1s=""; while read cn; do
    if [[ "$cn1s" == "" ]]; then cn1s+=${cn}; else cn1s+=","${cn}; fi
done < <(cut -f 1 ${tn1}_1.ind | fgrep -wnf - ${fn1}.ind | cut -d ':' -f 1)

if [[ "$cn1s" == "" ]]; then
    cp ${fn1}.geno ${tn1}_3.geno
    cp ${fn1}.ind ${tn1}_3.ind
    cp ${fn1}.snp ${tn1}_3.snp
else
    cut --complement -c ${cn1s} ${fn1}.geno > ${tn1}_3.geno  ## remove overlapping IDs
    cut -f 1 ${tn1}_1.ind | fgrep -wvf - ${fn1}.ind > ${tn1}_3.ind
    cp ${fn1}.snp ${tn1}_3.snp
fi

paste ${tn1}_3.geno ${tn1}_2.geno -d '' > ${of1}.geno
cp ${tn1}_3.snp ${of1}.snp
cat ${tn1}_3.ind ${tn1}_1.ind | awk '{OFS="\t"} {print $1,$2,$3}' > ${of1}.ind

rm ${tn1}_*

mv ${fn1}.* old_merged/

####################################################################################
## 10. Merge data with curated 1240K + HumanOrigins data last Avar masterdata     ##
####################################################################################

rootpt="/mnt/archgen/HistoGenes/" # root path

cd ${rootpt}analysis/genotypes/

pt1=($(pwd)"/")
fn1=${rootpt}"analysis/genotypes/1240KHO.poseidon2_merged.MD.Avar.220126"
fn2="HISTOGENES.batch26.1240K.ssMode"
fn3="HISTOGENES.batch27.1240K.ssMode"
fn4="HISTOGENES.batch28.1240K.ssMode"

of1="1240KHO.poseidon2_merged.MD.Avar.220420"
tn1="temp1_"${of1}

#### make par file for mergit

## Merge 1240K data into one
paste ${fn2}.geno ${fn3}.geno ${fn4}.geno -d '' > ${tn1}_1.geno
cat ${fn2}.ind ${fn3}.ind ${fn4}.ind > ${tn1}_1.ind
cp ${fn2}.snp ${tn1}_1.snp

awk '{print $1}' ${fn1}.snp > ${tn1}_2  ## Extract SNP IDs for HumanOrigins markers
paste ${tn1}_1.snp ${tn1}_1.geno | fgrep -wf ${tn1}_2 - | awk '{print $7}' > ${tn1}_2.geno  ## Extract HumanOrigins markers
paste ${tn1}_1.snp ${tn1}_1.geno | fgrep -wf ${tn1}_2 - | awk '{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' > ${tn1}_2.snp

## Mark overlapping IIDs to update data
cn1s=""; while read cn; do
    if [[ "$cn1s" == "" ]]; then cn1s+=${cn}; else cn1s+=","${cn}; fi
done < <(cut -f 1 ${tn1}_1.ind | fgrep -wnf - ${fn1}.ind | cut -d ':' -f 1)

if [[ "$cn1s" == "" ]]; then
    cp ${fn1}.geno ${tn1}_3.geno
    cp ${fn1}.ind ${tn1}_3.ind
    cp ${fn1}.original.ind ${tn1}_3.original.ind
    cp ${fn1}.2.ind ${tn1}_3.2.ind
    cp ${fn1}.3.ind ${tn1}_3.3.ind
    cp ${fn1}.snp ${tn1}_3.snp
else
    cut --complement -c ${cn1s} ${fn1}.geno > ${tn1}_3.geno  ## remove overlapping IDs
    cut -f 1 ${tn1}_1.ind | fgrep -wvf - ${fn1}.ind > ${tn1}_3.ind
    cut -f 1 ${tn1}_1.ind | fgrep -wvf - ${fn1}.original.ind > ${tn1}_3.original.ind
    cut -f 1 ${tn1}_1.ind | fgrep -wvf - ${fn1}.2.ind > ${tn1}_3.2.ind
    cut -f 1 ${tn1}_1.ind | fgrep -wvf - ${fn1}.3.ind > ${tn1}_3.3.ind
    cp ${fn1}.snp ${tn1}_3.snp
fi

paste ${tn1}_3.geno ${tn1}_2.geno -d '' > ${of1}.geno
cp ${tn1}_3.snp ${of1}.snp
cat ${tn1}_3.ind ${tn1}_1.ind | awk '{OFS="\t"} {print $1,$2,$3}' > ${of1}.ind
cat ${tn1}_3.2.ind ${tn1}_1.ind | awk '{OFS="\t"} {print $1,$2,$3}' > ${of1}.2.ind
cat ${tn1}_3.3.ind ${tn1}_1.ind | awk '{OFS="\t"} {print $1,$2,$3}' > ${of1}.3.ind
cat ${tn1}_3.original.ind ${tn1}_1.ind | awk '{OFS="\t"} {print $1,$2,$3}' > ${of1}.original.ind

rm ${tn1}_*

###########################################################
## 2. Run preliminary PCA for the new data set as a QC  ##
###########################################################
rootpt="/mnt/archgen/HistoGenes/"

cd ${rootpt}analysis/

mkdir -p ./PCA_220420/; cd ./PCA_220420/

pt1=($(pwd)"/")
fn1=${rootpt}"analysis/genotypes/1240KHO.poseidon2_merged.MD.Avar.220420"
fn2=${rootpt}"analysis/PCA_211018/PCA_211018"
of1="PCA_220420"
tn1="temp1_"${of1}

nset=($(ls ${fn2}_*.pops | wc -l))
for K in $(seq 1 $nset); do cp ${fn2}_${K}.pops ${of1}_${K}.pops; done

## Write down parameter files
nset=($(ls ${of1}_*.pops | wc -l))
for K in $(seq 1 $nset); do
    of2=${of1}"_"${K}; parf=${of2}".par"
    echo "genotypename: "${fn1}".geno" > ${parf}
    echo "snpname: "${fn1}".snp" >> ${parf}
    if [ "$K" -lt "$nset" ]; then
        echo "indivname: "${fn1}".original.ind" >> ${parf}
    else
        echo "indivname: "${fn1}".ind" >> ${parf}
    fi
    echo "evecoutname: "${pt1}${of2}".evec" >> ${parf}
    echo "evaloutname: "${pt1}${of2}".eval" >> ${parf}
    echo "poplistname: "${pt1}${of2}".pops" >> ${parf}
    echo -e "altnormstype: NO\nnumoutevec: 20\nnumoutlieriter: 0\nnumoutlierevec: 0" >> ${parf}
    echo -e "outliersigmathresh: 6.0\nnumthreads: 8\nqtmode: 0\nlsqproject: YES\nautoshrink: YES" >> ${parf}
done

## Run smartpca
for K in $(seq 1 $nset); do
    of2=${of1}"_"${K}
    qsub -cwd -N pca -b y -pe smp 8 -l h_vmem=64G smartpca -p ${of2}.par > ${of2}.log
done

## Reformat output and remove unnecessary files
/mnt/archgen/users/gnecchi/Jena/scripts/eigenstrat_output_trimmer_160222.py ${of1} ${nset}
rm ${of1}_*.log; rm pca.*



