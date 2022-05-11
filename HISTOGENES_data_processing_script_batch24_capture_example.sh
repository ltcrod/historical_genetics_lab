##########################################################
## 1. Process 1240K capture data of HistoGenes          ##
##########################################################

cd /mnt/archgen/HistoGenes/eager/

rootpt="/mnt/archgen/HistoGenes/" # root path
fqpt="/mnt/archgen/data_releases/2021/211222_K00233_0246_BHMTYCBBXY_SRdi_Jena0049_JR" ## path to sequencing run in Jena
fqptLPZ="/mnt/archgen/data_releases/2021/211222_K00233_0246_BHMTYCBBXY_SRdi_Jena0049_JR" ## path to sequencing run in Leipzig
nbtch=24 ## batch
libtyp="ssLib" ## dsLib or ssLib
seqtyp="capture" ## shotgun, capture or capture2
nmask=2 # numbers of bp to mask, 2 dsLib, 10 ssLib
listfile=known_anno_b24.txt


idf="HISTOGENES_data_batch"${nbtch}"_list.txt"
inp="HISTOGENES_data_batch"${nbtch}"_inputEAGER2"
outdir="HISTOGENES_Batch"${nbtch}"_"${seqtyp}
col=3

# ## DEPRECATED: Retrieve sample IDs and add it to the list file changed to work without a idf file existing prior
# DEPRECATED: ls ${fqpt} | awk '$1 ~ /FVD/ && $1 ~ /TF/ || $1 ~ /LEO/ && $1 ~ /TF/ || $1 ~ /KDD/ && $1 ~ /TF/ || $1 ~ /SJR/ && $1 ~ /TF/ || $1 ~ /TGH/ && $1 ~ /TF/' | sed s/".TF1.1"/""/g | sed s/".TF2.1"/""/g | cat - | cut -f 1 | tail -n +1 | sort -k1,1 | uniq > bbb1


#########################
# Prepare Pandora input #
#########################

#make sure that known_anno.txt is correctly formatted and the name is correct
#e.g. try: grep Histogenes known_anno.txt
#if not, copying from a text editor and pasing in nano/vim should be enough

cp ${fqpt}/known_anno.txt ${listfile}

#make a list of the three-letter sites that are present in the sequencing list
grep -e HistoGenes -e FORMOR -e Slavs ${listfile} | grep -e TF -e SG \
| cut -f 2 \
| cut -c 1-3 \
| sort \
| uniq > list_arcsites.txt

#remove the old tmps in the directory suppressing the error message
rm cmd.tmp.sh 2>/dev/null

#create the command that creates "bbb1"
#this replaces the manual step of substitution
echo \#\!/bin/bash >> cmd.tmp.sh
echo "fqpt=\$1" >> cmd.tmp.sh
echo >> cmd.tmp.sh
echo "ls \${fqpt} \\" >> cmd.tmp.sh

for i in $(cat list_arcsites.txt)
do
echo "| \$1 ~ /"$i"/ && \$1 ~ /TF/  |" \
>> cmd.tmp.sh
done
first=$(head -1 list_arcsites.txt)
last=$(tail -1 list_arcsites.txt)
sed -i "s/\$1\ \~\ \/${first}/awk\ \'\$1\ \~\ \/${first}/g" cmd.tmp.sh
sed -i "s/${last}\/\ \&\&\ \$1\ \~\ \/TF\/\ /${last}\/\ \&\&\ \$1\ \~\ \/TF\/\ \'/g" cmd.tmp.sh
# sed -i cmd.tmp.sh 's/\'\ \|\ awk\ \'/\ \|\ awk\ /g'

echo \
" sed s/\".TF1.1\"/\"\"/g | sed s/\".TF2.1\"/\"\"/g | cat - | cut -f 1 | tail -n +1 | sort -k1,1 | uniq > bbb1" \
>> cmd.tmp.sh
sed -i ':a;N;$!ba;s/|\n/|/g' cmd.tmp.sh

#launch the freshly made command, bbb1 will be created and no one will be hurt
bash cmd.tmp.sh ${fqpt}
sleep 1s
#remove the temporary script
rm cmd.tmp.sh

##############################################################################
##############################################################################

echo -e 'ID\tShotgun\tCapture' > bbb2
while read iid; do
    cid=($(ls ${fqpt} | awk -v iid="$iid" 'BEGIN {cid="NA"} {if ($1 ~ iid) cid=$1} END {print cid}'))
    echo -e ${iid}"\t"${iid}".SG1.1\t"${cid} >> bbb2
    echo -e ${cid} >> ${rootpt}eager2/${inp}.txt
done < bbb1

mv bbb2 ${idf}; rm bbb*

cd /mnt/archgen/HistoGenes/eager2/

bash /mnt/archgen/tools/pandora2eager/0.2.2/pandora2eager.sh ${inp}.txt  \
> ${inp}.tsv

#in case pandora2eager cannot assess udg tratement, use:
#sed -i 's/none/half/g' ${inp}.tsv
#sed -i 's/Unknown/half/g' ${inp}.tsv
#sed -i "s|${fqptLPZ}|${fqpt}|g" ${inp}.tsv
#for as long as we keep using only half-UDG, of course

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
screen -S eagerscreen_batch${nbtch}
cd /mnt/archgen/HistoGenes/eager2/


nextflow run nf-core/eager -r 2.3.2 \
-profile eva,archgen,big_data,humanTK_SE \
-c /mnt/archgen/HistoGenes/1-Scripts/human_eva_bwa.config \
-with-tower \
--input ${rootpt}'eager2/'${inp}'.tsv' \
--outdir ${rootpt}'eager2/'${outdir}'/' \
-w ${rootpt}'eager2/'${outdir}'/work/' \
--email luca_traverso@eva.mpg.de \
-name  ${outdir}

######################################
## once eager run its done proceed  ##
######################################


## Generate per sample directory and symlink the results from eager2 into per-sample folder structure of eager
cd /mnt/archgen/HistoGenes/eager/

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

#old command
#sbatch -p short -c 2 --mem 8000 TEMPLATE_data_processing_script_mapDamage_210104.sh ${nbtch} ${seqtyp} ${libtyp} ${rootpt}

qsub -cwd TEMPLATE_data_processing_script_mapDamage_210104.sh ${nbtch} ${seqtyp} ${libtyp} ${rootpt}
rm output_\$\{nbtch\}

###############################################################################
## 3 Run endogenous content on 1240K off target reads                        ##
###############################################################################

#old command
# sbatch -p short -c 2 --mem 8000 TEMPLATE_data_processing_script_endo_off_210521.sh ${nbtch} ${seqtyp} ${rootpt}

cd ${rootpt}analysis/offtar_endo/

rm output_\$\{nbtch\}
qsub -cwd TEMPLATE_data_processing_script_endo_off_210521.sh ${nbtch} ${seqtyp} ${rootpt}
rm output_\$\{nbtch\}


###############################
## 4. Mask 2-bp of each end  ##
###############################

cd ${rootpt}eager/

#sbatch -p short --mem 16000 TEMPLATE_Mask_210104.sh ${nbtch} ${seqtyp} ${nmask} ${rootpt}

qsub -cwd TEMPLATE_Mask_210104.sh ${nbtch} ${seqtyp} ${nmask} ${rootpt}
rm output

##########################################################################
## 5. Calculate coverage on 1240K sites & mtDNA using 2-bp masked data  ##
##########################################################################

cd ${rootpt}analysis/coverage/

rm output_\$\{nbtch\}
qsub -cwd TEMPLATE_run_coverage_210104.sh ${nbtch} ${seqtyp} ${nmask} ${rootpt}
rm output_\$\{nbtch\}


#######################################################
## 6. X-based contamination estimation using ANGSD   ##
#######################################################

cd ${rootpt}analysis/Xcont/

rm output_\$\{nbtch\}
qsub -cwd TEMPLATE_data_processing_script_Xcont_210104.sh ${nbtch} ${seqtyp} ${nmask} ${rootpt}
rm output_\$\{nbtch\}


###############################################
## 7. Run Schmutzi using 1240K capture data  ##
###############################################

cd ${rootpt}analysis/schmutzi/

smtz_lib="single" # double or single for schmutzi

qsub -cwd TEMPLATE_data_processing_script_SE_schmutzi_210104.sh ${nbtch} ${seqtyp} ${smtz_lib} ${rootpt}
rm output_\$\{nbtch\}
##############################################
## 8. Genotype calling using pileupCaller   ##
##############################################

cd ${rootpt}analysis/genotypes/

libtyp="ssLib" ## dsLib or ssLib, for the masked_all + unmasked_Tvs or --singleStrandMode option

qsub -cwd -V TEMPLATE_data_processing_GenCall_210104.sh ${nbtch} ${seqtyp} ${libtyp} ${rootpt}
rm output_\$\{nbtch\}_gt

######################################################
## 9. Determine Y haplogroup using Yhaplo program   ##
##    This time, use the majority calling mode      ##
######################################################

cd ${rootpt}analysis/Yhap/

libtyp="ssLib" ## dsLib or ssLib, for the masked_all + unmasked_Tvs or --singleStrandMode option

qsub -cwd -V TEMPLATE_data_processing_YGenCall_210106.sh ${nbtch} ${seqtyp} ${libtyp} ${rootpt}

cat ./*/haplogroups.HISTOGENES.batch${nbtch}.yHaplo.*.txt | sort -k1,1 | awk '{OFS="\t"} {if ($2 == $3) print $1,$4"("$2")"; else print $1,$4"("$2","$3")"}'
