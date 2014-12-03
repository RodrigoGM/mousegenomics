#!/bin/bash
#SBATCH --time=48:00:00


set -e                                                       ## kill if anything fails
source ~/.bashrc                                             ## ensure consistency of paths, modules, etc...
source $(find $GLOBALSCRATCH/tgv -name 'envar.sh' -type f)   ## load environment variables containing Grcm38 file paths
GATK=~/programs/GATK-latest/                                 ## path/to/GATK


# Input file at the comand line
# for i in $(ls *dd.bam) ; do ./gatk.qpp.sh $i & done
input=$1

echo "============== GATK Pipeline  =================="
echo "============== In/Del Realignment  =================="
echo "============== BasePair Recalibration  =================="
echo "============== 2nd. BasePair Recalibration  =================="
echo "============== Variant Discovery w/ Haplotype Caller  =================="
echo `date`
# -P = scatter parameter

java -Xmx500M -jar ${GATK}/Queue.jar \
    -S gatk.qpp.scala \
    -P 15 \
    -R $MM10 \
    -I $input \
    -known $INDELVCF \
    -N 6 -C 6 \
    -jobRunner Drmaa \
    -jobNative "--time=12:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4000" \
    -l DEBUG

##    -run

##    -startFromScratch -run

echo `date`
echo "============== Finished =================="


