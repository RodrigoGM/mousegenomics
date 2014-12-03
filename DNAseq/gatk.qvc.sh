#!/bin/bash
#SBATCH --time=48:00:00


set -e                                                       ## kill if anything fails
source ~/.bashrc                                             ## ensure consistency of paths, modules, etc...
source $(find $GLOBALSCRATCH/tgv -name 'envar.sh' -type f)   ## load environment variables containing Grcm38 file paths
GATK=~/programs/GATK-latest/                                 ## path/to/GATK


# Input file at the comand line
# for i in $(ls *dd.bam) ; do ./gatk.qpp.sh $i & done
input=$1

echo "============== GATK Variant Discovery              =================="
echo "============== Joint Genotyping of Genomic VCF     =================="
echo "============== 1st. Variant Quality Recalibration  =================="
echo "============== Variant Selection                   =================="
echo `date`
# -P = scatter parameter

java -Xmx500M -jar ${GATK}/Queue.jar \
    -S gatk.qvcf.scala \
    -P 15 \
    -R $MM10 \
    -I $input \
    -N 8 -C 8 \
    -o QStrainVariants.vcf.gz \
    -jobRunner Drmaa \
    -jobNative "--time=12:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=4000" \
    -l DEBUG \
    -run

##    -startFromScratch -run

echo `date`
echo "============== Finished =================="


