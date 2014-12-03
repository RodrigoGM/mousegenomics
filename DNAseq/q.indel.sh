#!/bin/bash
#SBATCH --time=48:00:00

set -e
source ~/.bashrc
source $(find $GLOBALSCRATCH/tgv -name 'envar.sh' -type f)
GATK=~/programs/GATK-latest/
#list of files
input=$1

echo "============== realigning in/dels  =================="
echo `date`
# -P = scatter parameter
#    -glm BOTH \ is forced within UG.scala
#    -l DEBUG \ is forced within UG.scala
#    -o $rawVCF \ is hard-encoded within UG.scala

java -Xmx500M -jar ${GATK}/Queue.jar \
    -S q.indel.scala \
    -P 35 \
    -R $MM10 \
    -I $input \
    -known $INDELVCF \
    -N 2 \
    -jobRunner Drmaa -jobNative "--time=48:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4000" \
    -startFromScratch -run


##    --filter_bases_not_stored \
##    -known $SV_NE \

# submitting next job
#sbatch gatk.brc.slq $GATK $input

echo `date`
echo "============== Finished =================="
