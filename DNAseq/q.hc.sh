#!/bin/bash
#SBATCH --time=48:00:00

set -e
source ~/.bashrc
source $(find $GLOBALSCRATCH/tgv -name 'envar.sh' -type f)
GATK=~/programs/GATK-latest/

input=$1

echo "# ============ Run Information ============ #"

echo "Input File is    : " $input
echo "Reference Genome : " $MM10
echo "SNP db           : " $DBSNP

echo "# ============ END Run Information ============ #"

echo "# ============ Queue running HC ============ #"

java -Xmx500M -jar ${GATK}/Queue.jar \
    -S q.hc.scala \
    -P 20 \
    -R $MM10 \
    -I $input \
    -maxAltAlleles 10 \
    --dbsnp $DBSNP \
    -jobRunner Drmaa \
    -jobNative "--time=48:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4000" \
    -C 6 


##    -run


## -startFromScratch
##    -gt_mode Discovery --pcr_indel_model CONSERVATIVE \
##    --variant_index_type LINEAR --variant_index_parameter 128000 \
##    -ERC GVCF \
##    -pairHMM VECTOR_LOGLESS_CACHING \
