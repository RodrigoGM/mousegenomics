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

java -Xmx1G -jar ${GATK}/Queue.jar \
    -S q.hc.scala \
    -P 64 \
    -R $MM10 \
    -I $input \
    -gt_mode Discovery --pcr_indel_model CONSERVATIVE \
    -maxAltAlleles 10 \
    -ERC GVCF \
    -pairHMM VECTOR_LOGLESS_CACHING \
    --dbsnp $DBSNP \
    -C 8

# -run
#   --variant_index_type LINEAR --variant_index_parameter 128000 \
