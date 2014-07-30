#!/bin/bash
#SBATCH --time=48:00:00

## kill if anything fails
set -e
## ensure consistency of paths, modules, etc...
source ~/.bashrc
## load environment variables containing Grcm38 file paths
source $(find $GLOBALSCRATCH/tgv -name 'envar.sh' -type f)
## path/to/GATK
GATK=~/programs/GATK-latest/

# Input file at the comand line
input=$1

echo "============== In/Del Realingment + Base Recalibration  =================="
echo `date`
# -P = scatter parameter

java -Xmx500M -jar ${GATK}/Queue.jar \
    -S q.indel.scala \
    -P 35 \
    -R $MM10 \
    -I $input \
    -known $INDELVCF \
    -knownSites $DBSNP \
    -N 2 -C 2 \
    -jobRunner Drmaa \
    -jobNative "--time=48:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4000" \
#    -run -startFromScratch

echo `date`
echo "============== Finished =================="
