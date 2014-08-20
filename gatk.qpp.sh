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
    -S gatk.qpp.scala \
    -P 5 \
    -R $MM10 \
    -I $input \
    -known $INDELVCF \
    -N 6 -C 6 \
    -jobRunner Drmaa \
    -jobNative "--time=12:00:00 --nodes=1 --ntasks-per-node=6 --mem-per-cpu=4000" \
    -run -startFromScratch

echo `date`
echo "============== Finished =================="


