#!/bin/bash

## generate Splice-Junction databse

for SJ in $(find $GLOBALSCRATCH/tgv/sequencing/RNAseq/ -name 'SJ.out.tab' | sort )
do
    awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' $SJ >> SJ.Pass1.merged.sjdb
done

sort -k1,1V -k2n SJ.Pass1.merged.sjdb | uniq > SJ.Pass1.sjdb