#!/bin/bash
## The awk line was obtained from the RNA-STAR google groups' thread :
## https://groups.google.com/forum/#!searchin/rna-star/awk$20$27begin$20%7Bofs$3D/rna-star/yvJ6C3h7OMk/BOiN0KUJYnYJ
## Written by Alex Dobin.  It's implementation to merge multiple SJ.out.tab is not discussed in this thread. 

## generate Splice-Junction databse

for SJ in $(find $GLOBALSCRATCH/tgv/sequencing/RNAseq/ -name 'SJ.out.tab' | sort )
do
    awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' $SJ >> SJ.Pass1.merged.sjdb
done

sort -k1,1V -k2n SJ.Pass1.merged.sjdb | uniq > SJ.Pass1.sjdb