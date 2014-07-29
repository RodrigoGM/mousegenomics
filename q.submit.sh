#!/bin/bash
for G in $(ls *.dd.bam)
do
    ./q.indel.sh $G &
done
