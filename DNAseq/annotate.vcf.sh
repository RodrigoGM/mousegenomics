#!/bin/bash
vcf-merge -c any -d -s -t $1 $2 | java -jar /opt/src/snpEff/SnpSift.jar annotate $3 | awk 'BEGIN{OFS="\t"}{if($3 == ".") {i++; $3="uk"i};print}' | bgzip -c > $4