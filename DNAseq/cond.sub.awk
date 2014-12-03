gzcat mgp.v3.snps.rsIDdbSNPv137.annot.vcf.gz | awk 'BEGIN{OFS="\t"}{if($3 == "."){i++;$3="uk"i};print }' | bgzip -c > mgp.v3.snps.rsIDdbSNPv137.UKannot.vcf.gz
