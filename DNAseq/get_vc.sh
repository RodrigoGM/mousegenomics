#!/bin/bash
set -e

DE_chr15=(
17:34719190-35222156
8:83358175-83861902
17:34706436-35209238
4:150753652-151294665
18:3016048-3587748
1:71335520-71903200
9:105983455-106498126
1:54802770-55322702
11:98133811-98634953
15:81552421-82076861
11:106050271-106551865
10:82980848-83587882
9:107046026-107552784
9:65044260-65580658
11:59949089-60472581
7:45455088-45960155
7:97162938-97667730
3:100201628-100739324
12:104088486-104596144
18:33187019-33714029
4:10906431-11424379
3:95628503-96132251
13:51451246-52043747
7:26585305-27093548
)
echo "Created array of B6.A-15 DE genes"

DE_chr17=(
15:82899644-83399794
9:119879351-120379553
19:39139556-39713075
19:39592660-40136769
6:47531624-48031725
7:19052721-19560045
4:132020080-132520213
)
echo "Created array of B6.A-17 DE genes"

DE_chr19=(
X:169435199-170255736
X:169435199-170255736
X:169435199-170255736
X:169435199-170255736
11:94686224-95203042
8:122019569-122522650
4:137218769-137820630
17:34420535-34969815
9:74611921-75141781
12:104088486-104596144
7:4675844-5194821
5:134453781-134997241
9:56615104-57149870
7:27055136-27587692
)
echo "Created array of B6.A-19 DE genes"

DE_chrX=(
6:84321414-84843908
1:54802770-55322702
19:37018743-37587852
16:10537936-11038655
7:46611203-47128142
1:72765074-73265909
12:104088486-104596144
1:73660231-74374449
)
echo "Created array of B6.A-X DE genes"

echo "selecting DE genes on Platypus vcf.gz"
tabix -H StrainVariants.plat.priv.vcf.gz > tmp.vcf
echo "##VariantSelectionCmd: \"echo ${DE_chr15[@]} ${DE_chr17[@]} ${DE_chr19[@]} ${DE_chrX[@]} |  xargs -n 1 -J % tabix StrainVariants.plat.priv.vcf.gz %\"" >> tmp.vcf
echo ${DE_chr15[@]} ${DE_chr17[@]} ${DE_chr19[@]} ${DE_chrX[@]} |  xargs -n 1 -J % tabix StrainVariants.plat.priv.vcf.gz % >> tmp.vcf
cat tmp.vcf | vcf-sort | java -Xmx6g -jar /opt/src/snpeff/snpEff.jar GRCm38.76 | bgzip -c > VC_cssDE.plat.priv.eff.vcf.gz && tabix VC_cssDE.plat.priv.eff.vcf.gz
rm tmp.vcf

echo "selecting DE genes on freebayes q20 vcf.gz"
tabix -H StrainVariants.freeb.q20.priv.vcf.gz > tmp.vcf
echo "##VariantSelectionCmd: \"echo ${DE_chr15[@]} ${DE_chr17[@]} ${DE_chr19[@]} ${DE_chrX[@]} |  xargs -n 1 -J % tabix StrainVariants.freeb.q20.priv.vcf.gz %\"" >> tmp.vcf
echo ${DE_chr15[@]} ${DE_chr17[@]} ${DE_chr19[@]} ${DE_chrX[@]} |  xargs -n 1 -J % tabix StrainVariants.freeb.q20.priv.vcf.gz % >> tmp.vcf
cat tmp.vcf | vcf-sort | java -Xmx6g -jar /opt/src/snpeff/snpEff.jar GRCm38.76 | bgzip -c > VC_cssDE.freeb.q20.priv.eff.vcf.gz && tabix VC_cssDE.freeb.q20.priv.eff.vcf.gz
rm tmp.vcf

#echo "selecting DE genes on GATK vqsr.99 vcf.gz"
tabix -H StrainVariants.gatk.vqsr.99.priv.vcf.gz > tmp.vcf
echo "##VariantSelectionCmd: \"echo ${DE_chr15[@]} ${DE_chr17[@]} ${DE_chr19[@]} ${DE_chrX[@]} |  xargs -n 1 -J % tabix StrainVariants.gatk.vqsr.99.priv.vcf.gz %\"" >> tmp.vcf
echo ${DE_chr15[@]} ${DE_chr17[@]} ${DE_chr19[@]} ${DE_chrX[@]} |  xargs -n 1 -J % tabix StrainVariants.gatk.vqsr.99.priv.vcf.gz % >> tmp.vcf
cat tmp.vcf | vcf-sort | java -Xmx6g -jar /opt/src/snpeff/snpEff.jar GRCm38.76 | bgzip -c > VC_cssDE.gatk.vqsr.99.priv.eff.vcf.gz && tabix VC_cssDE.gatk.vqsr.99.priv.eff.vcf.gz
rm tmp.vcf


#
#cat tmp.vcf | vcf-subset -c 4701.variant,4704.variant,4706.variant,4710.variant,4712.variant,6200.variant,7832.variant | vcf-sort | awk 'BEGIN{OFS="\t"} $3 ~ /ulg/' | uniq | java -Xmx6g -jar /opt/src/snpeff/snpEff.jar GRCm38.76 | bgzip -c > VC_cssDE.gatk.vqsr.99.eff.vcf.gz
#  366  gzcat VC_cssDE.freeb.q20.eff.vcf.gz | java -jar /opt/src/snpEff/SnpSift.jar filter "( EFF[*].IMPACT = 'HIGH' ) | ( EFF[*].IMPACT = 'MODERATE' ) | ( EFF[*].IMPACT = 'LOW' )" | bgzip -c > VC_cssDE.freeb.q20.eff.lmh_impact.vcf.gz
#  368  gzcat VC_cssDE.plat.eff.vcf.gz | java -jar /opt/src/snpEff/SnpSift.jar filter "( EFF[*].IMPACT = 'HIGH' ) | ( EFF[*].IMPACT = 'MODERATE' ) | ( EFF[*].IMPACT = 'LOW' )" | bgzip -c > VC_cssDE.plat.eff.lmh_impact.vcf.gz
