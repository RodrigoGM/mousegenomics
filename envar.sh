#!/bin/bash

# path to genome assembly and annotation files
# MMU = Genome refrence release; MM10 = Genome Reference Sequennce (FASTA)
# BWT2 = Bowtie2 index; GTF = genes.gtf for release ; REFFLAT = refFlat.txt (GenePred Format)
# LIBTYPE = RNA-seq library type (either by command line or written)

# genome fasta
RELEASE=Grcm38
MMU=${GLOBALSCRATCH}/mus_musculus/Ensembl/${RELEASE}
MM10=${MMU}/Sequence/WholeGenomeFasta/genome.fa

# path to indexes
BWA=${MMU}/Sequence/BWAIndex/genome.fa
BWT2=${MMU}/Sequence/Bowtie2Index/genome
STARG=${MMU}/Sequence/StarGenome/
STAR2PG=${MMU}/Sequence/Star2pass/

# path to gene annotations
GTF=${MMU}/Annotation/Genes/${RELEASE}.genes.gtf
REFFLAT=${MMU}/Annotation/Genes/${RELEASE}.refFlat.txt
GFF=${MMU}/Annotation/Genes/${RELEASE}.gff

# path to variant annotations
# snp
DBSNP=${MMU}/Annotation/Variation/Mus_musculus.sorted.vcf
SANGERSNP=${MMU}/Annotation/Variation/mgp.v3.snps.rsIDdbSNPv137.vcf.gz
SNP_MO=${MMU}/Annotation/Variation/Mus_musculus_SNP_MO.vcf
JAX_MDA1=${MMU}/Annotation/Variation/JAX-MDA1.vcf

# in/dels
INDELVCF=${MMU}/Annotation/Variation/Mus_musculus_InDel.vcf.gz
INDELMO=${MMU}/Annotation/Variation/Mus_musculus_InDel_MO.vcf.gz
SANGERINDEL=${MMU}/Annotation/Variation/mgp.v3.indels.rsIDdbSNPv137.vcf.gz

# scructural variants
SV_ALL=${MMU}/Annotation/Variation/Mus_musculus_StrucVar.vcf.gz
SV_NE=${MMU}/Annotation/Variation/Mus_musculus_StrucVarNoEmpty.vcf

LIBTYPE=fr-firststrand
