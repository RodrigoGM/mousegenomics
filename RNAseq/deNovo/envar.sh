#!/bin/bash

# path to genome assembly and annotation files
# MMU = Genome refrence release; MMUFA = Genome Reference Sequennce (FASTA)
# BWT2 = Bowtie2 index; GTF = genes.gtf for release ; REFFLAT = refFlat.txt (GenePred Format)
# LIBTYPE = RNA-seq library type (either by command line or written)

# genome fasta
GLOBALSCRATCH=${GLOBALSCRATCH:-/home/rjgularte/}
RELEASE=
MMU=$GLOBALSCRATCH/mus_musculus/Ensembl/Grcm38.p4/


MMUFA=${MMU}/Sequence/WholeGenomeFasta/genome.fa

MMU1tx=${MMU}/Sequence/TxFasta/MMU1tx.fa
MMU3tx=${MMU}/Sequence/TxFasta/MMU3tx.fa

# path to indexes
## MMU1tx, MMU3tx = multi-lib assemblies
## M195tx, M156tx, single mouse assembly
BWA=${MMU}/Sequence/BWAIndex/genome.fa
BWA_MMU1tx=${MMU}/Sequence/BWAIndex/MMU1tx.fa
BWA_MMU3tx=${MMU}/Sequence/BWAIndex/MMU3tx.fa
BWA_M195tx=${MMU}/Sequence/BWAIndex/M195tx.fa

BWT2=${MMU}/Sequence/Bowtie2Index/genome
BWT2_MMU1tx=${MMU}/Sequence/Bowtie2Index/MMU1tx
BWT2_MMU3tx=${MMU}/Sequence/Bowtie2Index/MMU3tx
BWT2_M195tx=${MMU}/Sequence/Bowtie2Index/M195tx

RSEM_MMU1tx=${MMU}/Sequence/RSEMIndex/MMU1tx.fa
RSEM_MMU3tx=${MMU}/Sequence/RSEMIndex/MMU3tx.fa

RTG=${MMU}/Sequence/RTGIndex/genomeSDF/
RTGaa=${MMU}/Sequence/RTGIndex/genesSDF/
RTG_MMU1tx=${MMU}/Sequence/RTGIndex/MMU1tx
RTG_MMU3tx=${MMU}/Sequence/RTGIndex/MMU3tx
RTG_M195tx=${MMU}/Sequence/RTGIndex/M195tx
RTG_M157tx=${MMU}/Sequence/RTGIndex/M157tx

STAR2PG=${MMU}/Sequence/Star2pass/

# path to gene annotations
FAA=${MMU}/Annotation/Genes/Predicted_gene.faa
GTF=
REFFLAT=
GFF=

# path to variant annotations
# snp
DBSNP=${MMU}/Annotation/Variation/

# in/dels
INDELVCF=${MMU}/Annotation/Variation/

# scructural variants
SV_ALL=${MMU}/Annotation/Variation/
SV_NE=${MMU}/Annotation/Variation/


