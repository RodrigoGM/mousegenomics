#!/bin/bash

## path to genome assembly and annotation files
## MMU = Genome refrence release; MM10 = Genome Reference Sequennce (FASTA)
## BWT2 = Bowtie2 index; GTF = genes.gtf for release ; REFFLAT = refFlat.txt (GenePred Format)
## LIBTYPE = RNA-seq library type (either by command line or written)
MMU=Grcm38
MM10=${GLOBALSCRATCH}/mus_musculus/Ensembl/${MMU}/Sequence/WholeGenomeFasta/genome.fa
BWA=${GLOBALSCRATCH}/mus_musculus/Ensembl/${MMU}/Sequence/BWAIndex/genome.fa
BWT2=${GLOBALSCRATCH}/mus_musculus/Ensembl/${MMU}/Sequence/Bowtie2Index/genome
STARG=${GLOBALSCRATCH}/mus_musculus/Ensembl/${MMU}/Sequence/StarGenome/
GTF=${GLOBALSCRATCH}/mus_musculus/Ensembl/${MMU}/Annotation/Genes/${MMU}.genes.gtf
REFFLAT=${GLOBALSCRATCH}/mus_musculus/Ensembl/${MMU}/Annotation/Genes/${MMU}.refFlat.txt
GFF=${GLOBALSCRATCH}/mus_musculus/Ensembl/${MMU}/Annotation/Genes/${MMU}.gff
LIBTYPE=fr-firststrand


