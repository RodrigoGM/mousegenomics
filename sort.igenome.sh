#!/bin/bash

echo "# === Sorting Annotation/Genes/* === #"
cd ./Annotation/Genes/
sort -V ChromInfo.txt > ChromInfo && mv ChromInfo ChromInfo.txt
sort -V genes.gtf > genes && mv genes genes.gtf
sort -k3 -V refGene.txt > refGene && mv refGene refGene.txt
if [ -f refFlat.txt.gz ]
then zcat refFlat.txt.gz | sort -k3 -V > refFlat.txt
else sort -k3 -V refFlat.txt > refFlat && mv refFlat refFlat.txt
fi

echo "# === Creating genome.fa file ===#"

cd ../../Sequence/Chromosomes
rm ../WholeGenomeFasta/*
for chr in $(cut -f 1 ../../Annotation/Genes/ChromInfo.txt)
do
    echo "Appending chr:  " $chr
    cat ${chr}.fa >> ../WholeGenomeFasta/genome.fa
done

echo "# === Generating genome.fa fasta index === #"
cd ../WholeGenomeFasta/
samtools fadix genome.fa

echo "# == Generating sequence ditionary === #"
CSD=$(find ~/ -name CreateSequenceDictionary.jar)
if [ -f genome.dict ]; then rm genome.dict; fi
java -jar $CSD R=genome.fa O=genome.dict

## created executable for consistancy with the original iGenome download
chmod +x genome.fa
chmod +x genome.fa.fai
chmod +x genome.dict

echo "# === Generationg BWA indexes in v6 === #"
cd ../BWAIndex/version0.6.0
bwa index genome.fa

echo "# == Generating Bowtie indexes === #"
if [ -x $HOME/bin/bowtie-build ]
then 
    cd ../../BowtieIndex
    bowtie-build -f genome.fa ./
else echo "Bowtie v1 not found, moving on"
fi

echo "# == Generating Bowtie2 indexes === #"
cd ../Bowtie2Index
bowtie2-build -f genome.fa ./

