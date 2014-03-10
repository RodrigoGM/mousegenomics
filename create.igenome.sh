#!/bin/bash
set -e

# copy directory structure from an existing iGenome.
# start at the distribution level, such that, the directory
# contains made contains the three standard directories 
# for the iGenome : ../Sequence ../Annotation ../GenomeStudio

echo " # ====   Running script to create a iGenome Like   ==== # "
echo " # ====       Reference for mapping NGS data        ==== # "
echo " # ====   THIS SCRIPT IS SPECIFIC FOR MOUSE GRCm38  ==== # "

## script requires a few arguments :
## REF = Name of the genome reference release
## CURRIG = name of current iGenome if you need to copy any file
## IGENOME = path to the location of where you want the iGenome
##    to be built
REF=Grcm38
CURRIG=/globalfs/rgularte/mus_musculus/Ensembl/NCBIM37
IGENOME=/globalfs/rgularte/mus_musculus/Ensembl/$REF

## this chunk of code copies directory structure from an exisitng 
## iGenome.  The iGenome uses Archives, and recursively links to 
## the archive-YYYY-MM-DD-hh-mm-ss to an archive-current, and
## the appropriate higher level directories.
## for simplicity, we will not use Archive linking.

#cd $CURRIG && find * -type d -exec bash -c 'mkdir -p ../Grcm38/$1' - {} \; && cd ../
#cd $IGENOME/Annotation
#mkdir Genes SmallRNA Variation

## Thus, we create the directory backbone
echo " Creating directory structure "
mkdir -p $IGENOME/Annotation/Genes $IGENOME/Annotation/SmallRNA $IGENOME/Annotation/Variation  $IGENOME/Sequence/Chromosomes  $IGENOME/Sequence/WholeGenomeFasta  $IGENOME/Sequence/AbundantSequences $IGENOME/Sequence/BowtieIndex  $IGENOME/Sequence/Bowtie2Index  $IGENOME/Sequence/BWAIndex  $IGENOME/Sequence/Squashed-Mus_musculus-Ensembl-${REF}

## download chromosomes from ensembl
echo " Donwloading individual chromosomes from Ensembl :"
echo " ftp.ensembl.org/pub/release-75/fasta/mus_musculus/dna/"
cd $IGENOME/Sequence/Chromosomes
rsync -avz rsync://ftp.ensembl.org/ensembl/pub/release-75/fasta/mus_musculus/dna/*dna.chromosome.*gz .

## check if file exists and remove
if [ -e ChromInfo ]
    then rm ChromInfo
fi

## create a Chromosome Information file in */Genes
for gz in $(ls *gz) 
do  
    zcat $gz | head | grep -e ">"  | cut -d : -f 4,6 | sed 's/:/\t/' >> ChromInfo
done
sort -V ChromInfo > $IGENOME/Annotation/Genes/ChromInfo.txt
rm ChromInfo

cd $IGENOME/Sequence/AbundantSequences
echo "Not sure what goes in $IGENOME/Sequence/AbundantSequences." > README
echo "However, please take a look at your favorite iGenome" >> README
echo "and download the necesary *.fa files.  You may donwload" >> README
echo "the sepecific sequences by searching for the GI numbers" >> README
echo "in NCBI or ensembl." >> README
echo "This file should contain : " >> README
echo "Sequence in FASTA format of mitochondrial genome 'MT.fa'" >> README
echo "Sequence in FASTA format from Adapters, names 'adapter_contam1.fa'">> README
echo "Sequence in FASTA format from Mus musculus Ribosomal RNA gi:38176281 ; BK000964.x">> README
echo "Sequence in FASTA format of Enterobacteria phage phiX174 gi:9626372; NC_001422.x">> README
echo "Sequence in FASTA format of Poly A (A repeat 100mer)">> README
echo "Sequence in FASTA format of Poly C (C repeat 99mer)">> README
#cp $CURRIG/Sequence/AbundantSequences/* .

echo "Downloading annotation files to Annotation/Genes from Ensembl"
echo "ftp.ensembl.org/ensembl/pub/release-75/gtf/mus_musculus/"
cd $IGENOME/Annotation/Genes
rsync -avz rsync://ftp.ensembl.org/ensembl/pub/release-75/gtf/mus_musculus/* .
echo "Extracting to ${REF}.genes.gtf"
zcat  Mus_musculus.GRCm38.75.gtf.gz | sort -V > ${REF}.genes.gtf

## Download refseq from UCSC (ensGene.txt.gz and .sql) 
## Remove leading 'chr' on the chromosome number, and sort by chr and bp
rsync -avz rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/ensGene* .
zcat ensGene.txt.gz | sed -e 's/chr//g' | sort -k3V -k4n > ensGene.txt


## need to create ref flat
if [ -e ~/bin/gtfToGenePred ]
	then gtfToGenePred ${REF}.genes.gtf ${REF}.refFlat.txt
	else echo "gtfToPred not found. Please download it using : "
		 echo "wget -nH -O $HOME/bin/ http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred"
		 exit 1;
fi

echo "Downloading Small RNA and Non-Coding RNA data from :"
echo "Ensembl : ftp.ensembl.org/pub/release-75/fasta/mus_musculus/ncrna/"
echo "miRBase : ftp://mirbase.org/pub/mirbase/CURRENT/"
cd $IGENOME/Annotation/SmallRNA
rsync -avz rsync://ftp.ensembl.org/ensembl/pub/release-75/fasta/mus_musculus/ncrna/ .
mv README README.ncrna
wget -nc ftp://mirbase.org/pub/mirbase/CURRENT/*gz 
wget -nc ftp://mirbase.org/pub/mirbase/CURRENT/README

gunzip -f *.gz
mv Mus_musculus.GRCm38.75.ncrna.fa ${REF}.ncrna.fa

echo "Downloading Annotation/Variation from Ensembl"
echo "ftp.ensembl.org/pub/release-75/variation/gvf/mus_musculus/"
cd $IGENOME/Annotation/Variation
rsync -avz rsync://ftp.ensembl.org/ensembl/pub/release-75/variation/gvf/mus_musculus/ .

echo "Create whole genome fasta file in Sequence/WholeGenomeFasta"
echo "Using file Annotation/Genes/ChromInfo.txt, as base to order"
echo "chromosomes.  Chromosomes sorted as 1:19, MT, X, Y"
cd $IGENOME/Sequence/Chromosomes

echo "Cleaning ../WholeGenomeFasta/"
if [ -e $(find $IGENOME/Sequence/WholeGenomeFasta/genome.dict) ]
    then rm $IGENOME/Sequence/WholeGenomeFasta/genome.fa
    then rm $IGENOME/Sequence/WholeGenomeFasta/genome.fa.fai
    then rm $IGENOME/Sequence/WholeGenomeFasta/genome.dict
    else " No documents found in ../WholeGenomeFasta/, continuing"
fi

echo " Creating Reference Genome for release $REF"
for gz in $(cut -f 1 $IGENOME/Annotation/Genes/ChromInfo.txt)
do
    echo "Appending chromsome : " $gz
    zcat *chromosome.${gz}.fa.gz | sed 's/dna:chromosome chromosome://'>> $IGENOME/Sequence/WholeGenomeFasta/genome.fa
done

echo "# === Generating genome.fa fasta index === #"
cd $IGENOME/Sequence/WholeGenomeFasta/
samtools faidx genome.fa

echo "# == Generating sequence ditionary === #"

if [ -e $(find ~/ -name CreateSequenceDictionary.jar) ]
    then CSD=$(find ~/ -name CreateSequenceDictionary.jar)
         if [ -e genome.dict ]; then rm genome.dict; fi
         java -jar $CSD R=genome.fa O=genome.dict
    else echo "picard-tools CreateSequenceDictionary.jar not found"
         exit 1
fi

## created executable for consistancy with the original iGenome download
## for consistency we make it executable
chmod +x $IGENOME/Sequence/WholeGenomeFasta/genome.fa
chmod +x $IGENOME/Sequence/WholeGenomeFasta/genome.fa.fai
chmod +x $IGENOME/Sequence/WholeGenomeFasta/genome.dict

echo "# === Generationg BWA indexes in v6 === #"
cd $IGENOME/Sequence/BWAIndex/
cp -rs $IGENOME/Sequence/WholeGenomeFasta/genome.fa ./genome.fa
bwa index genome.fa

echo "# == Generating Bowtie indexes === #"
if [ -e $HOME/bin/bowtie-build ]
then 
    cd $IGENOME/Sequence/BowtieIndex
    cp -rs $IGENOME/Sequence/WholeGenomeFasta/genome.fa ./genome.fa
    bowtie-build -f genome.fa ./
else echo "Bowtie v1 not found, moving on"
fi

echo "# == Generating Bowtie2 indexes === #"
cd  cd $IGENOME/Sequence/Bowtie2Index
cp -rs $IGENOME/Sequence/WholeGenomeFasta/genome.fa ./genome.fa
bowtie2-build -f genome.fa ./
