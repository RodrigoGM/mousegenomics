#!/bin/bash

## For seting up a Mac OS X genomics core, the best tool to use is Homebrew.
##  We will install this first, and then proceed with additional packages not
##  found in Homebrew.
## Hombrew requires git, please have git installed.

ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

## compilers and VC tools
brew install gcc cmake ant git subversion cvs

## other useful tools
brew install emacs python perl

## notable packages for bioinformatics / genomics / genetics
brew install bamtools bcftools bedtools blast blat bowtie2 bwa circos cufflinks  genometools htslib picard-tools plink primer3 rmblast repeatmasker sambamba samtools vcftools

## other tools

echo "downloading STAR"
git clone https://github.com/alexdobin/STAR.git
cd ../

echo "downloading freebayes"
git clone --recursive https://github.com/ekg/freebayes.git
cd freebayes 
make
cd vcflib
make
cd ../../

echo "downloading snpEff"
wget -O snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download
unzip snpEff_latest_core.zip

