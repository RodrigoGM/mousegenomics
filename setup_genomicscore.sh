#!/bin/bash

##########################################################################
## The purpose of this script is to download, or compile software for 
##   genomics.  Manily, for pipelines for RNA-seq differential expression,
##   exon bias, etc. Additionaly, for mapping genomes to the reference
##   genome.
##  -- Rodrigo Gularte Mérida
##########################################################################

set -e

if [ -f $HOME/bin ]
    then echo "$HOME/bin present"
    else mkdir -p $HOME/bin
         echo "Created $HOME/bin"
         export PATH=$HOME/bin:$PATH
         echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
fi

if [ -f $HOME/src ]
    then echo "$HOME/src present"
    else mkdir -p $HOME/src
         echo "Created $HOME/src"
fi

cd $HOME/src

if [ -x `which git` ]
    then echo "git found, continuing"
    else "'git' not found, please install git to continue http://git-scm.com/downloads "
fi

echo "downloading htslib"
[[ -x htslib ]] || git clone https://github.com/samtools/htslib.git
cd htslib
git pull origin master 
git checkout tags/1.3.2
make
cd ../

echo "downloading samtools"
[[ -x samtools ]] || git clone https://github.com/samtools/samtools.git
cd samtools
git pull origin master
git checkout tags/1.3.1
make
cd ../

echo "downloading bcftools"
[[ -x bcftools ]] || git clone https://github.com/samtools/bcftools.git
cd bcftools
git pull origin master
git checkout tags/1.3.1
make
cd ../

echo "downloading bedtools2"
[[ -x bedtools2 ]] || git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
git pull origin master
git checkout tags/v2.26.0
make
cd ../

echo "downloading bwa"
[[ -x bwa ]] || git clone https://github.com/lh3/bwa.git
cd bwa
git checkout tags/0.7.15
make
cd ../

echo "downloading vcftools"
[[ -x vcftools ]] || git clone https://github.com/vcftools/vcftools.git
cd vcftools
git pull origin master
./autogen.sh
./configure
make
cd ../

echo "downloading bowtie-1"
#wget -O https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download
wget -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download
unzip bowtie.zip
rm bowtie.zip

echo "downloading bowtie2"
#wget -O bowtie2.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip/download
wget -O bowtie2.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip/download
unzip bowtie2.zip
rm bowtie2.zip

echo "downloading tophat2"
#wget -O tophat2.tar.gz http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.14.Linux_x86_64.tar.gz
wget -O tophat2.tar.gz http://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz
tar xvf tophat2.tar.gz
rm tophat2.tar.gz

echo "downloading cufflinks"
wget -O cufflinks.tar.gz http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar xvf cufflinks.tar.gz
rm cufflinks.tar.gz

echo "downloading numpy and scipy"
[[ -x numpy ]] || git clone http://github.com/numpy/numpy.git
cd numpy
git pull origin master
python setup.py install --user
python runtests.py
cd ../

echo "downloading miso"
[[ -x MISO ]] || git clone https://github.com/yarden/MISO.git
cd MISO
git pull origin master
python setup.py install --user
cd ../

echo "downloading matplotlib"
wget -O matplotlib.tar.gz https://pypi.python.org/packages/source/m/matplotlib/matplotlib-1.3.1.tar.gz#md5=789699851de28a543f3d244cf657ff68
tar xvf matplotlib.tar.gz
rm matplotlib.tar.gz
cd matplotlib-1.3.1
python setup.py install --user
cd ../

echo "downloading HTSeq"
wget -O htseq.tar.gz https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1.tar.gz#md5=b7f4f38a9f4278b9b7f948d1efbc1f05
tar xvf htseq.tar.gz
rm htseq.tar.gz
cd HTSeq-0.6.1
python setup.py install --user
cd ../

echo "downloading FastQC"
wget -O fastqc.tar.gz http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip
unzip fastqc.tar.gz
rm fastqc.tar.gz

echo "downloading picard-tools"
#wget -O picard-tools.zip http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip/download
# git clone https://github.com/broadinstitute/picard.git
wget https://github.com/broadinstitute/picard/releases/download/2.7.1/picard.jar
#unzip picard-tools.zip
#rm picard-tools.zip
# cd picard-tools-1.119

echo "downloading ensembl-tools"
[[ -x ensembl-tools ]] || git clone https://github.com/Ensembl/ensembl-tools.git
cd ensembl-tools
git pull origin master
cd ../

echo "downloading freebayes"
[[ -x freebayes ]] || git clone --recursive https://github.com/ekg/freebayes.git
cd freebayes
git pull origin master
make
cd vcflib
make
cd ../../

echo "downloading snpEff"
wget -O snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip

echo "downloading STAR"
[[ -x STAR ]] || git clone https://github.com/alexdobin/STAR.git
cd STAR/
git pull origin master

cd source/
make STAR
cd ../../


echo "downloading Platytpus"
[[ -x Platypus ]] || git clone https://github.com/andyrimmer/Platypus.git
cd Platypus
git pull origin master
export C_INCLUDE_PATH=../htslib/:$C_INCLUDE_PATH
export LIBRARY_PATH=../htslib/:$LIBRARY_PATH
export LD_LIBRARY_PATH=../htslib/:$LD_LIBRARY_PATH
make
cd ../

echo "downloading Trinity assembler, TransDecoder, and Triannotate"
[[ -x trinityrnaseq ]] || git clone https://github.com/trinityrnaseq/trinityrnaseq.git
cd trinityrnaseq
git pull origin master
git checkout tags/v2.0.6
cd ../

[[ -x TransDecoder ]]|| https://github.com/TransDecoder/TransDecoder.git
cd TransDecoder
git pull origin master
git checkout tags/v2.0.1
cd ../

[[ -x Trianotate ]] || https://github.com/Trinotate/Trinotate.git
cd Trianotate
git pull origin master
git checkout tags/v2.0.2
cd ../
