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

echo "downloading htlib"
git clone https://github.com/samtools/htslib.git

echo "downloading samtols"
git clone https://github.com/samtools/samtools.git
cd samtools
git checkout tags/0.1.19
make
cd ../

echo "downloading bedtools2"
git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
git checkout tags/v2.19.0
make
cd ../

echo "downloading bcftools"
git clone https://github.com/samtools/bcftools.git
cd bcftools
make
cd ../

echo "downloading bwa"
git clone https://github.com/lh3/bwa.git
cd bwa
git checkout tags/0.7.7
cd ../

echo "downloading vcftools"
# svn checkout http://svn.code.sf.net/p/vcftools/code/trunk/ vcftools
wget -O vcftools.tar.gz http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12b.tar.gz/download
tar xvf vcftools.tar.gz
rm vcftools.tar.gz
cd vcftools_0.1.12b
make
cd ../


# echo "downloading bowtie-1"
# wget -O bowtie.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip/download
# unzip bowtie.zip
# rm bowtie.zip

echo "downloading bowtie2"
wget -O bowtie2.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip/download
unzip bowtie2.zip
rm bowtie2.zip

echo "downloading tophat2"
wget -O tophat2.tar.gz http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.14.Linux_x86_64.tar.gz
tar xvf tophat2.tar.gz
rm tophat2.tar.gz


echo "downloading cufflinks"
wget -O cufflinks.tar.gz http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar xvf cufflinks.tar.gz
rm cufflinks.tar.gz

echo "downloading miso"
git clone https://github.com/yarden/MISO.git
cd MISO
python setup.py install --user
cd ../

echo "downloading matplotlib"
wget -O matplotlib.tar.gz https://pypi.python.org/packages/source/m/matplotlib/matplotlib-1.3.1.tar.gz#md5=789699851de28a543f3d244cf657ff68
tar xvf matplotlib.tar.gz
rm matplotlib.tar.gz
cd matplotlib-1.3.1
python setup.py install --user
cd ../

echo "downloading numpy and scipy"
git clone http://github.com/numpy/numpy.git
cd numpy
python setup.py install --user
python runtests.py
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
wget -O picard-tools.zip http://sourceforge.net/projects/picard/files/picard-tools/1.119/picard-tools-1.119.zip/download
unzip picard-tools.zip
rm picard-tools.zip
# cd picard-tools-1.119

echo "downloading ensembl-tools"
git clone https://github.com/Ensembl/ensembl-tools.git

echo "downloading freebayes"
git clone --recursive https://github.com/ekg/freebayes.git
cd freebayes 
make
cd ../

echo "downloading snpEff"
wget -O snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip

echo "downloading STAR"
git clone https://github.com/alexdobin/STAR.git
cd STAR
make STAR
cd ../
