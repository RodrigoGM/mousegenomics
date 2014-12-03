mouseomics
==========

This repository contains code to help with the setup of a Genomics cluster without
 root privilages.  The information detailed here, as the name implies, is oriented
 to Mouse Genomics.  Thus, most scripts will be downloading mouse data from 
 Ensembl, UCSC or NCBI.

Software set up should be generic, though check to make sure you are downloading
  appropriate reference files for you own species of interest (not everyone is 
  interested in humans).  Pay special attention to the $PATH where things are
  located and that you're not deleting something important.

Analysis pipelines will most likely need to be adapted to your own style and 
  data structure. Though... take the ideas and best of luck with your work !

custom file extension
---
<>.slq  : SLURM sbatch file

<>.Rtxt : R session transcript, not a script

<>.awk  : awk command



Directories/
---

### DNAseq/

Contains scripts for the mapping and alignment of DNA high throughput sequencing
data using bwa.  In addition it contains script for variant calling using
four pipelines

1) GATK using HaplotypeCaller and VQSR

2) Freebayes

3) Platypus

4) SAMtools mpileup | bcftools call

### RNAseq/

Contains scripts for the mapping and alignment of RNA  high throughput sequencing
data using STAR aligner and TopHat.  In addition it contains scripts for counting
and quality assesment of the RNA sequences.

### QTLMapping/

Contains scripts for QTL mapping in mouse inbred strains.

### GeneExpression/

Contains scripts for the analysis of gene expression data.

## Scripts

create.igenome.sh
  This script creates an iGenome « LIKE » reference genome.  The iGenomes are Illumina's bundled referece and annotation packages to be used with GenomeStudio.  an iGenome has the following simplified directory tree structure
  
    Mus_musculus/                                                       ## Organisms/species
    └── Ensembl                                                 ## Format/Data source
        ├── Grcm38                                                      ## Reference Genome release
        │   ├── Annotation                                              ## Feature annotation file
        │   │   ├── Genes
        │   │   │   └── genes.gtf, refFlat.txt, refGene.txt             ## gene annotation files
        │   │   ├── SmallRNA
        │   │   │   └── mature.fa precursor.fa                  ## Fasta files from miRBase
        │   │   └── Variation
        │   │       └── Mus_musculus.gvf                                ## Genomic variation
        │   └── Sequence
        │       ├── AbundantSequences                           ## fasta files
        │       │   └── adapter_contam1.fa* MT.fa* musRibosomal.fa*
        |           |       phix.fa* polyA.fa* polyC.fa*                ## * refers to .fa.2bpb and .fa.vld index
        │       ├── Bowtie2Index
        │       │   └── genome.*.bt2 genome.fa -> ../path/to/ref        ## Bowtie2 indexes and link to genome.fa
        │       ├── BowtieIndex
        │       │   └── genome.*.ebwt genome.fa -> ../path/to/ref       ## Bowtie1 indexes and link to genome.fa
        │       ├── BWAIndex
        │       │   └── genome.fa.* genome.fa -> ../path/to/ref ## BWA indexes and link to genome.fa
        │       ├── Chromosomes
        │       │   └── 1.fa ... 19.fa X.fa Y.fa MT.fa          ## Fasta files from individual chromosomes
        │       ├── Squashed-Mus_musculus-Ensembl-Grcm38
        └──     └── WholeGenomeFasta
                └── genome.fa genome.fa.fai genome.dict         ## Genome Reference in Fasta, samtools faidx
                                                                ## and picard CreateGenomeDictionary.jar index
                                                                ## files
   
   The directory structured is simplified as in addition to Annotation and Sequence, directories, there is an additional GenomeStudio which contains symbolic to the appropriate locations, and each directory has an Archives or versionXX directory where the data is contained and then sym-linked to the files shown here.  The symbolic link strategy was not reconstructed, however, the essence of the data structure was preserved in this script shuch that it resembles the above.
The location and additional README files will be downloaded along with the data for reference.
