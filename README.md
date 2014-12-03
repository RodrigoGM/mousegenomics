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


directories and scripts
-----------------

DNAseq/
---
Contains scripts for the mapping and alignment of DNA high throughput sequencing
data using bwa.  In addition it contains script for variant calling using
four pipelines
1) GATK using HaplotypeCaller and VQSR
2) Freebayes
3) Platypus
4) SAMtools mpileup | bcftools call

RNAseq/
---
Contains scripts for the mapping and alignment of RNA  high throughput sequencing
data using STAR aligner and TopHat.  In addition it contains scripts for counting
and quality assesment of the RNA sequences.

QTLMapping/
---
Contains scripts for QTL mapping in mouse inbred strains.

GeneExpression/
---
Contains scripts for the analysis of gene expression data.