custom file extension
---
<>.slq  : SLURM sbatch file

<>.Rtxt : R session transcript, not a script

<>.awk  : awk command



Read mapping and alignment
---
1) BWA

   bwa.map.slq : Mapping of several bam files using BWA.  Requires SLURM 2.7 or greater as it uses ```--array=<>```


Variant calling
---

1) GATK using HaplotypeCaller and VQSR

   gatk.pp.slq    : Deprecated.  GATK Best Practices on the shell using arrays.
   
   gatk.qpp.scala : Scala script for GATK Best Practices using Queue.jar. Runs In/Del realignment, 2x base recalibration, and performs the variant discovery with the Haplotype Caller.
   
   gatk.qvc.scala : Scala script for GATK Best Practices using Queue.jar. Runs Joint Genotyping of *.g.vcf, runs variant recalibration.

   gatk.qpp.sh    : GATK <>.qpp.scala run shell script 
   
   gatk.qvc.sh	  : GATK <>.qvc.scala run shell script 
   
   
   q.hc.scala     : Scala script to run HaplotypeCaller using Queue.jar
   
   q.hc.sh	  : q.hc.scala run shell script
   

   q.indel.scala  : Scala script to run indelrealinger using Queue.jar
   
   q.indel.sh	  : q.indel.scala run shell script
   
   q.submit.sh	  : batch sumbition of many *bam files to SLURM

2) Freebayes

   freeb.vc.slq   : Runs variant discovery using freebayes-parallel, or optionally a single core freebayes.
   
3) Platypus

   plat.vc.slq    : Runs variant discovery using Platypus
   
4) SAMtools mpileup | bcftools call


Variant filtering, annotation, and processing
---
annotate.vcf.sh	  : Annotates using SnpSift.jar.  Takes <input>.vcf.gz <output>.vcf.gz and <dbsnp>.vcf.gz to annotate. Used to anotate the mouse genomes project vcf, thus all variants with no ID (i.e. ID = ".") are annotated as uk[i++]

cond.sub.awk      : awk script for conditional annotation of ID of a vcf.gz file

