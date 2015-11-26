Read counts
---
* Cufflinks

   star.cuffX.slq : Indexes Aligned.out.sort.bam using sambamba, runs cufflinks, and cuffquant on the Align.out.sort.bam
   
   cuff.analysis.slq : runs cuffmerge and cuffdiff on outputs of tophat.  
   

* HTSeq

   star.htseq.slq : Runs htseq-cout to count reads on the output of STAR i.e. Aligned.out.nsort.bam files
   
   idx_htseq.slq : idexes and runs htseq-cout to count reads on the output of Tophat i.e. accepted_hits.bam files
   

Differential Expression Analyses
---

* DESeq2

   deseq2.R : runs DESeq2 on multiple files comparing to a single control

   deseq2.plots.R : plots DESeq2 results

DESeq2 works well for two experimental conditions e.g. Treatment vs. Control, however, for Treatment 1 vs. Treatment 2 vs Treatment 3 vs Control | Male or Female , other packages such as EBSeq may be easier to work with. For scripts using EBSeq visit deNovo/diffExp/ in this repository

custom file extension
---
<>.slq  : SLURM sbatch file

<>.Rtxt : R session transcript, not a script

<>.awk  : awk command

