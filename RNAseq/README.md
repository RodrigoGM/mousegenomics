custom file extension
---
<>.slq  : SLURM sbatch file

<>.Rtxt : R session transcript, not a script

<>.awk  : awk command



Read mapping and alignment
---
1) STAR

   star.map.slq : Runs STAR aligner in <RNAseq>.fastq.gz files. In addition it performs a quality assesment using htseq-qa, sorts on coordinate and on name using samtools, and submits output for counting
   
   make.sjdb.sh : Merges the splice junction output from the first run of STAR (star.map.slq) into a single large splice junction db.
   
   star.2pass.genome.slq :  Uses the merged splice junction db and Ensembl Gene annotation to create a STAR genome index for running Engstrom 2pass method.
   
   star.2pass.slq : Runs STAR alinger using the 2pass STAR genome index.
   
   
3) TopHat
   tophat.map.slq : Runs TopHat to map <RNAseq>.fastq.gz files.  In addition it sorts on coordinate and on name using samtools.  Additionally it runs cufflinks on accepted_hits.bam and submits the output for counting using htseq


Read counts
---
1) Cufflinks
   star.cuffX.slq : Indexes Aligned.out.sort.bam using sambamba, runs cufflinks, and cuffquant on the Align.out.sort.bam
   
   cuff.analysis.slq : runs cuffmerge and cuffdiff on outputs of tophat.  
   

2) HTSeq

   star.htseq.slq : Runs htseq-cout to count reads on the output of STAR i.e. Aligned.out.nsort.bam files
   
   idx_htseq.slq : idexes and runs htseq-cout to count reads on the output of Tophat i.e. accepted_hits.bam files
   
   
