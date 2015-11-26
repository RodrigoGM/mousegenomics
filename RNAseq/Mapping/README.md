Read mapping and alignment
---
* STAR

   star.map.slq : Runs STAR aligner in <RNAseq>.fastq.gz files. In addition it performs a quality assesment using htseq-qa, sorts on coordinate and on name using samtools, and submits output for counting
   
   make.sjdb.sh : Merges the splice junction output from the first run of STAR (star.map.slq) into a single large splice junction db.
   
   star.2pass.genome.slq :  Uses the merged splice junction db and Ensembl Gene annotation to create a STAR genome index for running Engstrom 2pass method.
   
   star.2pass.slq : Runs STAR alinger using the 2pass STAR genome index.
   
   
* TopHat
   tophat.map.slq : Runs TopHat to map <RNAseq>.fastq.gz files.  In addition it sorts on coordinate and on name using samtools.  Additionally it runs cufflinks on accepted_hits.bam and submits the output for counting using htseq




custom file extension
---
<>.slq  : SLURM sbatch file

<>.Rtxt : R session transcript, not a script

<>.awk  : awk command


