Transcriptome assembly
---
The scripts here are to assemble a transcriptome for a non-reference species.  It is here, maybe a bit out of place, since the mouse is extremely well annotated, however, Trinity RNAseq assembler is used for assembling cancer genomes, and thus potentially useful.  My scripts are not related to cancer research, though take ideas if you find them useful.  There are four directories assembly/ annotation/ mapping/ diffExp/ . The vast majority of the scripts launch utilities and plugins that come with [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki).  In addition, we run [TransDecoder](https://transdecoder.github.io) and [Trinotate](https://trinotate.github.io), all the work to annotate a non-reference genome and carry out a differential expression analyis is being shown.

Other examples of scripts to assemble and annotate a genome are contained in the example/ directory of each of the programs.  The main difference from these script than the runMe.sh in example directories in the three programs is that these are made for a non prepared genomics cluster.  Thus, these scripts are the succesive trial and error as the work progressed.

assembly/
---

* For multiple libraries :
   multiLib.tx.assembly.slq : Script executes a single run of the Trinity assembler combining all the libraries (i.e. *R1.fastq.gz and *R2.fastq.gz) in a single directory containing both R1 and R2 in fastq format (passed as first argument).  Script will format the *R1.fatsq.gz and *R2.fastq.gz contained in a single directory to a single comma separated string.  Optionally, the option to use all the libraries can be changed to use a *R1.csv and *R2.csv files by changing the R1CSV and R2CSV environment variables.

* For individual/single library :
   byLib.tx.assembly.slq : Runs the Trinity assembler in all the libraries (i.e. *R1.fastq.gz and *R2.fastq.gz) contained in a single directory individualy.
   

* Identifying ORF and peptides
   run.transdecoder.sh : runs [TransDecoder](https://transdecoder.github.io)

annotation/
---
This scripts the programs suggested by [Trinotate](https://trinotate.github.io) package in an effort to make some sence out of all the transcripts assembled. Minor modifications were made due to portability and compatibility in various programs, these will be described briefly.  In addition, a sqlite database to contain all the annotation tables is created for querying the results.

assembly.annotation.slq : The script runs a local  [blastx, blastp](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), on the assembled transcripts.fasta and transcripts.pep on SwissProt, Uniref90 and optionally on a TrEMBL protein databases.  In addition, we searched for protein homologs on Pfam using [HMMER](http://hmmer.janelia.org), a signaling peptide prediction with [SignalP](http://www.cbs.dtu.dk/services/SignalP/); a trans membrane prediction with [TmHMM](http://www.cbs.dtu.dk/services/TMHMM/), and an rRNA prediction with [RNAmmer](http://www.cbs.dtu.dk/services/RNAmmer/). 

TrinotateDB.setup.sh : Script to import all the annotation and differential expression tables into a single database.


diffExp/
---

ebseq.slq : launch script for R scripts

ebseq.design.R : generate common data
ebseq.2cond.R  : runs EBTest on a two conditions 

ebseq.SexTRT.Rg.R : runs EBSeq on four comparisons

edgeR.DE.sh : runs the edgeR differential expression script from trinity at the "gene" and "transcript" levels.




custom file extension
---
<>.slq  : SLURM sbatch file

<>.Rtxt : R session transcript, not a script

<>.awk  : awk command


