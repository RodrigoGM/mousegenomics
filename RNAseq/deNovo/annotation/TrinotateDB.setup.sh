#!/bin/bash

TRINOTATE_HOME=/opt/src/Trinotate/
TX_ASSEMBLY=GAL4
DB=boilerplateDB.sqlite
EVALUE=1e-05
PFAM_CUTOFF='DNC'

## initialize db
## Trinotate <db>.sqlite init --gene_trans_map <>.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep transdecoder.pep
$TRINOTATE_HOME/Trinotate $DB init --gene_trans_map $TX_ASSEMBLY.fa.gene_to_trans_map  --transcript_fasta $TX_ASSEMBLY.fa --transdecoder_pep $TX_ASSEMBLY.pep
if [ `echo $?` -eq 0 ] ;  then touch $DB.init.ok ; else  echo "Error" && exit $? ;fi

## load protein hits
## Trinotate <db>.sqlite LOAD_swissprot_blastp blastp.outfmt6
$TRINOTATE_HOME/Trinotate $DB LOAD_swissprot_blastp $TX_ASSEMBLY.pep.sprot.blastp.out
if [ `echo $?` -eq 0 ] ;  then touch $DB.sprot.blastp.ok ; else  echo "Error" && exit $? ;fi

## load transcript hits
## Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
$TRINOTATE_HOME/Trinotate $DB LOAD_swissprot_blastx $TX_ASSEMBLY.fa.sprot.blastx.out
if [ `echo $?` -eq 0 ] ;  then touch $DB.sprot.blastx.ok ; else  echo "Error" && exit $? ;fi

## load protein hits
## Trinotate <db>.sqlite LOAD_trembl_blastp blastp.outfmt6 "uniref90 or trembl"
$TRINOTATE_HOME/Trinotate $DB LOAD_trembl_blastp $TX_ASSEMBLY.pep.u90.blastp.out
if [ `echo $?` -eq 0 ] ;  then touch $DB.u90.blastp.ok ; else  echo "Error" && exit $? ;fi

## load transcript hits
## Trinotate Trinotate.sqlite LOAD__blastx blastx.outfmt6
$TRINOTATE_HOME/Trinotate $DB LOAD_trembl_blastx $TX_ASSEMBLY.fa.u90.blastx.out
if [ `echo $?` -eq 0 ] ;  then touch $DB.u90.blastx.ok ; else  echo "Error" && exit $? ;fi

## load Pfam domain entries
## Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
$TRINOTATE_HOME/Trinotate $DB LOAD_pfam $TX_ASSEMBLY.pep.PFAM.hmm.out
if [ `echo $?` -eq 0 ] ;  then touch $DB.pfam.ok ; else  echo "Error" && exit $? ;fi

## load transmembrane domains
## Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
$TRINOTATE_HOME/Trinotate $DB LOAD_tmhmm $TX_ASSEMBLY.pep.tmhmm.out
if [ `echo $?` -eq 0 ] ;  then touch $DB.tmhmm.ok ; else  echo "Error" && exit $? ;fi

## load signal peptide predictions
##Trinotate Trinotate.sqlite LOAD_signalp signalp.out
$TRINOTATE_HOME/Trinotate $DB LOAD_signalp $TX_ASSEMBLY.pep.signalp.gff
if [ `echo $?` -eq 0 ] ;  then touch $DB.signalp.ok ; else  echo "Error" && exit $? ;fi

## load rRNA predictions
##Trinotate Trinotate.sqlite LOAD_rnammer <file>
$TRINOTATE_HOME/Trinotate $DB LOAD_rnammer $TX_ASSEMBLY.fa.rnammer.gff
if [ `echo $?` -eq 0 ] ;  then touch $DB.rnammer.ok ; else  echo "Error" && exit $? ;fi

## Generate Report -- Note: report is tab delimited, use .txt not .xls
## DNC DGC DTC ; SNC SGC STC = Domain or Sequence cutoffs
##Trinotate Trinotate.sqlite report -E 1e-5 --pfam_cutoff 'DNC|DGC|DTC|SNC|SGC|STC' > trinotate_annotation_report.xls
$TRINOTATE_HOME/Trinotate $DB report -E EVALUE --pfam_cutoff $PFAM_CUTOFF > ${TX_ASSEMBLY}_annotation_${EVALUE}_${PFAM_CUTOFF}.txt
if [ `echo $?` -eq 0 ] ;  then touch $DB.report.ok ; else  echo "Error" && exit $? ;fi

## load edgeR results --transcript_mode
## import the fpkm and DE analysis stuff
## $TRINOTATE_HOME/util/transcript_expression/import_expression_and_DE_results.pl --sqlite Trinotate.sqlite \
##       --samples_file samples_n_reads_described.txt  --count_matrix Trinity_trans.counts.matrix \
##       --fpkm_matrix Trinity_trans.counts.matrix.TMM_normalized.FPKM \
##       --DE_dir edgeR_trans/ < --transcript_mode | --component_mode >

## NOTE:  The above is based on using the transcript-level abundance estimates and corresponding DE analyses.
##        Abundance estimation and DE analysis can be performed at the gene level as well. Optionally, run
##        the above using --component_mode (gene mode) with the corresponding gene-based files. Both
##        transcript-level and gene-level data can be loaded into the single Trinotate.sqlite instance for analysis.

$TRINOTATE_HOME/util/transcript_expression/import_expression_and_DE_results.pl --sqlite $DB \
       --samples_file ../RSEM/conditions$TX_ASSEMBLY.txt --count_matrix ../RSEM/${TX_ASSEMBLY}_genes/genes.counts.matrix \
       --fpkm_matrix ../RSEM/${TX_ASSEMBLY}_isoforms/isoforms.TMM.fpkm.matrix --DE_dir ../RSEM/${TX_ASSEMBLY}_edgeR_tanscripts/ \
       --transcript_mode
if [ `echo $?` -eq 0 ] ;  then touch $DB.edgeR.transcript.ok ; else  echo "Error" && exit $? ;fi

## load edgeR results --component_mode / gene mode
$TRINOTATE_HOME/util/transcript_expression/import_expression_and_DE_results.pl --sqlite $DB \
       --samples_file ../RSEM/conditions$TX_ASSEMBLY.txt --count_matrix ../RSEM/${TX_ASSEMBLY}_genes/genes.counts.matrix \
       --fpkm_matrix ../RSEM/${TX_ASSEMBLY}_genes/genes.TMM.fpkm.matrix --DE_dir ../RSEM/${TX_ASSEMBLY}_edgeR_components/ \
       --component_mode
if [ `echo $?` -eq 0 ] ;  then touch $DB.edgeR.component.ok ; else  echo "Error" && exit $? ;fi

# import transcript clusters
## $TRINOTATE_HOME/util/transcript_expression/import_transcript_clusters.pl \
##    --group_name edgeR_DE_analysis \
##    --analysis_name edgeR_trans/diffExpr.P0.001_C2.matrix.R.all.RData.clusters_fixed_P_20 \
##    --sqlite Trinotate.sqlite \
##    edgeR_trans/diffExpr.P0.001_C2.matrix.R.all.RData.clusters_fixed_P_20/*matrix

$TRINOTATE_HOME/util/transcript_expression/import_transcript_clusters.pl --group_name ${TX_ASSEMBLY}_DE --analysis_name ../RSEM/${TX_ASSEMBLY}_edgeR/$TX_ASSEMBLY.DE.P0.001.FC2x.txt.matrix.RData.clusters_fixed_Kmeans_12 --sqlite $DB  ../RSEM/${TX_ASSEMBLY}_edgeR/$TX_ASSEMBLY.DE.P0.001.FC2x.txt.matrix.RData.clusters_fixed_Kmeans_12/*matrix
if [ `echo $?` -eq 0 ] ;  then touch $DB.edgeR.clusters.ok ; else  echo "Error" && exit $? ;fi

# Reload annotations
#$TRINOTATE_HOME/util/annotation_importer/import_transcript_names.pl Trinotate.sqlite Trinotate_report.xls
$TRINOTATE_HOME/util/annotation_importer/import_transcript_names.pl $DB   ${TX_ASSEMBLY}_annotation_${EVALUE}_${PFAM_CUTOFF}.txt
if [ `echo $?` -eq 0 ] ;  then touch $DB.reload.annot.ok ; else  echo "Error" && exit $? ;fi
## Loading sequence data
$TRINOTATE_HOME/util/trinotateSeqLoader/TrinotateSeqLoader.pl --sqlite $DB --gene_trans_map $TX_ASSEMBLY.fa.gene_to_trans_map --transcript_fasta $TX_ASSEMBLY.fa --transdecoder_pep $TX_ASSEMBLY.pep
if [ `echo $?` -eq 0 ] ;  then touch $DB.sequences.ok ; else  echo "Error" && exit $? ;fi

## cd $TRINOTATE_HOME/TrinotateWeb/ && ./run_mongoose_webserver.sh
