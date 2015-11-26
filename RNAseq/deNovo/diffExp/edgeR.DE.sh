#!/bin/bash

set -e -x

TRINOTATE_HOME=/opt/src/Trinotate/
TRINITY_HOME=/opt/src/trinityrnaseq/

## conditionsMMU4.txt
## tabl delimited file containg :
## <TRT1> <sample1>
## <TRT1> <sample2>
## <TRT1> ...
## <TRT1> <sampleN>
## <TRT..N> ...
## <CTRL> <sample1>
## <CTRL> <sample2>
## <CTRL> ...
## <CTRL> <sample1>

## run on RSEM/
$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix genes/genes.counts.matrix --method edgeR --samples_file conditionsMMU4tx.txt --min_rowSum_counts 100 --output edgeR_components/

$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix isoforms/isoforms.counts.matrix --method edgeR --samples_file conditionsMMU4tx.txt --min_rowSum_counts 100 --output edgeR_transcripts/

