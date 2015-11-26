#!/bin/bash
[[ $# -gt 2 ]] || {
echo "ERROR : requires a Trinity.fasta , a Pfam library, and a CDS training set "
echo "E.g."
echo "sbatch run.transdecoder.sh /path/to/Trinity.fasta ~/Pfam/ /path/to/known_species.cds.all.fa"
echo ""
  exit 1; }

source ../envar.sh

FASTA=$1
PFAM=$2
TRAIN_SET=$3

~/src/TransDecoder/TransDecoder.LongOrfs -t $FASTA
~/src/TransDecoder/TransDecoder.Predict -t Trinity.fasta --retain_pfam_hits $PFAM --train $TRAIN_SET
