#!/bin/sh
cp "$1-utr5.fasta" gene.fasta
perl miRNAsearch.pl +
perl miRNAexplore.pl
DIR=`pwd`
cd ../Deep\ Sequencing/Data/
python miRNAtiling.py $1 > $DIR/tiling.out
cd $DIR
mkdir $1
mv gene.fasta db.fasta rna $1
mv tiling.out $1
rm tmp.fasta
