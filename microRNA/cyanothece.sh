#!/bin/sh
blastn -query db.fasta -subject ../Genome/circular.fasta -word_size 7 -gapopen 4 -gapextend 2 -reward 2 -penalty -2
