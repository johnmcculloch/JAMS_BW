#!/bin/bash
accession=$1
module load perl
module load sratoolkit
fastq-dump --split-e "$accession"
splitfile="$accession"_1.fastq
if [ -f "$splitfile" ]
then
    mv "$accession"_1.fastq "$accession"_R1.fastq
    mv "$accession"_2.fastq "$accession"_R2.fastq
    rm "$accession".fastq
else
    mv "$accession".fastq "$accession"_SE.fastq
fi

module unload sratoolkit
module unload perl
