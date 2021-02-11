#!/bin/bash
accession=$1
module load perl
module load sratoolkit
fastq-dump --split-e "$accession"
mv "$accession"_1.fastq "$accession"_R1.fastq
mv "$accession"_2.fastq "$accession"_R2.fastq
rm "$accession".fastq
module unload sratoolkit
module unload perl
