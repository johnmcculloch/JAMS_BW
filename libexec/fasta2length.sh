#!/bin/bash
#Prepares fastq reads for classifying with kraken
#For use with JAMSalpha
#Input is a fasta file. Output is a TSV of sequence names and their length.
cat "${@}" | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
