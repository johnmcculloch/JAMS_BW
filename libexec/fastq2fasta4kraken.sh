#!/bin/bash
#Prepares fastq reads for classifying with kraken
#For use with JAMSalpha
#Input is a fastq file. Output is a fasta file with the headers renamed.
#cat "${@}" | sed -n '1~4s/^@/>/p;2~4p' | awk '/^>/{print ">NAss_read_" ++i; next}{print}' > NAss.fasta
cat "${@}" | sed '/^@/!d;s//>/;N' | awk '/^>/{print ">NAss_read_" ++i; next}{print}' > NAss.fasta
