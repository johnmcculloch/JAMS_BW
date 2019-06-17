#!/bin/bash
#Install JAMS dependencies
module load gcc/9.1.0
module load java
module load python
module load perl
cd $HOME
brew update
brew tap brewsci/bio
for JAMSdep in git wget pigz sratoolkit samtools trimmomatic bowtie2 kraken2 bedtools bedops prokka
do
    brew ls --versions "$JAMSdep" && brew upgrade "$JAMSdep" || brew install "$JAMSdep"
done

for JAMSassem in megahit spades
do
    brew install --ignore-dependencies "$JAMSassem"
done

module unload gcc
module load R
