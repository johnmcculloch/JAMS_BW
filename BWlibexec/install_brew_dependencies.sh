#!/bin/bash
#Install JAMS dependencies
module load gcc/8.2.0
module load java
module load python
module load perl
cd $HOME
brew update
brew tap brewsci/bio
brew install git wget pigz sratoolkit samtools trimmomatic bowtie2 kraken2 bedtools bedops prokka
brew install --ignore-dependencies megahit spades
module unload gcc
module load R
