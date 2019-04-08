install_biocon_deps<-function(){
bioconductordeps<-c("dada2", "metagenomeSeq", "ComplexHeatmap", "HybridMTest")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(bioconductordeps)
    deps <- c("pigz", "sratoolkit", "trimmomatic", "bowtie2", "megahit", "spades", "kraken2", "convert2bed", "bedtools", "prokka")
}
