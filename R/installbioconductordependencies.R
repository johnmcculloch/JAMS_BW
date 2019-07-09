install_biocon_deps <- function(){
bioconductordeps <- c("dada2", "metagenomeSeq", "ComplexHeatmap", "HybridMTest", "genefilter")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install(bioconductordeps)
}
