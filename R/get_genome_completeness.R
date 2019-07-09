#' get_genome_completeness<-function(pheno=pheno, list.data=NULL))
#'
#' Returns a data frame with the estimated number of complete genomes of given taxa in each sample
#' @export


get_genome_completeness <- function(pheno = pheno, list.data = NULL){
    #Get data for features
    Samples <- rownames(pheno)
    assobjects <- paste(Samples, "assemblystats", sep="_")
    assdoses <- list.data[assobjects]
    names(assdoses) <- Samples

    assall <- bind_rows(assdoses, .id = "id")
    assall[is.na(assall)] <- 0
    colnames(assall)[1] <- "Sample"
    assall <- assall[ , c("Sample", "Taxon", "ProbNumGenomes")]

    genomecompletenessdf <- assall %>% spread(Sample, ProbNumGenomes)
    genomecompletenessdf[is.na(genomecompletenessdf)] <- 0

    rownames(genomecompletenessdf) <- genomecompletenessdf$Taxon
    genomecompletenessdf$Taxon <- NULL

    return(genomecompletenessdf)
}
