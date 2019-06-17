#' make_antibiogram_experiment(mgseqresfinder=NULL)
#'
#' Agglomerates data from a resfinder mgSeq object to return an antibiogram-style mgSeq object
#' @export

make_antibiogram_experiment <- function(mgseqresfinder = NULL){

    ftt <- fData(mgseqresfinder)
    antibiotics <- colnames(ftt)[!(colnames(ftt) %in% c("Accession", "Description", "Class", "PMID","ResistancePhenotype", "Mechanism"))]
    atbmat <- ftt[,c("Accession", antibiotics)]
    #Transform R into 1 and S into 0
    for (colm in colnames(atbmat)){
        atbmat[,colm][which(atbmat[,colm] == "S")] <- 0
        atbmat[,colm][which(atbmat[,colm] == "R")] <- 1
    }

    countmat <- MRcounts(mgseqresfinder)
    #Remove resfinder_none
    notatbgene <- paste("resfinder", "none")
    rowsToKeep <- which(!(rownames(countmat) %in% notatbgene))
    countmat <- countmat[rowsToKeep, ]
    countmat <- as.data.frame(countmat)
    countmat$Accession <- rownames(countmat)
    atbcounts <- left_join(countmat,atbmat, by = "Accession")
    atbcounts$Accession <- NULL
    atbcounts[] <- lapply(atbcounts, as.numeric)

    Samples <- rownames(pData(mgseqresfinder))
    antibiogram <- as.data.frame(matrix(ncol = length(Samples), nrow = length(antibiotics), data = 0))
    colnames(antibiogram) <- Samples
    rownames(antibiogram) <- antibiotics
    for(sample in Samples){
        antibiogram[,sample] <- unname(sapply(antibiotics, function(x) { sum(atbcounts[,x]*atbcounts[,sample]) } ))
    }

    ptb <- pData(mgseqresfinder)
    cts <- antibiogram[ , Samples]
    atbftt <- data.frame(Accession = rownames(cts), Description = rep("antibiotic", nrow(cts)))
    rownames(atbftt) <- atbftt$Accession

    ##Create metagenomeSeq MRexperiment
    phenotypeData <- AnnotatedDataFrame(ptb)
    ttdata <- AnnotatedDataFrame(atbftt)
    mgseq <- newMRexperiment(cts, phenoData=phenotypeData, featureData=ttdata)

    attr(mgseq, "analysis") <- "antibiogram"

    return(mgseq)
}
