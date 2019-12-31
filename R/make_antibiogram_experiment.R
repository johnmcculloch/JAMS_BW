#' make_antibiogram_experiment(mgseqresfinder = NULL, SEobjresfinder = NULL)
#'
#' Agglomerates data from a resfinder SEobj or mgSeq object to return an antibiogram-style mgSeq or SummarizedExperiment object
#' @export

make_antibiogram_experiment <- function(mgseqresfinder = NULL, SEobjresfinder = NULL){

    if (is.null(SEobjresfinder)){
        ptb <- pData(mgseqresfinder)
        ftt <- fData(mgseqresfinder)
        countmat <- MRcounts(mgseqresfinder)
    } else {
        ptb <- as.data.frame(colData(SEobjresfinder))
        ftt <- as.data.frame(rowData(SEobjresfinder))
        countmat <- as.matrix(assays(SEobjresfinder)$BaseCounts)
    }
    Samples <- rownames(ptb)

    antibiotics <- colnames(ftt)[!(colnames(ftt) %in% c("Accession", "Description", "Class", "PMID","ResistancePhenotype", "Mechanism"))]
    atbmat <- ftt[,c("Accession", antibiotics)]
    #Transform R into 1 and S into 0
    for (colm in colnames(atbmat)){
        atbmat[, colm][which(atbmat[, colm] == "S")] <- 0
        atbmat[, colm][which(atbmat[, colm] == "R")] <- 1
    }

    #Remove resfinder_none
    notatbgene <- paste("resfinder", "none")
    rowsToKeep <- which(!(rownames(countmat) %in% notatbgene))
    countmat <- countmat[rowsToKeep, ]
    countmat <- as.data.frame(countmat)
    countmat$Accession <- rownames(countmat)
    atbcounts <- suppressMessages(left_join(countmat,atbmat, by = "Accession"))
    atbcounts$Accession <- NULL
    atbcounts[] <- lapply(atbcounts, as.numeric)

    antibiogram <- as.data.frame(matrix(ncol = length(Samples), nrow = length(antibiotics), data = 0))
    colnames(antibiogram) <- Samples
    rownames(antibiogram) <- antibiotics
    for(sample in Samples){
        antibiogram[, sample] <- unname(sapply(antibiotics, function(x) { sum(atbcounts[,x] * atbcounts[, sample]) } ))
    }

    cts <- antibiogram[ , Samples]
    cts <- as.matrix(cts)
    atbftt <- data.frame(Accession = rownames(cts), Description = rep("antibiotic", nrow(cts)))
    rownames(atbftt) <- atbftt$Accession

    if (is.null(SEobjresfinder)){
        ##Create metagenomeSeq MRexperiment
        phenotypeData <- AnnotatedDataFrame(ptb)
        ttdata <- AnnotatedDataFrame(atbftt)
        atbobj <- newMRexperiment(cts, phenoData=phenotypeData, featureData=ttdata)
        attr(atbobj, "analysis") <- "antibiogram"
    } else {
        ##Create SummarizedExperiment
        TotalBasesSequenced <- metadata(SEobjresfinder)$TotalBasesSequenced
        assays <- list(cts)
        names(assays) <- "BaseCounts"

        atbobj <- SummarizedExperiment(assays, rowData = atbftt, colData = ptb)
        metadata(atbobj)$TotalBasesSequenced <- TotalBasesSequenced
        attr(atbobj, "analysis") <- "antibiogram"
    }

    return(atbobj)
}
