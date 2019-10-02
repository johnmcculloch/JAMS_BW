#' randomize_experiment(mgseqobj = NULL, newname = "Noise", shuffleclassifieds = FALSE, shufflemetadata = TRUE, iterations = 10)
#'
#' Shuffles data (counts and metadata) of a metagenomeSeq MRexperiment object. Useful for sanity checking.
#' @export

randomize_experiment <- function (mgseqobj = NULL, newname = "Noise", shuffleclassifieds = FALSE, shufflemetadata = TRUE, iterations = 10){

    analysis <- attr(mgseqobj, "analysis")

    #Get counts matrix and metadata
    cts <- MRcounts(mgseqobj)
    pdata <- pData(mgseqobj)

    #Define shuffling functions
    shuff_row <- function (x = NULL, matr = NULL){
        shuffvec <- sample(matr[x, ], replace = TRUE)

        return(shuffvec)
    }

    shuff_colm <- function (x = NULL, matr = NULL){
        shuffvec <- sample(matr[ , x], replace = TRUE)

        return(shuffvec)
    }

    flog.info(paste("Randomizing", analysis, "MRexperiment"))

    for (itn in 1:iterations){
        flog.info(paste0("Shuffling matrix horizontally and vertically ", itn, "/", iterations, " times."))

        #Start by shuffling sample values by row
        for (feat in 1:nrow(cts)){
            cts[feat, ] <- shuff_row(x = feat, matr = cts)
        }

        #Shuffle row values by sample now. But don't touch unclassifieds, if requested.
        if (shuffleclassifieds == FALSE){
            dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
            classifiedrows <- which(!(rownames(cts) %in% dunno))
            unclasscts <- cts[(rownames(cts) %in% dunno), ]
            if(class(unclasscts) == "numeric"){
                unclasscts <- t(unclasscts)
                rownames(unclasscts) <- rownames(cts)[(rownames(cts) %in% dunno)]
            }
            cts <- cts[classifiedrows, ]
        }

        for (colm in 1:ncol(cts)){
            cts[ , colm] <- shuff_colm(x = colm, matr = cts)
        }

        if (shuffleclassifieds == FALSE){
            cts <- rbind(cts, unclasscts)
        }

        if (shufflemetadata == TRUE){
            for (colm in 1:ncol(pdata)){
                pdata[ , colm] <- sample(pdata[, colm], replace=FALSE)
                rownames(pdata) <- pdata$Sample
            }
        }
    }

    #Rebuild metagenomeseq experiment and rename it.
    sampleorder <- rownames(pdata)
    ftt <- fData(mgseqobj)
    ftt <- ftt[rownames(cts), ]
    cts <- cts[, sampleorder]

    ##Create metagenomeSeq MRexperiment
    phenotypeData <- AnnotatedDataFrame(pdata)
    ttdata <- AnnotatedDataFrame(ftt)
    randmgseqobj <- newMRexperiment(cts, phenoData = phenotypeData, featureData = ttdata)

    attr(randmgseqobj, "analysis") <- newname

    return(randmgseqobj)
}
