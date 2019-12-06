#' export_expvec_to_XL(expvec = NULL, usefulexp = NULL, filename = NULL, asPPM = TRUE, asPA = FALSE, mgSeqnorm = FALSE, includemetadata = TRUE)
#'
#' Exports counts and featuredata in an experiment vector into a single spreadsheet.
#' @export

export_expvec_to_XL <- function(expvec = NULL, usefulexp = NULL, filename = NULL, asPPM = TRUE, asPA = FALSE, mgSeqnorm = FALSE, includemetadata = TRUE){

    if (!(is.null(usefulexp))){
        usefulexp <- names(expvec)[(names(expvec) %in% usefulexp)]
        expvec2 <- expvec[usefulexp]
    } else {
        expvec2 <- expvec
    }

    pt <- as.data.frame(pData(expvec2[[1]]))

    countvec <- NULL
    countvec <- vector("list", length = 1000)
    cvn <- 1
    if (includemetadata == TRUE){
        print("Exporting metadata.")
        countvec[[cvn]] <- pt
        names(countvec)[cvn] <- "Metadata"
        cvn <- cvn + 1
    }

    if (asPPM == TRUE){
        countunits <- "PPM"
    } else {
        countunits <- "NumBases"
    }

    if (asPA == TRUE) {
        countunits <- "Presence_Absence"
    }

    #Get counts
    for (x in 1:length(expvec2)){
        print(paste("Exporting", names(expvec2)[x]))

        exp_filt <- filter_experiment(mgseqobj = expvec2[[x]], asPA = asPA, asPPM = asPPM, mgSeqnorm = mgSeqnorm)
        cts <- MRcounts(exp_filt)
        #Eliminate empty rows
        ctsnotzero <- which(rowSums(cts) > 0)
        cts <- cts[ctsnotzero, ]
        ctsname <- paste(names(expvec2)[x], countunits, sep="_")
        countvec[[cvn]] <- as.data.frame(cts)
        names(countvec)[cvn] <- ctsname
        cvn <- cvn + 1

        feattbl <- as.data.frame(fData(expvec2[[x]]))
        ftsname <- paste(names(expvec2)[x], "featuretable", sep="_")
        countvec[[cvn]] <- feattbl
        names(countvec)[cvn] <- ftsname
        cvn <- cvn + 1
    }

    countvec <- countvec[sapply(countvec, function(x){ !(is.null(x)) })]

    print(paste("Saving spreadsheet as", filename))
    write.xlsx(countvec, file = filename, asTable = TRUE, rowNames = TRUE, colNames = TRUE, borders = "surrounding", colWidths = "auto")

}
