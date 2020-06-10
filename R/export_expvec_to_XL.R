#' export_expvec_to_XL(expvec = NULL, usefulexp = NULL, filename = NULL, asPPM = TRUE, asPA = FALSE, PPM_normalize_to_bases_sequenced = FALSE, includemetadata = TRUE)
#'
#' Exports counts and featuredata in a SummarizedExperiment vector into a single spreadsheet.
#' @export

export_expvec_to_XL <- function(expvec = NULL, usefulexp = NULL, filename = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, includemetadata = TRUE){

    if (!(is.null(usefulexp))){
        usefulexp <- names(expvec)[(names(expvec) %in% usefulexp)]
        expvec2 <- expvec[usefulexp]
    } else {
        expvec2 <- expvec
    }

    pt <- as.data.frame(colData(expvec2[[1]]))

    countvec <- NULL
    countvec <- vector("list", length = 1000)
    cvn <- 1
    if (includemetadata == TRUE){
        flog.info("Exporting metadata.")
        countvec[[cvn]] <- pt
        names(countvec)[cvn] <- "Metadata"
        cvn <- cvn + 1
    }

    if (asPPM == TRUE){
        countunits <- "PPM"
    } else {
        countunits <- "NumBases"
    }

    #Get counts
    for (x in 1:length(expvec2)){
        flog.info(paste("Exporting", names(expvec2)[x]))
        exp_filt <- filter_experiment(ExpObj = expvec2[[x]], asPPM = asPPM, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced)
        cts <- assays(exp_filt)$BaseCounts
        ctsname <- paste(names(expvec2)[x], countunits, sep="_")
        countvec[[cvn]] <- as.data.frame(cts)
        names(countvec)[cvn] <- ctsname
        cvn <- cvn + 1

        if ("GenomeCompleteness" %in% names(assays((expvec2)[[x]]))){
            cts <- assays(exp_filt)$GenomeCompleteness
            ctsname <- paste(names(expvec2)[x], "GenomeCompleteness", sep="_")
            countvec[[cvn]] <- as.data.frame(cts)
            names(countvec)[cvn] <- ctsname
            cvn <- cvn + 1
        }

        if ("PctFromCtgs" %in% names(assays((expvec2)[[x]]))){
            cts <- assays(exp_filt)$PctFromCtgs
            ctsname <- paste(names(expvec2)[x], "PctFromCtgs", sep="_")
            countvec[[cvn]] <- as.data.frame(cts)
            names(countvec)[cvn] <- ctsname
            cvn <- cvn + 1
        }

        feattbl <- as.data.frame(rowData(expvec2[[x]]))
        ftsname <- paste(names(expvec2)[x], "featuretable", sep="_")
        countvec[[cvn]] <- feattbl
        names(countvec)[cvn] <- ftsname
        cvn <- cvn + 1
    }

    countvec <- countvec[sapply(countvec, function(x){ !(is.null(x)) })]

    print(paste("Saving spreadsheet as", filename))
    write.xlsx(countvec, file = filename, asTable = TRUE, rowNames = TRUE, colNames = TRUE, borders = "all", colWidths = "auto")
}
