#' export_expvec_to_XL(expvec = NULL, usefulexp = NULL, filename = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, includemetadata = TRUE, returncounts = FALSE)
#'
#' Exports counts and featuredata in a SummarizedExperiment vector into a single spreadsheet.
#' @export

export_expvec_to_XL <- function(expvec = NULL, usefulexp = NULL, filename = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, includemetadata = TRUE, returncounts = FALSE){

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

    if (any(c((asPPM == TRUE), !is.null(applyfilters), !is.null(featcutoff), !is.null(GenomeCompletenessCutoff), !is.null(PctFromCtgscutoff)))){
        flog.info("Counts will be exported as PPM")
        countunits <- "PPM"
    } else {
        flog.info("Counts will be exported as raw number of bases")
        countunits <- "NumBases"
    }

    #Get counts
    for (x in 1:length(expvec2)){
        analysis <- metadata(expvec2[[x]])$analysis
        flog.info(paste("Exporting", analysis))

        presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

        exp_filt <- filter_experiment(ExpObj = expvec2[[analysis]], featcutoff = presetlist$featcutoff, asPPM = asPPM, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

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

    flog.info(paste("Saving spreadsheet as", filename))
    write.xlsx(countvec, file = filename, asTable = TRUE, rowNames = TRUE, colNames = TRUE, borders = "all", colWidths = "auto")
    if (returncounts){
        return(countvec)
    }
}
