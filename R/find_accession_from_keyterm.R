#' find_accession_from_keyterm(ExpObj = NULL, keyterm = NULL)
#'
#' Returns the accession number(s) of a feature(s) matching a general expression present in a SummarizedExperiment or metagenomeSeq object.
#' @export

find_accession_from_keyterm <- function(ExpObj = NULL, keyterm = NULL){

    feattable <- rowData(ExpObj)
    accessionvector <- NULL
    feattable$Feature <- paste(feattable$Accession, feattable$Description, sep="-")

    #Get info by general expression matching
    for (kt in keyterm){
        rowsIwant <- grep(kt, feattable$Description, ignore.case = TRUE)
        features <- feattable$Feature[rowsIwant]
        print(features)
        accessionvector <- c(accessionvector, as.vector(feattable$Accession[rowsIwant]))
    }

    return(accessionvector)
}
