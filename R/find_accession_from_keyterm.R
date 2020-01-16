#' find_accession_from_keyterm(ExpObj = NULL, keyterm = NULL)
#'
#' Returns the accession number(s) of a feature(s) matching a general expression present in a SummarizedExperiment or metagenomeSeq object.
#' @export

find_accession_from_keyterm<-function(ExpObj = NULL, keyterm = NULL){

    if (as.character(class(ExpObj)) == "SummarizedExperiment"){
        feattable <- rowData(ExpObj)
    } else {
        feattable <- fData(ExpObj)
    }

    feattable$Feature <- paste(feattable$Accession, feattable$Description, sep="-")
    #Get info by general expression matching
    rowsIwant <- grep(keyterm, feattable$Description, ignore.case = TRUE)
    feature <- feattable$Feature[rowsIwant]
    print(feature)
    accession <- as.vector(feattable$Accession[rowsIwant])

    return(accession)
}
