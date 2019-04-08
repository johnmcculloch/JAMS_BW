#' find_accession_from_keyterm(mgseqobj=NULL, keyterm=NULL)
#'
#' Returns the accession number(s) of a feature(s) matching a general expression present in a metagenomeSeq object. 
#' @export

find_accession_from_keyterm<-function(mgseqobj=NULL, keyterm=NULL){
    obj<-mgseqobj

    #Make dictionary
    feattable<-fData(obj)
    feattable$Feature<-paste(feattable$Accession, feattable$Description, sep="-")
    #Get info by general expression matching
    rowsIwant<-grep(keyterm, feattable$Description, ignore.case = TRUE)
    feature<-feattable$Feature[rowsIwant]
    print(feature)
    accession<-as.vector(feattable$Accession[rowsIwant])
    
    return(accession)
}
