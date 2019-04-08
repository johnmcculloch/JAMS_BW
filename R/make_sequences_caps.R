#' make_sequences_caps(sequence=NULL)
#'
#' Transforms sequinr sequences to uppercase
#' @export

make_sequences_caps<-function(sequence=NULL){
    seqnames<-names(sequence)
    sequence<-lapply(1:length(sequence), function(x){ toupper(sequence[[x]]) })
    names(sequence)<-seqnames

    return(sequence)
}
