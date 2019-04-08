#' filter_sequence_by_length(sequence=NULL, minlength=0, maxlength=Inf)
#'
#' Filters sequences by size
#' @export

filter_sequence_by_length<-function(sequence=NULL, minlength=0, maxlength=Inf){
    sequence_lengths<-lapply(1:length(sequence), function(x){ length(sequence[[x]]) })
    seqsIwant<-which(sequence_lengths > minlength & sequence_lengths < maxlength)
    sequence<-sequence[seqsIwant]

    return(sequence)
}
