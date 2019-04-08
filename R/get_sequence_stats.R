#' get_contig_stats
#'
#' JAMSalpha function
#' @export

get_sequence_stats<-function(sequences=NULL){
    sequencenames<-names(sequences)
    sequencelengths<-getLength(sequences)
    sequenceGCs<-sapply(1:length(sequences), function(x){ GC(sequences[[x]]) })
    
    seqstats<-data.frame(Sequence=sequencenames, Length=sequencelengths, GC=sequenceGCs)

    return(seqstats)
}
