#' filter_sequence_by_name(input_sequences=NULL, sequencenames=NULL, keep=TRUE)
#'
#' Filters sequinr sequences by name either keeping or discarding specified sequence names
#' @export

filter_sequence_by_name<-function(input_sequences=NULL, sequencenames=NULL, keep=TRUE){
    input_seqnames<-names(input_sequences)

    if(!(is.null(sequencenames))){
        #Filter by read name
        if(keep==TRUE){
            seqtoKeep<-which(input_seqnames %in% sequencenames)
            #Check that all sequences you asked for are available
            #if(length(input_sequences) != length(seqtoKeep)){
                #print("WARNING: Number of sequences requested is different to number of sequences available in input file.")
            #}
        } else {
            seqtoKeep<-which(!(input_seqnames %in% sequencenames))
            #Check that all sequences you asked to eliminate existed in the input file.
            if(!(all(sequencenames %in% input_seqnames))){
                print("WARNING: Not all sequences you requested to exclude were available in input file.")
            }
        }
        output_sequences<-input_sequences[seqtoKeep]
    } else {
        output_sequences<-input_sequences
    }

    output_sequences<-make_sequences_caps(output_sequences)

    return(output_sequences)
}

