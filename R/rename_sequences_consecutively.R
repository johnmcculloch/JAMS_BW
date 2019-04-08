#' rename_sequences_consecutively(sequence=NULL, sequence=NULL, headerprefix="ctg")
#'
#' Renames sequences numerically consecutively
#' @export

rename_sequences_consecutively<-function(sequence=NULL, headerprefix="ctg"){
    sequence_names<-paste(headerprefix, formatC(1:length(sequence), width=nchar(length(sequence)), flag="0"), sep="_")
 
    #Fix attribute when loaded using sequinr load.fasta with set.attributes being TRUE and ensure CAPS
    for(t in 1:length(sequence)){
        tmps<-sequence[[t]]
        attr(tmps, "name")<-sequence_names[t]
        tmps<-toupper(tmps)
        sequence[[t]]<-tmps
    }
    names(sequence)<-sequence_names

    return(sequence)
}
