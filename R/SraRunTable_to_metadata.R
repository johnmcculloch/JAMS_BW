#' SraRunTable_to_metadata(SraRunTable=NULL, Sample=NULL, discrete=NULL, subsettable=NULL, continuous=NULL)
#'
#' Converts an SRA Run Table to JAMS-style metadata (phenotable, phenolabels).
#' @export
SraRunTable_to_metadata<-function(SraRunTable=NULL, Sample="Run", discrete=NULL, subsettable=NULL, continuous=NULL){
    allvariables <- c(Sample, subsettable, discrete, continuous)
    phenotable <- SraRunTable[ , allvariables]
    phenolabels <- data.frame(Var_label=allvariables, Var_type=c("Sample", rep("subsettable", length(subsettable)), rep("discrete", length(discrete)), rep("continuous", length(continuous))), stringsAsFactors = FALSE)

    metadata <- list()
    metadata[[1]] <- phenotable
    metadata[[2]] <- phenolabels
    names(metadata) <- c("phenotable", "phenolabels")

    return(metadata)
}
