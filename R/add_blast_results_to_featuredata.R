#' add_interpro_to_featuredata
#'
#' JAMSalpha function
#' @export

add_blast_results_to_featuredata <- function(opt = opt, blastanalyses = NULL){

    if (!(is.null(blastanalyses))){
        #Aggregate accessions serially
        blastanalysislist <- lapply(blastanalyses, function(x) get_feature_to_blast_result_table(opt = opt, blastanalysis = x))
        names(blastanalysislist) <- blastanalyses

        #Redefine blast list to contain only elements with results
        blastanalysislist <- blastanalysislist[sapply(blastanalysislist, function(x) { !(is.null(x)) } )]

        for (blastanalysis in names(blastanalysislist)){
            opt$featuredata <- left_join(opt$featuredata, blastanalysislist[[blastanalysis]], by = "Feature")
            opt$featuredata[, blastanalysis] <- as.character(opt$featuredata[, blastanalysis])
            opt$featuredata[, blastanalysis][is.na(opt$featuredata[, blastanalysis])] <- "none"
        }
    }

    return(opt)
}
