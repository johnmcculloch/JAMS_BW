#' get_feature_to_blast_result_table
#'
#' JAMSalpha function
#' @export

get_feature_to_blast_result_table <- function(opt = NULL, blastanalysis = NULL){

    #Check first if analysis is available.
    if (blastanalysis %in% names(opt)){
        flog.info(paste("Adding", blastanalysis, "results to featuredata."))
        blastanalysisinterest <- opt[[blastanalysis]]
        featsIwant <- sapply(unique(blastanalysisinterest[ , "Feature"]), function (x) { paste0(sort(unique(blastanalysisinterest[which(blastanalysisinterest[ , "Feature"] == x), "Accession"])), collapse = "|")} )
        feat2acc <- data.frame(Feature = names(featsIwant), Accession = unname(featsIwant), stringsAsFactors = FALSE)
        colnames(feat2acc)[2] <- blastanalysis
    } else {
        flog.info(paste(blastanalysis, "results not found."))
        feat2acc <- NULL
    }

    return(feat2acc)
}
