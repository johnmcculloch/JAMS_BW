#' agglomerate_features(ExpObj = NULL, glomby = NULL)
#'
#' Agglomerates features in a JAMS SummarizedExperiment object safely
#' @export

agglomerate_features <- function(ExpObj = NULL, glomby = NULL){

    #Find out what kind of an object it is
    analysis <- metadata(ExpObj)$analysis

    if (analysis == "LKT"){
        ftt <- as.data.frame(rowData(ExpObj))
        targetfeats <- unique(ftt[, glomby])

        cts <- as.matrix(assays(ExpObj)$BaseCounts)

    } else {
        flog.warn(paste("The current version of JAMS is unable to aggregate", analysis, "SummarizedExperiment objects. Sorry for the inconvenience."))
        glomExpObj <- ExpObj
    }

    return(glomExpObj)
}
