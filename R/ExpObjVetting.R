#' ExpObjVetting(ExpObj = NULL)
#'
#' Performs vetting of a SummarizedExperiment object for use in several functions
#' @export

ExpObjVetting <- function(ExpObj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, glomby = NULL, variables_to_fix = NULL, class_to_ignore = NULL){

        #Get appropriate object to work with
        if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
            stop("This function can only take a SummarizedExperiment object as input.")
        }
        obj <- ExpObj

        if (!(is.null(glomby))){
            obj <- agglomerate_features(ExpObj = obj, glomby = glomby)
        }

        #Exclude samples and features if specified
        if (!(is.null(samplesToKeep))){
            samplesToKeep <- samplesToKeep[samplesToKeep %in% colnames(obj)]
            obj <- obj[, samplesToKeep]
        }

        if (!(is.null(featuresToKeep))){
            featuresToKeep <- featuresToKeep[featuresToKeep %in% rownames(obj)]
            obj <- obj[featuresToKeep, ]
        }

        obj <- suppressWarnings(filter_sample_by_class_to_ignore(SEobj = obj, variables = variables_to_fix, class_to_ignore = class_to_ignore))

    return(obj)
}
