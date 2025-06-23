#' filter_sample_by_class_to_ignore(SEobj = NULL, variables = NULL, class_to_ignore = NULL)
#'
#' Filters a SummarizedExperiment object by several criteria.
#' @export

filter_sample_by_class_to_ignore <- function(SEobj = NULL, variables = NULL, class_to_ignore = NULL){

    ptb <- as.data.frame(colData(SEobj))
    allmetadata <- metadata(SEobj)
    Samples <- rownames(ptb)

    if (!is.null(class_to_ignore)){
        #Start off with all of them
        valid_samples <- rownames(ptb)
        variables <- colnames(ptb)[colnames(ptb) %in% variables]
        #Keep omitting samples which do not fit the criteria, as long as variables are contained in current metadata
        if (all(c((!is.null(variables)), (length(variables) > 0)))){
            for (v in 1:length(variables)){
                valid_samples <-  valid_samples[valid_samples %in% (rownames(ptb)[!(ptb[ , variables[v]] %in% class_to_ignore)])]
            }
        }
        if (length(valid_samples) < 1){
            flog.warn("There are no samples matching the criteria. Returning original object with all samples")
        } else {
            omitted_samples <- rownames(ptb)[!(rownames(ptb) %in% valid_samples)]
            flog.info(paste("A total of", length(omitted_samples), "samples were omitted for containing", paste0(class_to_ignore, collapse = ", "), "within metadata variables", paste0(variables, collapse = ", ")))
            SEobj <- SEobj[ , valid_samples]
            ptb <- as.data.frame(colData(SEobj))
        }
    }

    #Fix categories of metadata within Experiment object
    for (colm in 1:ncol(ptb)){
        numtest <- length(which(is.na(as.numeric(ptb[, colm])) == TRUE))
        if (numtest == 0){
            ptb[, colm] <- as.numeric(ptb[, colm])
        }
    }

    ExpObj <- SummarizedExperiment(assays = assays(SEobj), rowData = rowData(SEobj), colData = ptb)
    metadata(ExpObj) <- allmetadata

    return(ExpObj)
}
