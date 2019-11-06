#' filter_sample_by_class_to_ignore(mgseqobj = NULL, variables = NULL, class_to_ignore = "N_A")
#'
#' Filters a metagenomeSeq object by several criteria.
#' @export

filter_sample_by_class_to_ignore <- function(mgseqobj = NULL, variables = NULL, class_to_ignore = NULL){

    if (!is.null(class_to_ignore)){
        #Start off with all of them
        valid_samples <- rownames(pData(mgseqobj))
        #Keep omitting samples which do not fit the criteria
        for (v in 1:length(variables)){
            valid_samples <-  valid_samples[valid_samples %in% (rownames(pData(obj))[!(pData(mgseqobj)[ , variables[v]] %in% class_to_ignore)])]
        }
        if (length(valid_samples) < 1){
            flog.warn("There are no samples matching the criteria. Returning original object with all samples")
        } else {
            omitted_samples <- rownames(pData(mgseqobj))[!(rownames(pData(mgseqobj)) %in% valid_samples)]
            flog.info(paste("A total of", length(omitted_samples), "samples were omitted for containing", paste0(class_to_ignore, collapse = ", "), "within metadata variables", paste0(variables, collapse = ", ")))
            mgseqobj <- mgseqobj[ , valid_samples]
        }
    }

    #Fix categories of metadata within mgseqobj
    pt <- pData(mgseqobj)
    for (colm in 1:ncol(pt)){
        numtest <- length(which(is.na(as.numeric(pt[, colm])) == TRUE))
        if (numtest == 0){
            pt[, colm] <- as.numeric(pt[, colm])
        }
    }
    pData(mgseqobj) <- pt

    return(mgseqobj)
}
