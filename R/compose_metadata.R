#' compose_metadata
#'
#' Returns JAMS-style metadata for a vector of prefixes based on data available in a metadata farm (folder with one or several metadata files)
#' This function is experimental. Use with caution.
#' @export


compose_metadata <- function(metadatafolder = NULL, recursive = FALSE, prefixsubstitutiontable = NULL, restraintoprefixes = NULL){

    require(plyr)
    require(openxlsx)

    #load all xlsx from folder
    metadatafns <- list.files(path = metadatafolder, pattern = ".xlsx$", recursive = recursive)
    list.metadata <- lapply(metadatafns, function (x) { read.xlsx(x) } )

    #Bind metadata and leave missing data with NA
    allmetadata <- rbind.fill(list.metadata)

    if (!is.null(prefixsubstitutiontable)){
        samplerowstosub <- which(allmetadata$Sample %in% prefixsubstitutiontable$OldPrefix == TRUE)
        allmetadata$Sample[samplerowstosub] <- unlist(sapply(samplerowstosub, function (x) { prefixsubstitutiontable$NewPrefix[which(prefixsubstitutiontable$OldPrefix == allmetadata$Sample[x])] }))
    }

    if (!is.null(restraintoprefixes)){
        allmetadata <- subset(allmetadata, Sample %in% restraintoprefixes)
    }

    phenotable <- allmetadata
    phenolabels <- data.frame(Var_label = colnames(phenotable), Var_type = rep("Unknown", ncol(phenotable)), stringsAsFactors = FALSE)
    phenolabels$Var_type[which(phenolabels$Var_label == "Sample")] <- "Sample"

    varstoguess <- which(phenolabels$Var_label != "Sample")

    for (varb in varstoguess){

        colinterest <- phenolabels$Var_label[varb]
        vartypeguess <- class(phenotable[ , which(colnames(phenotable) == colinterest)][!is.na(phenotable[ , which(colnames(phenotable) == colinterest)])])
        if (vartypeguess == "character"){
            vartype <- "discrete"
        } else if (vartypeguess %in% c("numeric", "integer")){
            vartype <- "continuous"
        }
        phenolabels$Var_type[varb] <- vartype
    }

    metadata <- list()
    metadata[[1]] <- phenotable
    metadata[[2]] <- phenolabels
    names(metadata) <- c("phenotable", "phenolabels")

    return(metadata)

}
