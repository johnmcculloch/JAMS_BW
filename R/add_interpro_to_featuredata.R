#' add_interpro_to_featuredata
#'
#' JAMSalpha function
#' @export

add_interpro_to_featuredata <- function(opt = NULL, doinparallel = FALSE){

    #Define which interpro analyses I have
    iproanalyses <- sort(unique(opt$interproscanoutput$Analysis))
    iproanalyses <- c(iproanalyses, "Interpro", "GO")
    flog.info(paste("Found the following Interproscan analyses to harvest:", paste0(iproanalyses, collapse = ", ")))

    if (doinparallel == TRUE){
        iproanalysislist <- get_feature_to_accession_table_parallel(opt = opt)
    } else {
        flog.info("Adding Interproscan analyses signatures to featuredata. Please be patient.")
        #Aggregate accessions serially
        iproanalysislist <- lapply(iproanalyses, function(x) get_feature_to_accession_table(opt = opt, iproanalysis = x))
        names(iproanalysislist) <- iproanalyses
    }

    for (iproanalysis in names(iproanalysislist)){
        opt$featuredata <- left_join(opt$featuredata, iproanalysislist[[iproanalysis]])
        opt$featuredata[, iproanalysis] <- as.character(opt$featuredata[, iproanalysis])
        opt$featuredata[, iproanalysis][is.na(opt$featuredata[, iproanalysis])] <- "none"
    }

    return(opt)
}
