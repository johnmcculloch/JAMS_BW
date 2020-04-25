#' add_interpro_to_featuredata
#'
#' JAMSalpha function
#' @export

add_interpro_to_featuredata <- function(opt = NULL, doinparallel = FALSE){

    #Define which interpro analyses I have
    iproanalyses <- sort(unique(opt$interproscanoutput$Analysis))
    if ("IproAcc" %in% colnames(opt$interproscanoutput)){
        iproanalyses <- c(iproanalyses, "Interpro")
    }
    if ("GOterms" %in% colnames(opt$interproscanoutput)){
        iproanalyses <- c(iproanalyses, "GO")
    }
    if ("Pathways" %in% colnames(opt$interproscanoutput)){
        iproanalyses <- c(iproanalyses, "MetaCyc")
    }
    flog.info(paste("Found the following Interproscan analyses to harvest:", paste0(iproanalyses, collapse = ", ")))

    if (doinparallel == TRUE){

        appropriatenumcores <- max(1 , (min((opt$threads - 2), length(iproanalyses))))
        flog.info(paste("Adding Interproscan analyses signatures to featuredata with", appropriatenumcores, "CPUs. Please be patient."))
        iproanalysislist <- mclapply(iproanalyses, function (x) { get_feature_to_accession_table(opt = opt, iproanalysis = x) }, mc.cores = appropriatenumcores)
        names(iproanalysislist) <- iproanalyses

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
