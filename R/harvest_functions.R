#' harvest_functions
#'
#' JAMSalpha function
#' @export

harvest_functions<-function(opt = opt, noninterproanalyses = c("FeatType", "ECNumber", "Product", "resfinder", "plasmidfinder", "napdos", "serofinderH", "serofinderO", "vfdb", "abricate")){

    data(ECdescmap)
    data(GOtermdict)
    data(MetaCycAccession2Description)

    flog.info("Harvesting functional data.")
    #Harvest common functions
    basicanalyses <- noninterproanalyses[noninterproanalyses %in% colnames(opt$featuredata)]

    featurenumbaseslist <- lapply(basicanalyses, function(x) { compute_signature_numbases(featuredata = opt$featuredata, columnname =x ) })
    featurenumbaseslist <- plyr::ldply(featurenumbaseslist, rbind)
    featurenumbaseslist$`.id` <- NULL
    rownames(featurenumbaseslist) <- featurenumbaseslist$Accession

    #Harvest interpro functions, if applicable
    if(opt$skipipro != TRUE){
        opt <- fix_interproscanoutput(opt=opt)
    }
    interpronumbaseslist <- NULL
    if("interproscanoutput" %in% names(opt)){
        opt <- add_interpro_to_featuredata(opt = opt, doinparallel = TRUE)

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

        flog.info("Creating counts table for Interpro signatures.")
        interpronumbaseslist <- lapply(iproanalyses, function(x) { compute_signature_numbases(featuredata = opt$featuredata, columnname = x) })
        names(interpronumbaseslist) <- iproanalyses
        interpronumbaseslist <- plyr::ldply(interpronumbaseslist, rbind)
        interpronumbaseslist$`.id` <- NULL
        rownames(interpronumbaseslist) <- interpronumbaseslist$Accession
    }

    opt$featuredose <- rbind(featurenumbaseslist, interpronumbaseslist)

    #Add descriptions to featuredose
    LKTcols <- colnames(opt$featuredose)[4:ncol(opt$featuredose)]
    opt$featuredose$Description <- rep("none", nrow(opt$featuredose))
    #rearrange
    opt$featuredose <- opt$featuredose[, c("Analysis", "Accession", "Description", "NumBases", LKTcols)]
    #Add EC numbers
    descriptions <- ECdescmap$Description[match(opt$featuredose$Accession, ECdescmap$Accession)]
    opt$featuredose$Description[which(!(is.na(descriptions)))] <- descriptions[which(!(is.na(descriptions)))]

    #Add interpro descriptions, if applicable
    if("interproscanoutput" %in% names(opt)){
        #Add GO descriptions
        descriptions <- GOtermdict$Description[match(opt$featuredose$Accession, GOtermdict$Accession)]
        opt$featuredose$Description[which(!(is.na(descriptions)))] <- descriptions[which(!(is.na(descriptions)))]
        #Add MetaCyc descriptions
        MetaCycdescriptions <- MetaCycAccession2Description$Description[match(opt$featuredose$Accession, MetaCycAccession2Description$Accession)]
        opt$featuredose$Description[which(!(is.na(MetaCycdescriptions)))] <- MetaCycdescriptions[which(!(is.na(MetaCycdescriptions)))]

        #Make description dictionary
        dictaccessions <- opt$interproscanoutput$Accession
        dictdescriptions <- opt$interproscanoutput$Description
        dictaccessions <- append(dictaccessions, opt$interproscanoutput$IproAcc, after = length(dictaccessions))
        dictdescriptions <- append(dictdescriptions, opt$interproscanoutput$IproDesc, after = length(dictdescriptions))
        acc2desc <- data.frame(Accession = dictaccessions, Description = dictdescriptions, stringsAsFactors = FALSE)
        acc2desc <- acc2desc[!(duplicated(acc2desc)), ]
        #Add interpro descriptions to featuredose
        descriptions <- acc2desc$Description[match(opt$featuredose$Accession, acc2desc$Accession)]
        opt$featuredose$Description[which(!(is.na(descriptions)))] <- descriptions[which(!(is.na(descriptions)))]
    }

    return(opt)
}
