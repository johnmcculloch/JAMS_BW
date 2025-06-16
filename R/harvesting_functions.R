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
        iproanalysislist <- lapply(iproanalyses, function(x) { get_feature_to_accession_table(opt = opt, iproanalysis = x)} )
        names(iproanalysislist) <- iproanalyses
    }

    for (iproanalysis in names(iproanalysislist)){
        opt$featuredata <- left_join(opt$featuredata, iproanalysislist[[iproanalysis]], by = "Feature")
        opt$featuredata[, iproanalysis] <- as.character(opt$featuredata[, iproanalysis])
        opt$featuredata[, iproanalysis][is.na(opt$featuredata[, iproanalysis])] <- "none"
    }

    return(opt)
}


#' add_blast_results_to_featuredata
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


#' get_feature_to_accession_table
#'
#' JAMSalpha function
#' @export

get_feature_to_accession_table <- function(opt = NULL, iproanalysis = NULL){

    flog.info(paste("Adding", iproanalysis, "signatures to featuredata."))
    if (iproanalysis == "GO"){
        #get rid of information without GO terms
        iprointerest <- subset(opt$interproscanoutput, !(GOterms %in% c("none", "", "-")))[ , c("Feature", "GOterms")]
        #Make non-redundant data frame
        iprointerest <- tidyr::separate_rows(iprointerest, all_of("GOterms"), sep = fixed("\\|"))
        feat2acc <- iprointerest %>% group_by(Feature) %>% summarize(Accession = str_c(GOterms, collapse = "|"))
        feat2acc <- as.data.frame(feat2acc)

    } else if (iproanalysis == "MetaCyc"){
        #get rid of information without Pathways terms
        iprointerest <- subset(opt$interproscanoutput, !(Pathways %in% c("none", "", "-")))[ , c("Feature", "Pathways")]
        #Make non-redundant data frame
        iprointerest <- tidyr::separate_rows(iprointerest, all_of("Pathways"), sep = fixed("\\|"))
        iprointerest <- iprointerest[grep("MetaCyc", iprointerest$Pathways), ]
        iprointerest$Pathways <- gsub("MetaCyc: ", "", iprointerest$Pathways)
        feat2acc <- iprointerest %>% group_by(Feature) %>% summarize(Accession = str_c(Pathways, collapse = "|"))
        feat2acc <- as.data.frame(feat2acc)
    } else if (iproanalysis == "Interpro"){
        iprointerest <- subset(opt$interproscanoutput, !(IproAcc %in% c("none", "", "-")))[ , c("Feature", "IproAcc")]
        iprointerest <- tidyr::separate_rows(iprointerest, all_of("IproAcc"), sep = fixed("\\|"))
        feat2acc <- iprointerest %>% group_by(Feature) %>% summarize(Accession = str_c(IproAcc, collapse = "|"))
        feat2acc <- as.data.frame(feat2acc)
    } else {
        #Get appropriate analysis space and get rid of information without Accession terms
        iprointerest <- subset(opt$interproscanoutput, Analysis == iproanalysis)
        iprointerest <- subset(iprointerest, !(Accession %in% c("none", "", "-")))[ , c("Feature", "Accession")]
        #Make non-redundant data frame
        iprointerest <- tidyr::separate_rows(iprointerest, all_of("Accession"), sep = fixed("\\|"))
        feat2acc <- iprointerest %>% group_by(Feature) %>% summarize(Accession = str_c(Accession, collapse = "|"))
        feat2acc <- as.data.frame(feat2acc)
    }

    #Remove any non-informative information
    feat2acc <- feat2acc[which(!(feat2acc$Accession %in% c("none", "", "-"))), ]
    colnames(feat2acc)[2] <- iproanalysis

    return(feat2acc)
}


#' harvest_functions
#'
#' JAMSalpha function
#' @export

harvest_functions <- function(opt = opt, noninterproanalyses = c("FeatType", "ECNumber", "Product", "resfinder", "plasmidfinder", "napdos", "serofinderH", "serofinderO", "vfdb", "abricate"), taxonomic_spaces = c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin"), doinparallel = TRUE, check_ipro_jobs_status = TRUE){

    data(ECdescmap)
    data(GOtermdict)
    #data(MetaCycAccession2Description)

    flog.info("Harvesting functional data.")

    #Due to alternative taxonomic spaces in JAMS2, a feature abundance table for each different taxonomic space will be generated and exported to a named list called featuredoses (note the plural).
    featuredoses <- list()
    interprodoses <- list()

    valid_taxonomic_spaces <- taxonomic_spaces[taxonomic_spaces %in% colnames(opt$featuredata)]
    flog.info(paste("Taxonomic spaces available for stratifying functional features by taxa are", paste0(valid_taxonomic_spaces, collapse = ", ")))

    #Harvest common functions
    basicanalyses <- noninterproanalyses[noninterproanalyses %in% colnames(opt$featuredata)]
    for (taxonomic_space in valid_taxonomic_spaces){
        featurenumbaseslist <- lapply(basicanalyses, function(x) { compute_signature_numbases(featuredata = opt$featuredata, taxonomic_space = taxonomic_space, columnname = x ) })
        featurenumbaseslist <- plyr::ldply(featurenumbaseslist, rbind)
        featurenumbaseslist$`.id` <- NULL
        rownames(featurenumbaseslist) <- featurenumbaseslist$Accession
        featuredoses[[taxonomic_space]] <- featurenumbaseslist
    }

    #Harvest interpro functions, if applicable
    if (opt$skipipro != TRUE){
        opt <- fix_interproscanoutput(opt = opt, check_ipro_jobs_status = check_ipro_jobs_status)
    }
    interpronumbaseslist <- NULL

    if ("interproscanoutput" %in% names(opt)){

        opt <- add_interpro_to_featuredata(opt = opt, doinparallel = FALSE)

        iproanalyses <- sort(unique(opt$interproscanoutput$Analysis))
        if ("IproAcc" %in% colnames(opt$interproscanoutput)){
            iproanalyses <- c(iproanalyses, "Interpro")
        }
        if ("GOterms" %in% colnames(opt$interproscanoutput)){
            iproanalyses <- c(iproanalyses, "GO")
        }
        #Pathways were phased out
        #if ("Pathways" %in% colnames(opt$interproscanoutput)){
        #    iproanalyses <- c(iproanalyses, "MetaCyc")
        #}

        flog.info("Creating counts table for Interpro signatures.")

        for (taxonomic_space in valid_taxonomic_spaces){

            interpronumbaseslist <- lapply(iproanalyses, function(x) { compute_signature_numbases(featuredata = opt$featuredata, taxonomic_space = taxonomic_space, columnname = x) })
            names(interpronumbaseslist) <- iproanalyses
            interpronumbaseslist <- plyr::ldply(interpronumbaseslist, rbind)
            interpronumbaseslist$`.id` <- NULL

            #Remove a "-" which is annoyingly being added as an accession for GO and Interpro.
            interpronumbaseslist <- subset(interpronumbaseslist, Accession != "-")

            rownames(interpronumbaseslist) <- interpronumbaseslist$Accession
            interprodoses[[taxonomic_space]] <- interpronumbaseslist
        }
    }

    for (taxonomic_space in valid_taxonomic_spaces){
        ##Add descriptions to featuredose
        Taxoncols <- colnames(featuredoses[[taxonomic_space]])[4:ncol(featuredoses[[taxonomic_space]])]

        #Add EC numbers
        featuredoses[[taxonomic_space]] <- left_join(featuredoses[[taxonomic_space]], ECdescmap[ , c("Accession", "Description")], by = "Accession")

        #rearrange
        featuredoses[[taxonomic_space]] <- featuredoses[[taxonomic_space]][, c("Analysis", "Accession", "Description", "NumBases", Taxoncols)]

        #Add interpro descriptions, if applicable
        if ("interproscanoutput" %in% names(opt)){
 
            #Make description dictionary from the interproscan output
            dictaccessions <- opt$interproscanoutput$Accession
            dictdescriptions <- opt$interproscanoutput$Description
            dictaccessions <- append(dictaccessions, opt$interproscanoutput$IproAcc, after = length(dictaccessions))
            dictdescriptions <- append(dictdescriptions, opt$interproscanoutput$IproDesc, after = length(dictdescriptions))
            acc2desc <- data.frame(Accession = dictaccessions, Description = dictdescriptions, stringsAsFactors = FALSE)
            acc2desc <- acc2desc[!(duplicated(acc2desc)), ]
            #Add GO descriptions
            acc2desc <- rbind(acc2desc, GOtermdict[ , colnames(acc2desc)])
            #Add interpro descriptions to featuredose
            interprodoses[[taxonomic_space]] <- left_join(interprodoses[[taxonomic_space]], acc2desc[ , c("Accession", "Description")], by = "Accession")
            #rearrange
            interprodoses[[taxonomic_space]] <- interprodoses[[taxonomic_space]][ , c("Analysis", "Accession", "Description", "NumBases", Taxoncols)]
            #Bank to featuredoses
            featuredoses[[taxonomic_space]] <- rbind(featuredoses[[taxonomic_space]], interprodoses[[taxonomic_space]])
        }

        #Fix NAs and rownames
        featuredoses[[taxonomic_space]]$Description[which(is.na(featuredoses[[taxonomic_space]]$Description))] <- "none"
        rownames(featuredoses[[taxonomic_space]]) <- featuredoses[[taxonomic_space]]$Accession
    }

    opt$abundances$functional <- featuredoses

    return(opt)
}


#' compute_signature_numbases
#'
#' JAMSalpha function
#' @export

compute_signature_numbases <- function (featuredata = NULL, columnname = NULL, blastanalyses = c("abricate", "resfinder", "plasmidfinder", "vfdb"), taxonomic_space = "LKT"){

    if (columnname %in% blastanalyses){
        numbasesdf <- featuredata %>% group_by(across(all_of(c(columnname, taxonomic_space)))) %>% dplyr::summarise(NumBases = sum(as.integer(NumBases)), .groups = "keep")
    } else {
        numbasesdf <- tidyr::separate_rows(featuredata, all_of(columnname), sep = fixed("\\|")) %>% group_by(across(all_of(c(columnname, taxonomic_space)))) %>% dplyr::summarise(NumBases = sum(as.integer(NumBases)), .groups = "keep")
    }
    colnames(numbasesdf) <- c("Accession", "Taxon", "NumBases")
    numbasesdf <- numbasesdf[ , c("Accession", "Taxon", "NumBases")]
    numbasesdf <- numbasesdf %>% dplyr::arrange(-NumBases) %>% pivot_wider(names_from = Taxon, values_from = NumBases, values_fill = 0)
    numbasesdf <- as.data.frame(numbasesdf)
    numbasesdf[is.na(numbasesdf)] <- 0
    taxa <- colnames(numbasesdf)[2:ncol(numbasesdf)]

    if(length(taxa) > 1){
        numbasesdf$NumBases <- rowSums(numbasesdf[ , 2:ncol(numbasesdf)])
    } else {
        numbasesdf$NumBases <- (numbasesdf[ , taxa])
    }
    numbasesdf$Analysis <- rep(columnname, nrow(numbasesdf))
    numbasesdf <- numbasesdf[ , c("Analysis", "Accession", "NumBases", taxa)]
    numbasesdf$Accession[(which(numbasesdf$Accession == "none"))] <- rep(paste(columnname, "none", sep = "_"), length(which(numbasesdf$Accession == "none")))

    return(numbasesdf)
}


#' fix_interproscanoutput(opt = NULL, check_ipro_jobs_status = TRUE)
#'
#' JAMSalpha function
#' @export

fix_interproscanoutput <- function(opt = NULL, check_ipro_jobs_status = TRUE){

    if ("iprodir" %in% names(opt)){

        if (check_ipro_jobs_status){
            #Get interprojob
            opt$iprojob <- system2('cat', args=file.path(opt$iprodir, "ipro.job"), stdout = TRUE, stderr = FALSE)
            #See if job finished
            iprojobstatus <- system2('sacct', args=c("-j", opt$iprojob), stdout = TRUE)
            #eliminate header
            iprojobstatus <- iprojobstatus[3:length(iprojobstatus)]
            iprojobstatus <- iprojobstatus[grep("quick", iprojobstatus)]
            totaljobs <- length(iprojobstatus)
            completedjobs <- length(grep("COMPLETED", iprojobstatus))
            ratiojobscomplete <- completedjobs/totaljobs
            runningjobs <- length(grep("RUNNING", iprojobstatus))

            #Delay if there still are jobs to complete
            nattempt <- 1
            if (opt$analysis %in% c("metagenome", "metatranscriptome")){
                requirediprocompleteness <- 0.9
            } else {
                requirediprocompleteness <- 0.97
            }

            #Wait until jobs are close to complete
            while ((ratiojobscomplete < requirediprocompleteness) && nattempt < 50) {
                flog.info("Interproscan analysis of proteome is still incomplete.")
                flog.info(paste("There are", runningjobs, " Interpro jobs running for this sample."))
                flog.info(paste("The proportion of Interpro jobs complete is currently", round(ratiojobscomplete, 2)))
                flog.info(paste("Will check again in 5 minutes time. This is attempt", nattempt, "of 50 before giving up."))
                Sys.sleep(300)
                nattempt <- nattempt + 1
                #See if job finished
                iprojobstatus <- system2('sacct', args=c("-j", opt$iprojob), stdout = TRUE)
                #eliminate header
                iprojobstatus <- iprojobstatus[3:length(iprojobstatus)]
                iprojobstatus <- iprojobstatus[grep("quick", iprojobstatus)]
                totaljobs <- length(iprojobstatus)
                completedjobs <- length(grep("COMPLETED", iprojobstatus))
                ratiojobscomplete <- completedjobs/totaljobs
                runningjobs <- length(grep("RUNNING", iprojobstatus))
            }
        }

        flog.info("Harvesting and integrating Interproscan data.")

        interprotsvs <- file.path(opt$iprodir, list.files(path = opt$iprodir, pattern=".tsv"))
        appropriatenumcores <- max(1 , (min((opt$threads - 2), length(interprotsvs))))

        #load tsvs into a single object in memory
        readipro <- function(x){
            fread(file=x, sep="\t", header=FALSE, quote="", fill=TRUE, colClasses = "character", col.names=c("Feature","MD5","AALength","Analysis","Accession","Description","Start","Stop","Score","Status","Date","IproAcc","IproDesc","GOterms","Pathways"))
        }

        alliprotsvs <- mclapply(interprotsvs, readipro, mc.cores = appropriatenumcores)
        ipro <- plyr::ldply(alliprotsvs, rbind)
        ipro_spaces_to_keep <- c("CDD", "Pfam", "SUPERFAMILY", "NCBIfam", "ProSitePatterns", "ProSiteProfiles", "SMART")

        #Clean-up datafrane
        ipro <- subset(ipro, Analysis %in% ipro_spaces_to_keep)
        ipro$MD5 <- NULL
        ipro$Start <- NULL
        ipro$Stop <- NULL
        ipro$Start <- NULL
        #ipro$Score[is.na(ipro$Score)] <- 0
        #ipro$Score <- as.numeric(ipro$Score)
        ipro <- subset(ipro, Score <= 0.001)
        ipro <- subset(ipro, Status == "T")
        ipro$Score <- NULL
        ipro$Status <- NULL
        ipro$Date <- NULL
        ipro$AALength <- as.numeric(ipro$AALength)

        #Phase out pathways because there ater way too many entries, likely false positives
        ipro$Pathways <- NULL

        #Fill in the blanks
        ipro <- ipro %>% mutate(Description = ifelse(Description == "", "No_description", Description))
        ipro <- ipro %>% mutate(IproAcc = ifelse(IproAcc == "", "none", IproAcc))
        ipro <- ipro %>% mutate(IproDesc = ifelse(IproDesc == "", "No_description", IproDesc))
        ipro <- ipro %>% mutate(GOterms = ifelse(GOterms == "", "none", GOterms))
        #ipro <- ipro %>% mutate(Pathways = ifelse(Pathways == "", "none", Pathways))
        ipro$GOterms <- gsub("^-$", "none", ipro$GOterms)
        ipro$IproAcc <- gsub("^-$", "none", ipro$IproAcc)
        ipro[which(ipro$IproAcc == "none"), "IproDesc"] <- "No_description"

        #Remove the "(InterPro)" tag added to the new interproscan output. Need to remove to match to GO dictionary.
        ipro$GOterms <- gsub("\\(InterPro\\)", "", ipro$GOterms)

        #remove eventual duplicates
        ipro <- ipro[!duplicated(ipro), ]
        opt$interproscanoutput <- ipro

    }

    return(opt)
}
