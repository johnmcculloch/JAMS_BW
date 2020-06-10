#' add_interpro_to_featuredata
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
        #subset ipro to contain only applicable analysis
        if (!(iproanalysis %in% c("Interpro", "GO", "MetaCyc"))){
            iprointerest <- subset(opt$interproscanoutput, Analysis == iproanalysis)
            accessioncol <- "Accession"
            descriptioncol <- "Description"
        } else if (iproanalysis == "MetaCyc"){
            iprointerest <- subset(opt$interproscanoutput, Pathways != "none")
            accessioncol <- "Accession"
            descriptioncol <- "Description"
        } else {
            iprointerest <- subset(opt$interproscanoutput, IproAcc != "none")
            accessioncol <- "IproAcc"
            descriptioncol <- "IproDesc"
        }

        #Declare useful funcitons
        fish_pathway <- function (FeatInterest = NULL, iprointerest = NULL, pathwayspace = NULL) {
            nr_annots <- unique(unlist(strsplit(iprointerest[which(iprointerest[ , "Feature"] == FeatInterest), "Pathways"], "|", fixed = TRUE)))
            nr_annots <- nr_annots[grep(pathwayspace, nr_annots)]
            if (length(nr_annots) < 1){
                nr_annots <- "none"
            }
            nr_annots <- gsub("MetaCyc: ", "", nr_annots)
            nr_annots <- paste0(sort(unique(nr_annots)), collapse = "|")
            return(nr_annots)
        }

        if (iproanalysis == "GO"){
            #get rid of information without GO terms
            iprointerest <- subset(iprointerest, GOterms != "none")
            #If looking for GO terms, split up accessions in GOterms column, and get sorted, non-reundant list.
            featsIwant <- sapply(unique(iprointerest[ , "Feature"]), function (x) { paste0(sort(unique(unlist(strsplit(iprointerest[which(iprointerest[ , "Feature"] == x), "GOterms"], "|", fixed = TRUE)))), collapse = "|")} )
        } else if (iproanalysis == "MetaCyc"){
             #get rid of information without GO terms
            iprointerest <- subset(iprointerest, Pathways != "")
            #If looking for pathway terms, split up annotations in Pathways column, and get sorted, non-reundant list.
            featsIwant <- sapply(unique(iprointerest[ , "Feature"]), function (x) { fish_pathway(FeatInterest = x, iprointerest = iprointerest, pathwayspace = iproanalysis) } )
        } else {
            featsIwant <- sapply(unique(iprointerest[, "Feature"]), function (x) { paste0(sort(unique(iprointerest[which(iprointerest[,"Feature"] == x), accessioncol])), collapse = "|")} )
        }

    feat2acc <- data.frame(Feature = names(featsIwant), Accession = unname(featsIwant), stringsAsFactors = FALSE)
    #Remove any non-informative information
    feat2acc <- feat2acc[which(!(feat2acc$Accession %in% c("none", ""))), ]
    colnames(feat2acc)[2] <- iproanalysis

    return(feat2acc)
}


#' harvest_functions
#'
#' JAMSalpha function
#' @export

harvest_functions <- function(opt = opt, noninterproanalyses = c("FeatType", "ECNumber", "Product", "resfinder", "plasmidfinder", "napdos", "serofinderH", "serofinderO", "vfdb", "abricate")){

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
        opt <- fix_interproscanoutput(opt = opt)
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


#' compute_signature_numbases
#'
#' JAMSalpha function
#' @export

compute_signature_numbases <- function (featuredata = NULL, columnname = NULL){
    numbasesdf <- tidyr::separate_rows(featuredata, columnname, sep = fixed("\\|")) %>% dplyr::group_by(get(columnname),LKT) %>% dplyr::summarise(NumBases = sum(as.integer(NumBases)))
    colnames(numbasesdf) <- c("Accession", "LKT", "NumBases")
    numbasesdf <- numbasesdf[ , c("Accession", "LKT", "NumBases")]
    numbasesdf <- numbasesdf %>% dplyr::arrange(-NumBases) %>% tidyr::spread(LKT,NumBases)
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


#' fix_interproscanoutput
#'
#' JAMSalpha function
#' @export

fix_interproscanoutput <- function(opt = NULL){

    if ("iprodir" %in% names(opt)){
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

        flog.info("Harvesting and integrating Interproscan data.")

        interprotsvs<-file.path(opt$iprodir, list.files(path=opt$iprodir, pattern=".tsv"))
        #load tsvs into a single object in memory
        readipro<-function(x){
            read.table(file=x, sep="\t", header=FALSE, quote="", skipNul=FALSE, fill=TRUE, colClasses = "character", col.names=c("Feature","MD5","AALength","Analysis","Accession","Description","Start","Stop","Score","Status","Date","IproAcc","IproDesc","GOterms","Pathways"))
        }
        alliprotsvs <- lapply(interprotsvs, readipro)
        ipro<-plyr::ldply(alliprotsvs, rbind)

        #Clean-up datafrane
        ipro<-subset(ipro, Analysis !="MobiDBLite")
        ipro$MD5<-NULL
        ipro$Start<-NULL
        ipro$Stop<-NULL
        ipro$Start<-NULL
        ipro$Score<-as.numeric(ipro$Score)
        ipro$Score[is.na(ipro$Score)]<- 0
        ipro<-subset(ipro, Score < 0.001)
        ipro<-subset(ipro, Status == "T")
        ipro$Score<-NULL
        ipro$Status<-NULL
        ipro$Date<-NULL
        ipro$AALength<-as.numeric(ipro$AALength)

        #Fill in the blanks
        ipro <- ipro %>% mutate(Description = ifelse(Description == "", "none", Description))
        ipro <- ipro %>% mutate(IproAcc = ifelse(IproAcc == "", "none", IproAcc))
        ipro <- ipro %>% mutate(IproDesc = ifelse(IproDesc == "", "none", IproDesc))
        ipro <- ipro %>% mutate(GOterms = ifelse(GOterms == "", "none", GOterms))
        ipro <- ipro %>% mutate(Pathways = ifelse(Pathways == "", "none", Pathways))

        #remove eventual duplicates
        ipro <- ipro[!duplicated(ipro), ]
        opt$interproscanoutput<-ipro
    }

    return(opt)
}
