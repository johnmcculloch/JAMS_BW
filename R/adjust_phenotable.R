#Prepare pheno object for metagenomeSeq package
#' adjust_phenotable(opt = NULL, list.data = NULL, addtaxlevelstoisolates = NULL)
#'
#' Adjusts the metadata table into a format which can be used by the MetagenomeSeq package.
#' @export

adjust_phenotable <- function(opt = NULL, list.data = NULL, addtaxlevelstoisolates = NULL, class_to_ignore = "N_A"){

    #flog.info("Adjusting phenotable classes by type of variable.")
    varlist <- define_kinds_of_variables(phenolabels = opt$phenolabels, phenotable = opt$phenotable, maxclass = 10000, maxsubclass = 10000, class_to_ignore = class_to_ignore, verbose = FALSE)

    #Rename column to Sample, if it isn't already
    colmtosub <- which(colnames(opt$phenotable) == varlist$sample)
    colnames(opt$phenotable)[colmtosub] <- "Sample"
    #make unique if there was already another column called Sample
    if (length(which(colnames(opt$phenotable) == "Sample")) > 1){
        colnames(opt$phenotable) <- make.unique(colnames(opt$phenotable), sep = "_")
        colnames(opt$phenotable)[colmtosub] <- "Sample"
        opt$phenotable$Var_label <- colnames(phenotable)
    }
    rownames(opt$phenotable) <- opt$phenotable$Sample

    Samples <- rownames(opt$phenotable)

    #add information regarding sample type
    if (!is.null(list.data)){
        projdata <- as.data.frame(matrix(data = "unknown", nrow = (length(Samples)), ncol = 3))
        projdata[] <- lapply(projdata, as.character)
        colnames(projdata) <- c("Sample", "JAMS_Run_type", "JAMS_Process", "JAMS_Kdb_Version")
        projdata$Sample <- Samples

        #Fetch data pertaining to each sample
        for (s in 1:length(Samples)){
            projstats <- NULL
            projstats <- list.data[[paste(Samples[s], "projinfo", sep = "_")]]
            rownames(projstats) <- projstats$Run_info
            projstats <- projstats[c("Run_type", "Process"), ]
            projstats[] <- lapply(projstats, as.character)
            projdata[which(projdata$Sample == Samples[s]), which(colnames(projdata) == "JAMS_Run_type")] <- projstats["Run_type", "Run_value"]
            projdata[which(projdata$Sample == Samples[s]), which(colnames(projdata) == "JAMS_Process")] <- projstats["Process", "Run_value"]
            projdata[which(projdata$Sample == Samples[s]), which(colnames(projdata) == "JAMS_Kdb_Version")] <- projstats["Process", "Run_value"]
        }

        opt$phenotable <- left_join(opt$phenotable, projdata, by = "Sample")
        phenolabels_projinfo <- data.frame(Var_label = c("JAMS_Run_type", "JAMS_Process", "JAMS_Kdb_Version"), Var_type = c("discrete", "discrete", "discrete"), stringsAsFactors = FALSE)
        opt$phenolabels <- rbind(opt$phenolabels, phenolabels_projinfo)

        #add taxonomic information to isolates, if there are any
        if ((!is.null(addtaxlevelstoisolates)) && length(which(projdata$JAMS_Run_type == "isolate") > 0)){
            taxdata <- as.data.frame(matrix(data = "not_isolate", nrow = (length(Samples)), ncol = (length(addtaxlevelstoisolates) + 1 )))
            taxdata[] <- lapply(taxdata, as.character)
            newtaxlvlnames <- paste("Isolate", addtaxlevelstoisolates, sep = "_")
            colnames(taxdata) <- c("Sample", newtaxlvlnames)
            taxdata$Sample <- Samples
            isolatesamples <- taxdata$Sample[which(projdata$JAMS_Run_type == "isolate")]
            for (i in 1:length(isolatesamples)){
                taxstats <- NULL
                taxstats <- list.data[[paste(isolatesamples[i], "LKTdose", sep = "_")]]
                taxstats <- taxstats[, c(addtaxlevelstoisolates, "NumBases")]
                taxstats[] <- lapply(taxstats, as.character)
                taxstats$NumBases <- as.numeric(taxstats$NumBases)
                taxdata[which(taxdata$Sample == isolatesamples[i]), newtaxlvlnames] <- as.character(as.matrix(unname(taxstats[which(taxstats$NumBases == max(taxstats$NumBases)), addtaxlevelstoisolates])))
            }

            opt$phenotable <- left_join(opt$phenotable, taxdata, by = "Sample")
            phenolabels_taxon <- data.frame(Var_label = newtaxlvlnames, Var_type = rep("discrete", length(newtaxlvlnames)), stringsAsFactors = FALSE)
            opt$phenolabels <- rbind(opt$phenolabels, phenolabels_taxon)
        }
    }

    #Add data from readdata if desired
    if (!(is.null(opt$readdata))){
        dfr <- opt$readdata
        rownames(dfr) <- dfr$Sample
        dfr <- dfr[opt$phenotable$Sample, ]
        dfr$GbNAHS <- round((dfr$NonHost_bases / 1000000000), 2)
        dfr$GbTrim <- round((dfr$Trim_bases / 1000000000), 2)
        opt$phenotable$GbNAHS <- rep(0, nrow(opt$phenotable)) #Account for the fact that pheno Samples may be missing in readdata
        opt$phenotable$PctAss <- rep(0, nrow(opt$phenotable)) #Account for the fact that pheno Samples may be missing in readdata
        opt$phenotable$GbNAHS <- dfr$GbNAHS[match(opt$phenotable$Sample, dfr$Sample)]
        opt$phenotable$PctAss <- dfr$PctAss[match(opt$phenotable$Sample, dfr$Sample)]
    }

    rownames(opt$phenotable) <- opt$phenotable$Sample

    return(opt)
}
