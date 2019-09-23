#' make_readdata(Samples=NULL, list.data=NULL)
#'
#' This function makes a data frame with the sequencing statistics.
#' @export

make_readdata <- function(Samples = NULL, list.data = NULL){
    #Consider all samples if list not given
    if (missing(Samples)){
        Samples <- name_samples(list.data = list.data)
    }

    #Create empty data frame to hold all read data
    readdata <- as.data.frame(matrix(data = 0, nrow = (length(Samples)), ncol = 6))
    colnames(readdata) <- c("Sample", "Raw_bases", "Trim_bases", "NonHost_bases", "Subsampled_bases", "Assembled_bases", "Proj_type")
    readdata$Sample <- Samples

    #Fetch data pertaining to each sample
    for (s in 1:length(Samples)){
        readstats <- NULL
        projstats <- NULL
        #Try and find reads stats for that sample.
        readstats <- list.data[[paste(Samples[s], "readstats", sep="_")]]
        if (is.null(readstats)){
            readstats <- data.frame(Read_type=c("Raw", "Trimmed", "NonHost", "Subsampled", "Assembled"), Base_counts=c(0,0,0,0), stringsAsFactors = FALSE)
        }
        rownames(readstats) <- readstats$Read_type
        readstats$Read_type <- NULL

        projstats <- list.data[[paste(Samples[s], "projinfo", sep="_")]]
        rownames(projstats) <- projstats$Run_info

        readstats[] <- lapply(readstats, as.numeric)
        readdata[which(readdata$Sample == Samples[s]), which(colnames(readdata) == "Raw_bases")] <- readstats["Raw", "Base_counts"]
        readdata[which(readdata$Sample == Samples[s]), which(colnames(readdata) == "Trim_bases")] <- readstats["Trimmed", "Base_counts"]
        readdata[which(readdata$Sample == Samples[s]), which(colnames(readdata) == "NonHost_bases")] <- readstats["NonHost", "Base_counts"]
        readdata[which(readdata$Sample == Samples[s]), which(colnames(readdata) == "Subsampled")] <- readstats["Subsampled", "Base_counts"]
        readdata[which(readdata$Sample == Samples[s]), which(colnames(readdata) == "Assembled_bases")] <- readstats["Assembled", "Base_counts"]
        readdata[which(readdata$Sample == Samples[s]), which(colnames(readdata) == "Proj_type")] <- as.character(projstats["Process", "Run_value"])
    }

    #If there are samples without Host information (i.e. Host == none), then consider non-host as trimmed.
    if (length(readdata$NonHost_bases[which(is.na(readdata$NonHost_bases) == TRUE)]) > 0){
        #Report the occurrence
        print(paste("Samples", paste0(readdata$Sample[which(is.na(readdata$NonHost_bases) == TRUE)], collapse = ", "), "are missing non-host data. Will consider that as being Trimmed bases."))
        readdata$NonHost_bases[which(is.na(readdata$NonHost_bases) == TRUE)] <- readdata$Trim_bases[which(is.na(readdata$NonHost_bases) == TRUE)]
    }

    #Calculate statistics
    readdata$PctQCpass <- round((readdata$Trim_bases/readdata$Raw_bases) * 100, 1)
    readdata$PctHost <- 100 - (round((readdata$NonHost_bases/readdata$Trim_bases) * 100, 1))
    readdata$PctAss <- (round((readdata$Assembled_bases/readdata$NonHost_bases) * 100, 1))

    #Adjust for samples which were not assembled
    notassembledsamples <- subset(readdata, Proj_type != "Assemble_from_reads")[]$Sample
    readdata$PctQCpass[which(readdata$Sample %in% notassembledsamples)] <- rep(100, length(notassembledsamples))
    readdata$PctHost[which(readdata$Sample %in% notassembledsamples)] <- rep(0, length(notassembledsamples))
    readdata$PctAss[which(readdata$Sample %in% notassembledsamples)] <- rep(100, length(notassembledsamples))

    #Protect against values which are not numeric
    readdata[is.na(readdata)] <- 0
    rownames(readdata) <- readdata$Sample

    return(readdata)
}
