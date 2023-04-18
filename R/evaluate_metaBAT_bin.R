#' evaluate_metaBAT_bin(opt=NULL)
#'
#' Calculates the genome completeness and assembly stats for metaBAT bins.
#' @export

evaluate_metaBAT_bin <- function(opt = NULL){

    flog.info("Evaluating quality of MAGs.")

    ######################################
    ##   Declare functions that must    ##
    ## be declared within this funciton ##
    ######################################

    getN50 <- function(contigs2length){
        quantile50 <- (sum(as.numeric(contigs2length$Length))) * 0.5
        N50 <- max(subset(contigs2length, Lcum > quantile50)$Length)
        return(N50)
    }

    getL50 <- function(contigs2length){
        quantile50 <- (sum(as.numeric(contigs2length$Length))) * 0.5
        L50 <- length(subset(contigs2length, Lcum > quantile50)$Length)
        return(L50)
    }

    getN90 <- function(contigs2length){
        quantile90 <- (sum(as.numeric(contigs2length$Length))) * 0.9
        N90 <- max(subset(contigs2length, Lcum > quantile90)$Length)
        return(N90)
    }

    getcontigs2length <- function(opt = NULL, bininterest = NULL){
        contigsdata <- opt$contigsdata
        contigs2length <- subset(contigsdata, contigsdata[ , "MetaBATbin"] == bininterest)
        contigs2length <- contigs2length[, c("MetaBATbin", "Contig", "Length", "LKT")]
        contigs2length <- contigs2length[order(contigs2length$Length, decreasing = TRUE), ]
        contigs2length$Lcum <- cumsum(contigs2length$Length)
        return(contigs2length)
    }

    derive_probable_taxonomy <- function(contigs2length = NULL){
        N50 <- getN50(contigs2length)
        LKT_table <- table(subset(contigs2length, Length >= N50)[]$LKT)
        LKT_table <- LKT_table[order(LKT_table, decreasing = TRUE)]
        MPT <- names(LKT_table)[1]

        return(MPT)
    }

    find16SrRNA <- function(contigs2length = NULL, opt = NULL){
        featuredata <- opt$featuredata
        #Find how many 16S rRNAs which have been classified can be found in the contigs were looking at.
        num16S <- length(which((subset(featuredata, Contig %in% contigs2length$Contig)[]$SixteenSid) != "none"))

        return(num16S)
    }

    get_assembly_stats_by_bin <- function(opt = NULL, bininterest = NULL){
        contigs2length <- getcontigs2length(opt=opt, bininterest = bininterest)
        assemblystats <- data.frame(MetaBATbin = bininterest, MAG_MPT = derive_probable_taxonomy(contigs2length), NumContigs = nrow(contigs2length), ContigSum = sum(contigs2length$Length), LargestContigSize = max(contigs2length$Length), N50 = getN50(contigs2length), L50 = getL50(contigs2length), N90 = getN90(contigs2length), Num16S = find16SrRNA(contigs2length = contigs2length, opt = opt))

        taxinteresttype <- gsub("^LKT__", "", assemblystats$MAG_MPT)
        expgensize <- as.integer(JAMSMedian_Genome_Sizes[taxinteresttype, "Genome_Size_Median"])
        if ((length(expgensize) < 1) | (is.na(expgensize)) | (taxinteresttype %in% c("Unclassified", "Missing", "Unclassified__Unclassified"))){
            expgensize <- as.integer(JAMSMedian_Genome_Sizes["k__Bacteria", "Genome_Size_Median"])
        }
        probnumgen <- round((sum(contigs2length$Length) / expgensize), 2)
        assemblystats$ProbNumGenomes <- probnumgen

        return(assemblystats)
    }

    get_assembly_stats_for_MAGs <- function(opt = NULL){
        contigsdata <- opt$contigsdata
        #Get a vector of bin names
        named_bins <- unique(opt$contigsdata$MetaBATbin)
        named_bins <- named_bins[which(named_bins != "none")]

        assemblystats_MAGs <- lapply(named_bins, function(t) { get_assembly_stats_by_bin(opt = opt, bininterest = t) })
        assemblystats_MAGs <- plyr::ldply(assemblystats_MAGs, rbind)
        assemblystats_MAGs <- assemblystats_MAGs[order(assemblystats_MAGs$ProbNumGenomes, decreasing = TRUE), ]

        return(assemblystats_MAGs)
    }

    #Load expected genome size data
    JAMSMedian_Genome_Sizes_file <- file.path(opt$workingkrakendb, "JAMSMedian_Genome_Sizes.rda")
    if (file.exists(JAMSMedian_Genome_Sizes_file)){
        load(JAMSMedian_Genome_Sizes_file)
    } else {
        if (file.exists(file.path(opt$workingkrakendb, "JAMS_Median_Genome_Sizes.tsv"))){
            JAMSMedian_Genome_Sizes <- fread(file.path(opt$workingkrakendb, "JAMS_Median_Genome_Sizes.tsv"), data.table = FALSE, header = TRUE, stringsAsFactors = FALSE)
        } else {
            #Fall back on generic taxonomy table and warn user
            flog.info("JAMS Median Genome Sizes table not found. Falling back on generic JAMS Median Genome Sizes table.")
            data(JAMSMedian_Genome_Sizes)
        }
    }
    rownames(JAMSMedian_Genome_Sizes) <- JAMSMedian_Genome_Sizes$Taxon

    opt$assemblystats_MAGs <- get_assembly_stats_for_MAGs(opt=opt)

    return(opt)
}
