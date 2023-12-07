#' evaluate_metaBAT_bin_CheckM(opt=NULL)
#'
#' Calculates the genome completeness and assembly stats for metaBAT bins using CheckM (must be installed).
#' @export

evaluate_metaBAT_bin_CheckM <- function(opt = NULL){

    flog.info("Evaluating quality of MAGs with CheckM.")

    ######################################
    ##   Declare functions that must    ##
    ## be declared within this funciton ##
    ######################################

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


    #Get a vector of bin names
    opt$tempCheckMbinfolder <- file.path(opt$sampledir, "CheckM_bins")
    dir.create(opt$tempCheckMbinfolder)
    opt$tempCheckMresultsfolder <- file.path(opt$sampledir, "CheckM_results")
    dir.create(opt$tempCheckMresultsfolder)

    named_bins <- unique(opt$contigsdata$MetaBATbin)
    named_bins <- named_bins[which(named_bins != "none")]

    for (MAG in named_bins){
        Wanted_contigs <- NULL
        currtaxfn <- paste(MAG, "fna", sep = ".")
        currtaxfn <- file.path(opt$tempCheckMbinfolder, currtaxfn)
        Wanted_contigs <- subset(opt$contigsdata, MetaBATbin == MAG)[]$Contig
        flog.info(paste("Writing", MAG))
        write.fasta(sequences = opt$NHcontigs_sequence[Wanted_contigs], names = Wanted_contigs, nbchar = 80, file.out = currtaxfn)
    }

    #Run CheckM
    system2("checkm", args = c("lineage_wf", "-t", opt$threads, opt$tempCheckMbinfolder, opt$tempCheckMresultsfolder), stdout = FALSE, stderr = FALSE)
    #Load results from CheckM into opt


    bin_stats_ext <- fread(file.path(opt$tempCheckMresultsfolder, "storage", "bin_stats_ext.tsv"), header = FALSE, data.table = FALSE)
    bin_stats_analyze <- fread(file.path(opt$tempCheckMresultsfolder, "storage", "bin_stats.analyze.tsv"), header = FALSE, data.table = FALSE)


    #Clean up


    opt$assemblystats_MAGs <- get_assembly_stats_for_MAGs(opt=opt)

    return(opt)
}
