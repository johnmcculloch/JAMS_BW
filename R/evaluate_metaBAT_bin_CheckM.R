#' evaluate_metaBAT_bin_CheckM(opt=NULL)
#'
#' Calculates the genome completeness and assembly stats for metaBAT bins using CheckM (must be installed).
#' @export

evaluate_metaBAT_bin_CheckM <- function(opt = NULL, verbose = FALSE){

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
        if (verbose){
            flog.info(paste("Writing", MAG))
        }
        write.fasta(sequences = opt$NHcontigs_sequence[Wanted_contigs], names = Wanted_contigs, nbchar = 80, file.out = currtaxfn)
    }

    #Set the CheckM taxonomy database path
    system2("checkm", args = c("data", "setRoot", opt$CheckMdb))
    #Run CheckM
    system2("checkm", args = c("lineage_wf", "-t", opt$threads, opt$tempCheckMbinfolder, opt$tempCheckMresultsfolder), stdout = FALSE, stderr = FALSE)
    system2("checkm", args = c("qa", "-f", file.path(opt$tempCheckMresultsfolder, "CheckM_table.tsv"), "--tab_table", "-t", opt$threads, file.path(opt$tempCheckMresultsfolder, "lineage.ms"), opt$tempCheckMresultsfolder), stdout = FALSE, stderr = FALSE)

    #Load results from CheckM into opt
    CheckMout <- fread(file.path(opt$tempCheckMresultsfolder, "CheckM_table.tsv"), header = TRUE, data.table = FALSE, sep = "\t")
    colnames(CheckMout) <- c("MetaBATbin", "Marker_lineage", "Num_ref_genomes", "Num_ref_markers", "Num_ref_marker_sets", "C0", "C1", "C2", "C3", "C4", "C5_or_more", "Completeness", "Contamination", "Strain heterogeneity")

    opt$CheckMout <- CheckMout[order(CheckMout$Completeness, decreasing = TRUE), ]

    opt$assemblystats_MAGs <- get_assembly_stats_for_MAGs(opt=opt)

    return(opt)
}
