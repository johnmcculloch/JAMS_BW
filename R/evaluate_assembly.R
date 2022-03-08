#' evaluate_assembly(opt=NULL)
#'
#' Calculates the genome completeness and assembly stats for each taxonomic level in a contigsdata dataframe.
#' @export

evaluate_assembly <- function(opt = NULL){

    flog.info("Evaluating completeness of genomes.")

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

    estimate_genome_completeness <- function(contigs2length = NULL){
        taxinteresttype <- gsub("^LKT__", "", unique(contigs2length$Taxon))
        expgensize <- as.integer(JAMSMedian_Genome_Sizes[taxinteresttype, "Genome_Size_Median"])
        if ((length(expgensize) < 1) | (is.na(expgensize)) | (taxinteresttype %in% c("Unclassified", "Missing", "Unclassified__Unclassified"))){
            expgensize <- as.integer(JAMSMedian_Genome_Sizes["k__Bacteria", "Genome_Size_Median"])
        }
        probnumgen <- round((sum(contigs2length$Length) / expgensize), 2)
        return(probnumgen)
    }

    find16SrRNA <- function(contigs2length = NULL, opt = NULL){
        featuredata = opt$featuredata
        #Find how many 16S rRNAs which have been classified can be found in the contigs were looking at.
        num16S <- length(which((subset(featuredata, Contig %in% contigs2length$Contig)[]$SixteenSid) != "none"))

        return(num16S)
    }

    getcontigs2length <- function(opt = NULL, taxlevel = NULL, taxinterest = NULL){
        contigsdata <- opt$contigsdata
        contigs2length <- subset(contigsdata, contigsdata[,taxlevel] == taxinterest)
        contigs2length <- contigs2length[, c(taxlevel, "Contig", "Length")]
        colnames(contigs2length) <- c("Taxon", "Contig", "Length")
        contigs2length <- contigs2length[order(contigs2length$Length, decreasing = TRUE), ]
        contigs2length$Lcum <- cumsum(contigs2length$Length)
        return(contigs2length)
    }

    get_assembly_stats_by_taxon <- function(opt = NULL, taxlevel = NULL, taxinterest = NULL){
        contigs2length <- getcontigs2length(opt=opt, taxlevel=taxlevel, taxinterest = taxinterest)
        assemblystats <- data.frame(TaxLevel = taxlevel, Taxon = taxinterest, NumContigs = nrow(contigs2length), ContigSum = sum(contigs2length$Length), LargestContigSize = max(contigs2length$Length), N50 = getN50(contigs2length), L50 = getL50(contigs2length), N90 = getN90(contigs2length), ProbNumGenomes = estimate_genome_completeness(contigs2length), Num16S = find16SrRNA(contigs2length = contigs2length, opt = opt))
        return(assemblystats)
    }

    get_assembly_stats_by_taxlevel <- function(opt = NULL, taxlevel = NULL){
        contigsdata <- opt$contigsdata
        #Get a list of non-redundant tax entities
        Taxa <- unique(contigsdata[ , taxlevel])
        assemblystats_taxlevel <- lapply(Taxa, function(t) { get_assembly_stats_by_taxon(opt = opt, taxlevel = taxlevel, taxinterest = t) })
        assemblystats_taxlevel <- plyr::ldply(assemblystats_taxlevel, rbind)
        assemblystats_taxlevel <- assemblystats_taxlevel[order(assemblystats_taxlevel$ProbNumGenomes, decreasing = TRUE), ]
        return(assemblystats_taxlevel)
    }

    #Load expected genome size data
    JAMSMedian_Genome_Sizes_file <- file.path(opt$workingkrakendb, "JAMSMedian_Genome_Sizes.rda")
    if (file.exists(JAMSMedian_Genome_Sizes_file)){
        load(JAMSMedian_Genome_Sizes_file)
    } else {
        #make compatioble with old database versions
        if (file.exists(file.path(opt$workingkrakendb, "JAMS_Median_Genome_Sizes.tsv"))){
            JAMSMedian_Genome_Sizes <- fread(file.path(opt$workingkrakendb, "JAMS_Median_Genome_Sizes.tsv"), data.table = FALSE, header = TRUE, stringsAsFactors = FALSE)
        } else {
            #Fall back on generic taxonomy table and warn user
            flog.info("JAMS Median Genome Sizes table not found. Falling back on generic JAMS Median Genome Sizes table.")
            data(JAMSMedian_Genome_Sizes)
        }
    }
    rownames(JAMSMedian_Genome_Sizes) <- JAMSMedian_Genome_Sizes$Taxon

    #export to global env for the benefit of legacy versions
    assign("JAMSMedian_Genome_Sizes", JAMSMedian_Genome_Sizes, .GlobalEnv)

    taxlvls <- c("LKT", "IS1", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain")
    assemblystats_alllevels <- mclapply(taxlvls, function (x) { get_assembly_stats_by_taxlevel(opt = opt, taxlevel = x) }, mc.cores = max(1, min(length(taxlvls), (opt$threads - 2))))

    assemblystats_alllevels <- plyr::ldply(assemblystats_alllevels, rbind)

    opt$assemblystats <- assemblystats_alllevels

    return(opt)
}
