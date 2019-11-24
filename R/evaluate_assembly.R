#' evaluate_assembly(opt=NULL)
#'
#' Calculates the genome completeness and assembly stats for each taxonomic level in a contigsdata dataframe.
#' @export

evaluate_assembly<-function(opt=NULL){

    flog.info("Evaluating completeness of genomes.")

    #Load expected genome size data
    JAMSMedian_Genome_Sizes_file <- file.path(opt$workingkrakendb, "JAMS_Median_Genome_Sizes.tsv")
    if (file.exists(JAMSMedian_Genome_Sizes_file)){
        JAMSMedian_Genome_Sizes <- read.table(file = JAMSMedian_Genome_Sizes_file, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, quote = NULL, header = TRUE, colClasses = c("character", "numeric", "numeric", "numeric", "character"))
    } else {
        #Fall back on generic taxonomy table and warn user
        flog.info("JAMS Median Genome Sizes table not found. Falling back on generic JAMS Median Genome Sizes table.")
        data(JAMSMedian_Genome_Sizes)
    }
    rownames(JAMSMedian_Genome_Sizes) <- JAMSMedian_Genome_Sizes$Taxon

    taxlvls <- c("LKT", "IS1", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain")
    assemblystats_alllevels <- mclapply(taxlvls, function (x) { get_assembly_stats_by_taxlevel(opt = opt, taxlevel = x) }, mc.cores = max(1, min(length(taxlvls), (opt$threads - 2))))

    assemblystats_alllevels <- plyr::ldply(assemblystats_alllevels, rbind)

    opt$assemblystats <- assemblystats_alllevels

    return(opt)
}
