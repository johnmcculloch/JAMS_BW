#' get_feature_stats
#'
#' JAMSalpha function
#' @export

get_feature_stats <- function(opt = NULL, ucoindex = "freq", doContigs = TRUE, doFeatures = FALSE){

    flog.info("Calculating Tetranucleotide frequencies for all features")

    #Cap threads to 32 because of memory issues
    appropriatenumcores <-  min(max(1, (opt$threads - 2)), 32)

    get_TNF <- function(sequence = NULL){

        TNFdf <- as.data.frame(Biostrings::oligonucleotideFrequency(Biostrings::DNAString(paste(as.character(unlist(sequence)), collapse = "")), 4), stringsAsFactors = FALSE)
        colnames(TNFdf) <- "Freq"
        TNFdf$Tetranucleotide <- rownames(TNFdf)
        TNFdf <- TNFdf[ , c("Tetranucleotide", "Freq")]

        return(TNFdf)
    }

    if (doContigs){ 
        flog.info("Calculating Tetranucleotide frequencies for all contigs")
        TNF_list_contigs <- mclapply(1:length(opt$NHcontigs_sequence), function (x) { get_TNF(opt$NHcontigs_sequence[[x]]) }, mc.cores = appropriatenumcores)
        names(TNF_list_contigs) <- names(opt$NHcontigs_sequence)
        TNF_contigs <- plyr::ldply(TNF_list_contigs, rbind)
        colnames(TNF_contigs)[1] <- "Contig"
        TNF_contigs <- as.data.table(TNF_contigs)
        TNF_contigs <- pivot_wider(data = TNF_contigs, names_from = "Contig", values_from = "Freq", values_fill = 0)
        Tetranucleotides <- TNF_contigs$Tetranucleotide
        TNF_contigs$Tetranucleotide <- NULL
        contig_names <- colnames(TNF_contigs)
        TNF_contigs <- transpose(TNF_contigs)
        TNF_contigs <- as.data.frame(TNF_contigs)
        colnames(TNF_contigs) <- Tetranucleotides
        rownames(TNF_contigs) <- contig_names
        opt$TNF_contigs <- as.data.frame(TNF_contigs)
    }

    if (doFeatures){ 
        TNF_list_features <- mclapply(1:length(opt$genes), function (x) { get_TNF(opt$genes[[x]]) }, mc.cores = appropriatenumcores)
        names(TNF_list_features) <- names(opt$genes)
        TNF_features <- plyr::ldply(TNF_list_features, rbind)
        colnames(TNF_features)[1] <- "Feature"
        TNF_features <- as.data.table(TNF_features)
        TNF_features <- tidyr::pivot_wider(data = TNF_features, names_from = "Feature", values_from = "Freq", values_fill = 0)
        Tetranucleotides <- TNF_features$Tetranucleotide
        TNF_features$Tetranucleotide <- NULL
        feature_names <- colnames(TNF_features)
        TNF_features <- transpose(TNF_features)
        TNF_features <- as.data.frame(TNF_features)
        colnames(TNF_features) <- Tetranucleotides
        rownames(TNF_features) <- feature_names
        opt$TNF_features <- as.data.frame(TNF_features)

        flog.info("Calculating codon usage frequencies for all features")
        ucolist <- list()
        ucolist <- lapply(1:length(opt$genes), function (x) { uco(opt$genes[[x]], index = ucoindex, NA.rscu = 0) })
        names(ucolist) <- names(opt$genes)

        ucobias <- plyr::ldply(ucolist, rbind)
        colnames(ucobias)[which(colnames(ucobias) == ".id")] <- "Feature"
        rownames(ucobias) <- ucobias$Feature

        opt$ucobias <- ucobias
    }

    return(opt)
}
