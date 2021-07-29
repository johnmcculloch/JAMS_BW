#' get_feature_stats
#'
#' JAMSalpha function
#' @export

get_feature_stats <- function(opt = NULL, ucoindex = "freq"){

    flog.info("Calculating Tetranucleotide frequencies for all features")

    get_TNF <- function(sequence = NULL){

        TNFdf <- as.data.frame(oligonucleotideFrequency(DNAString(paste(as.character(unlist(sequence)), collapse = "")), 4), stringsAsFactors = FALSE)
        colnames(TNFdf) <- "Freq"
        TNFdf$Tetranucleotide <- rownames(TNFdf)
        TNFdf <- TNFdf[ , c("Tetranucleotide", "Freq")]

        return(TNFdf)
    }

    TNF_list_features <- lapply(1:length(opt$genes), function (x) { get_TNF(opt$genes[[x]]) } )
    names(TNF_list_features) <- names(opt$genes)

    TNF_features <- plyr::ldply(TNF_list, rbind)
    colnames(TNF_features)[1] <- "Feature"

    TNF_features <- as.data.frame(pivot_wider(data = TNF_features, names_from = "Feature", values_from = "Freq", values_fill = 0))
    rownames(TNF_features) <- TNF_features$Tetranucleotide
    TNF_features$Tetranucleotide <- NULL
    TNF_features <- t(TNF_features)
    opt$TNF_features <- TNF_features

    flog.info("Calculating Tetranucleotide frequencies for all contigs")
    TNF_list_contigs <- lapply(1:length(opt$NHcontigs_sequence), function (x) { get_TNF(opt$NHcontigs_sequence[[x]]) } )
    names(TNF_list_contigs) <- names(opt$NHcontigs_sequence)

    TNF_contigs <- plyr::ldply(TNF_list_contigs, rbind)
    colnames(TNF_contigs)[1] <- "Contig"

    TNF_contigs <- as.data.frame(pivot_wider(data = TNF_contigs, names_from = "Contig", values_from = "Freq", values_fill = 0))
    rownames(TNF_contigs) <- TNF_contigs$Tetranucleotide
    TNF_contigs$Tetranucleotide <- NULL
    TNF_contigs <- t(TNF_contigs)
    opt$TNF_contigs <- TNF_contigs

    flog.info("Calculating codon usage frequencies for all features")
    ucolist <- list()
    ucolist <- lapply(1:length(opt$genes), function (x) { uco(opt$genes[[x]], index = ucoindex, NA.rscu = 0) })
    names(ucolist) <- names(opt$genes)

    ucobias <- plyr::ldply(ucolist, rbind)
    colnames(ucobias)[which(colnames(ucobias) == ".id")] <- "Feature"
    rownames(ucobias) <- ucobias$Feature

    opt$ucobias <- ucobias

    return(opt)
}
