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
