#' get_feature_stats
#'
#' JAMSalpha function
#' @export

get_feature_stats <- function(opt = NULL, ucoindex = "freq"){

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
