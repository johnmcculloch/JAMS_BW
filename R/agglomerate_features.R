#' agglomerate_features(mgseqobj=NULL, glomby=NULL)
#'
#' Agglomerates features in a JAMS mgSeq object safely
#' @export

agglomerate_features <- function(mgseqobj = NULL, glomby = NULL){

    #Find out what kind of an object it is
    analysis <- attr(mgseqobj, "analysis")
    glommgseqobj <- aggTax(mgseqobj, lvl = glomby, out = 'MRexperiment', norm = FALSE)

    if (analysis != "LKT"){
        cts <- MRcounts(glommgseqobj)
        ftt <- fData(glommgseqobj)
        ftt$Accession <- rownames(ftt)
        ftt$Description <- rep(glomby, nrow(ftt))
        ftt <- ftt[, c("Accession", "Description", glomby)]
        ftt <- ftt[rownames(cts), ]
        fData(glommgseqobj) <- ftt
        attr(glommgseqobj, "analysis") <- paste(analysis, glomby, sep = "_")
    }

    return(glommgseqobj)
}
