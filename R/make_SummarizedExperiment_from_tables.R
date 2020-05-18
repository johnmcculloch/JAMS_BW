#' make_SummarizedExperiment_from_tables(pheno = NULL, pheno_fn= NULL, counttable = NULL, counttable_fn = NULL, featuretable = NULL, featuretable_fn = NULL, onlysamples = NULL, analysisname = NULL, restricttoLKTs = NULL)
#'
#' Safely makes a JAMS-compatible SummarizedExperiment object from a phenotable, a counttable and a featuretable.
#' @export

make_SummarizedExperiment_from_tables <- function(pheno = NULL, pheno_fn= NULL, counttable = NULL, counttable_fn = NULL, featuretable = NULL, featuretable_fn = NULL, onlysamples = NULL, analysisname = NULL, restricttoLKTs = NULL){

    md_list <- list(pheno, counttable, featuretable)
    fn_list <- list(pheno_fn, counttable_fn, featuretable_fn)

    if (any(sapply(1:length(md_list), function(x) { is.null(md_list[[x]])}))){
        missing_element <- which(sapply(1:length(md_list), function(x) { is.null(md_list[[x]])}) == TRUE)
        for (miss in missing_element){
            md_list[[miss]] <- read.table(file = fn_list[[miss]], sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
        }
    }

    #Fix phenotable
    names(md_list) <- c("pheno", "cts", "ft")
    rownames(md_list$pheno) <- md_list$pheno$Sample
    if (!is.null(onlysamples)){
        md_list$pheno <- md_list$pheno[(rownames(rownames(md_list$pheno)) %in% onlysamples), ]
    }

    #Fix countable
    #Test if first column is of feature names
    if (!is.numeric(md_list$cts[ , 1])){
        colnames(md_list$cts)[1] <- "Accession"
        rownames(md_list$cts) <- md_list$cts$Accession
        md_list$cts$Accession <- NULL
    }

    #make numeric
    for (colm in 1:ncol(md_list$cts)){
        md_list$cts[ , colm] <- as.numeric(md_list$cts[ , colm])
    }
    md_list$cts <- md_list$cts[order(rowSums(md_list$cts), decreasing = TRUE), ]
    featureorder <- rownames(md_list$cts)

    print(featureorder[1:100])

    #Fix feat table
    #Assume first column are feature names
    colnames(md_list$ft)[1] <- "Accession"
    rownames(md_list$ft) <- md_list$ft$Accession
    #md_list$ft$Accession <- NULL
    md_list$ft <- md_list$ft[(rownames(md_list$ft) %in% featureorder), ]
    sampleorder <- rownames(md_list$pheno)
    md_list$ft <- as.matrix(md_list$ft)
    md_list$ft <- md_list$ft[featureorder, ]
    md_list$cts <- as.matrix(md_list$cts)
    md_list$cts <- md_list$cts[, sampleorder]

    #Register the total number of bases (or counts) sequenced for each sample
    TotBasesSamples <- colSums(md_list$cts)
    TotalBasesSequenced <- t(as.matrix(TotBasesSamples))
    rownames(TotalBasesSequenced) <- "NumBases"

    assays <- list(md_list$cts)
    names(assays) <- "BaseCounts"

    SEobj <- SummarizedExperiment(assays = assays, rowData = md_list$ft, colData = md_list$pheno)
    metadata(SEobj)$TotalBasesSequenced <- TotalBasesSequenced
    metadata(SEobj)$TotalBasesSequencedinAnalysis <- TotalBasesSequenced #Have to assume this as info is missing.
    metadata(SEobj)$analysis <- analysisname

    return(SEobj)
}
