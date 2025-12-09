#' retrieve_features_by_taxa(FuncExpObj = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0)
#'
#' Returns a long form data frame of stratification by taxa of the relative abundance or number of bases wanted of functional features in wanted samples, given allfeaturesbytaxa_matrix and allfeaturesbytaxa_index metadata present in a JAMS SummarizedExperiment functional object.
#' @export

retrieve_features_by_taxa <- function(FuncExpObj = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0){

    if (assay_for_matrix == "GeneCounts"){
        sparsematrix_name <- "allfeaturesbytaxa_GeneCounts_matrix"
        sparsematrix_index_name <- "allfeaturesbytaxa_GeneCounts_index"
    } else {
        sparsematrix_name <- "allfeaturesbytaxa_matrix"
        sparsematrix_index_name <- "allfeaturesbytaxa_index"
    }

    #Determine first if FuncExpObj is a JAMS version > 2.1 object
    if ("SparseIndex" %in% names(assays(FuncExpObj))){
        JAMS2_object <- TRUE
        index_matrix <- as.data.frame(assays(FuncExpObj)[["SparseIndex"]])
        #Subset to only wanted features and samples
        index_matrix <- index_matrix[wantedfeatures, wantedsamples, drop = FALSE]
        rowsinterestdf <- as.data.frame(index_matrix) %>% rownames_to_column("Accession") %>% pivot_longer(cols = -Accession, names_to = "Sample", values_to = "RowNumber")
        rowsinterestdf <- as.data.frame(rowsinterestdf)

    } else {
        JAMS2_object <- FALSE
        allfeaturesbytaxa_index <- metadata(FuncExpObj)[[sparsematrix_index_name]]
        #Get appropriate rows
        rowsinterestdf <- subset(allfeaturesbytaxa_index, Sample %in% wantedsamples)
        rowsinterestdf <- subset(rowsinterestdf, Accession %in% wantedfeatures)
    }

    #Ensure correct class
    rowsinterestdf$RowNumber <- as.integer(rowsinterestdf$RowNumber)
    allfeaturesbytaxa_interest <- as.matrix(metadata(FuncExpObj)[[sparsematrix_name]][rowsinterestdf$RowNumber, , drop = FALSE])
    #Prune empty taxa, i.e. taxa which have 0 counts for all of the features of interest
    allfeaturesbytaxa_interest <- as.data.frame(allfeaturesbytaxa_interest[, which(colSums(allfeaturesbytaxa_interest) != 0)])
    allfeaturesbytaxa_interest$RowNumber <- as.integer(rownames(allfeaturesbytaxa_interest))
    allfeaturesbytaxa_interest <- left_join(allfeaturesbytaxa_interest, rowsinterestdf, by = "RowNumber")
    allfeaturesbytaxa_interest$RowNumber <- NULL
    allfeaturesbytaxa_interest <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", (sort(colnames(allfeaturesbytaxa_interest)[which(!colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession"))])))]

    if (asPPM){
        taxsplit <- allfeaturesbytaxa_interest
        #Transform to PPM
        if (PPM_normalize_to_bases_sequenced == TRUE){
            totbases <- "TotalBasesSequenced"
        } else {
            totbases <- "TotalBasesSequencedinAnalysis"
        }
        numbases2sampl <- as.data.frame(t(metadata(FuncExpObj)[[totbases]]))
        numbases2sampl$Sample <- rownames(numbases2sampl)
        taxsplit <- left_join(taxsplit, numbases2sampl, by = "Sample")
        LKTcolumns <- colnames(taxsplit)[!(colnames(taxsplit) %in% c("Sample", "Accession", "NumBases"))]

        #Transform to PPM
        for (colm in LKTcolumns){
            taxsplit[ , colm] <- round(((taxsplit[ , colm] / taxsplit$NumBases) * 1000000), 0)
        }

        #Denoise
        LKTsMaxima <-sapply(LKTcolumns, function(x) {max(taxsplit[ , x])})
        LKTsToKeep <- names(which(LKTsMaxima > PPMthreshold))
        sample2metadata <- as.data.frame(colData(FuncExpObj))
        nonSamplecolms <- colnames(sample2metadata)[colnames(sample2metadata) != "Sample"]
        taxsplit <- left_join(taxsplit, sample2metadata, by = "Sample")
        taxsplit <- taxsplit[ , c("Sample", "Accession", "NumBases", nonSamplecolms, LKTsToKeep)]
        taxsplit$NumBases <- NULL

        return(taxsplit)

    } else {

        return(allfeaturesbytaxa_interest)

    }
}
